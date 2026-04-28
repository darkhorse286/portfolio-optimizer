// SPDX-License-Identifier: MIT

#include "quantum/benchmark_runner.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <nlohmann/json.hpp>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "data/data_loader.hpp"

namespace portfolio
{
    namespace quantum
    {

        BenchmarkRunner::BenchmarkRunner(const portfolio::MarketData& data,
                                         const portfolio::backtest::BacktestParams& backtest_params)
            : data_(data), backtest_params_(backtest_params)
        {
        }

        // Hardcoded sector mapping for the 10-asset universe
        static std::unordered_map<std::string, std::string> get_sector_mapping()
        {
            return {
                {"AAPL", "Technology"},
                {"MSFT", "Technology"},
                {"GOOGL", "Technology"},
                {"JPM", "Financials"},
                {"BAC", "Financials"},
                {"JNJ", "Healthcare"},
                {"PFE", "Healthcare"},
                {"XOM", "Energy"},
                {"CVX", "Energy"},
                {"WMT", "Consumer"}
            };
        }

        // Compute sector weights from portfolio weights
        static std::unordered_map<std::string, double> compute_sector_weights(
            const std::vector<std::string>& universe,
            const Eigen::VectorXd& weights)
        {
            auto sector_map = get_sector_mapping();
            std::unordered_map<std::string, double> sector_weights;
            for (size_t i = 0; i < universe.size(); ++i)
            {
                std::string sector = sector_map[universe[i]];
                sector_weights[sector] += weights[i];
            }
            return sector_weights;
        }

        void BenchmarkRunner::add_classical_solver(
            const std::string& name,
            std::shared_ptr<portfolio::optimizer::MeanVarianceOptimizer> solver)
        {
            RegisteredSolver rs;
            rs.name = name;
            rs.solver_type = "classical";
            rs.execution_backend = "osqp";  // assuming
            rs.classical = solver;
            solvers_.push_back(rs);
        }

        void BenchmarkRunner::add_quantum_solver(
            const std::string& name,
            const std::string& solver_type,
            std::shared_ptr<portfolio::quantum::QuantumOptimizer> solver)
        {
            RegisteredSolver rs;
            rs.name = name;
            rs.solver_type = solver_type;
            rs.execution_backend = solver->execution_backend();
            rs.quantum = solver;
            solvers_.push_back(rs);
        }

        BenchmarkResult BenchmarkRunner::run(int num_runs)
        {
            BenchmarkResult result;

            // Generate run_id
            auto now = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::tm tm = *std::localtime(&time_t);
            std::ostringstream oss;
            oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
            result.run_id = oss.str();

            result.problem_sizes.push_back(static_cast<int>(data_.num_assets()));

            // Run each solver
            for (const auto& s : solvers_)
            {
                result.results.push_back(run_solver(s, num_runs));
            }

            // Find classical baseline (first classical)
            double classical_sharpe = 1.0;
            for (const auto& r : result.results)
            {
                if (r.solver_type == "classical")
                {
                    classical_sharpe = r.mean_sharpe;
                    break;
                }
            }

            // Compute solution_quality_vs_classical
            for (auto& r : result.results)
            {
                if (classical_sharpe == 0.0)
                    r.solution_quality_vs_classical = 0.0;
                else
                    r.solution_quality_vs_classical = r.mean_sharpe / classical_sharpe;
            }

            return result;
        }

        SolverResult BenchmarkRunner::run_solver(const RegisteredSolver& s, int num_runs) const
        {
            SolverResult result;
            result.solver_name = s.name;
            result.solver_type = s.solver_type;
            result.execution_backend = s.execution_backend;
            result.problem_size = static_cast<int>(data_.num_assets());
            result.num_runs = num_runs;
            result.successful_runs = 0;

            std::vector<double> sharpes;
            std::vector<double> returns;
            std::vector<double> volatilities;
            std::vector<double> solve_times;
            std::vector<double> circuit_times;
            std::vector<std::vector<double>> all_nav_series;
            std::vector<portfolio::backtest::TradeSummary> trade_summaries;
            std::vector<Eigen::VectorXd> all_weights;

            for (int i = 0; i < num_runs; ++i)
            {
                try
                {
                    portfolio::backtest::BacktestEngine engine(backtest_params_);
                    portfolio::backtest::BacktestResult bt_result;

                    auto t_start = std::chrono::high_resolution_clock::now();
                    if (s.classical)
                    {
                        bt_result = engine.run(data_, *s.classical);
                    }
                    else if (s.quantum)
                    {
                        bt_result = engine.run(data_, *s.quantum);
                    }
                    auto t_end = std::chrono::high_resolution_clock::now();
                    double wall_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

                    if (bt_result.success)
                    {
                        auto analytics = bt_result.compute_analytics(backtest_params_.risk_free_rate, 252);

                        sharpes.push_back(analytics.sharpe_ratio());
                        returns.push_back(analytics.total_return());
                        volatilities.push_back(analytics.annualized_volatility());

                        if (s.quantum)
                        {
                            // quantum solver reports per-rebalance solve time internally
                            solve_times.push_back(s.quantum->solve_time_ms());
                            circuit_times.push_back(s.quantum->circuit_execution_us());
                        }
                        else
                        {
                            solve_times.push_back(wall_ms);
                            circuit_times.push_back(-1.0);
                        }

                        result.successful_runs++;
                        all_nav_series.push_back(bt_result.nav_series);
                        trade_summaries.push_back(bt_result.trade_summary);

                        // Compute average weights across snapshots
                        Eigen::VectorXd avg_weights = Eigen::VectorXd::Zero(data_.num_assets());
                        for (const auto& snap : bt_result.snapshots)
                        {
                            avg_weights += snap.weights;
                        }
                        if (!bt_result.snapshots.empty())
                        {
                            avg_weights /= bt_result.snapshots.size();
                        }
                        all_weights.push_back(avg_weights);

                        // Dates are identical across runs; keep from first successful run
                        if (result.date_series.empty())
                            result.date_series = bt_result.dates;
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: Run " << i << " for " << s.name << " failed: " << e.what() << std::endl;
                }
            }

            // Element-wise mean of nav_series so the equity curve matches the averaged metrics
            if (!all_nav_series.empty())
            {
                size_t len = all_nav_series[0].size();
                result.nav_series.assign(len, 0.0);
                for (const auto& nav : all_nav_series)
                {
                    size_t usable = std::min(len, nav.size());
                    for (size_t k = 0; k < usable; ++k)
                        result.nav_series[k] += nav[k];
                }
                double n = static_cast<double>(all_nav_series.size());
                for (double& v : result.nav_series)
                    v /= n;
            }

            if (result.successful_runs > 0)
            {
                auto mean_val = [](const std::vector<double>& v) {
                    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
                };

                // Solve time and circuit time: mean across runs (unchanged)
                result.mean_solve_time_ms = mean_val(solve_times);
                result.mean_circuit_execution_us = mean_val(circuit_times);

                // Sharpe, return, volatility: recompute from the averaged NAV series
                // so these metrics are self-consistent with the equity curve.
                const auto& nav = result.nav_series;
                if (nav.size() >= 2 && nav.front() > 0.0)
                {
                    std::vector<double> daily_rets;
                    daily_rets.reserve(nav.size() - 1);
                    for (size_t k = 1; k < nav.size(); ++k)
                        if (nav[k - 1] > 0.0)
                            daily_rets.push_back((nav[k] - nav[k - 1]) / nav[k - 1]);

                    double total_ret = nav.back() / nav.front() - 1.0;
                    double ann_ret   = std::pow(1.0 + total_ret, 252.0 / nav.size()) - 1.0;
                    double mean_r    = std::accumulate(daily_rets.begin(), daily_rets.end(), 0.0)
                                       / daily_rets.size();
                    double sq_sum    = std::inner_product(daily_rets.begin(), daily_rets.end(),
                                                          daily_rets.begin(), 0.0);
                    double ann_vol   = std::sqrt(sq_sum / daily_rets.size() - mean_r * mean_r)
                                       * std::sqrt(252.0);

                    result.mean_portfolio_return     = total_ret;
                    result.mean_portfolio_volatility = ann_vol;
                    result.mean_sharpe               = (ann_vol > 0.0)
                        ? (ann_ret - backtest_params_.risk_free_rate) / ann_vol
                        : 0.0;
                    result.std_sharpe = 0.0; // std dev across runs kept as 0; mean is from averaged curve
                }
                else
                {
                    result.mean_sharpe               = mean_val(sharpes);
                    result.std_sharpe                = 0.0;
                    result.mean_portfolio_return     = mean_val(returns);
                    result.mean_portfolio_volatility = mean_val(volatilities);
                }
            }
            else
            {
                result.mean_sharpe = 0.0;
                result.std_sharpe = 0.0;
                result.mean_solve_time_ms = 0.0;
                result.mean_circuit_execution_us = -1.0;
                result.mean_portfolio_return = 0.0;
                result.mean_portfolio_volatility = 0.0;
            }

            // Aggregate trade summary across successful runs
            if (!trade_summaries.empty())
            {
                portfolio::backtest::TradeSummary& agg = result.trade_summary;
                for (const auto& ts : trade_summaries)
                {
                    agg.total_trades += ts.total_trades;
                    agg.buy_trades += ts.buy_trades;
                    agg.sell_trades += ts.sell_trades;
                    agg.total_notional += ts.total_notional;
                    agg.total_costs += ts.total_costs;
                    agg.rebalance_count += ts.rebalance_count;
                    agg.turnover += ts.turnover;
                }
                int n = trade_summaries.size();
                agg.avg_cost_per_trade = agg.total_costs / std::max(1, agg.total_trades);
                agg.total_trades /= n;
                agg.buy_trades /= n;
                agg.sell_trades /= n;
                agg.total_notional /= n;
                agg.total_costs /= n;
                agg.rebalance_count /= n;
                agg.turnover /= n;
            }

            // Compute average weights and sector weights
            if (!all_weights.empty())
            {
                Eigen::VectorXd avg_weights = Eigen::VectorXd::Zero(data_.num_assets());
                for (const auto& w : all_weights)
                {
                    avg_weights += w;
                }
                avg_weights /= all_weights.size();

                result.sector_weights = compute_sector_weights(data_.get_tickers(), avg_weights);

                // Store per-asset weights and tickers for per-solver attribution in the viz
                result.tickers = data_.get_tickers();
                result.avg_weights.resize(avg_weights.size());
                for (int i = 0; i < avg_weights.size(); ++i)
                    result.avg_weights[i] = avg_weights[i];

                // For benchmark, use equal weight across sectors (simplified)
                // In a real implementation, this would be the benchmark's sector weights
                std::unordered_map<std::string, double> benchmark_weights;
                auto sector_map = get_sector_mapping();
                std::unordered_map<std::string, int> sector_counts;
                for (const auto& ticker : data_.get_tickers())
                {
                    std::string sector = sector_map[ticker];
                    sector_counts[sector]++;
                }
                for (const auto& [sector, count] : sector_counts)
                {
                    benchmark_weights[sector] = 1.0 / sector_counts.size();  // equal weight
                }
                result.benchmark_sector_weights = benchmark_weights;
            }

            return result;
        }

        std::string BenchmarkResult::to_json() const
        {
            auto compute_max_drawdown = [](const std::vector<double>& nav) -> double {
                if (nav.empty()) return 0.0;
                double peak = nav.front();
                double max_dd = 0.0;
                for (double v : nav) {
                    if (v > peak) peak = v;
                    double dd = (peak > 0.0) ? (peak - v) / peak : 0.0;
                    if (dd > max_dd) max_dd = dd;
                }
                return -max_dd; // negative convention: e.g. -0.225
            };

            nlohmann::json j;
            j["schema_version"] = "1.0";
            j["run_id"] = run_id;
            j["problem_sizes"] = problem_sizes;
            j["runs"] = nlohmann::json::array();

            for (const auto& res : results)
            {
                nlohmann::json run;
                run["solver_name"] = res.solver_name;
                run["solver_type"] = res.solver_type;
                run["execution_backend"] = res.execution_backend;
                run["problem_size"] = res.problem_size;
                run["performance"]["sharpe_ratio"] = res.mean_sharpe;
                run["performance"]["total_return"] = res.mean_portfolio_return;
                run["performance"]["annualized_return"] = res.mean_portfolio_return;
                run["performance"]["annualized_volatility"] = res.mean_portfolio_volatility;
                run["performance"]["max_drawdown"] = compute_max_drawdown(res.nav_series);
                run["solve_time_ms"] = res.mean_solve_time_ms;
                run["circuit_execution_us"] = res.mean_circuit_execution_us;
                run["solution_quality_vs_classical"] = res.solution_quality_vs_classical;

                // nav_series with dates for the viz — cap at 1000 entries
                nlohmann::json nav_arr = nlohmann::json::array();
                {
                    const auto& nav = res.nav_series;
                    size_t total = nav.size();
                    constexpr size_t kMaxPoints = 1000;

                    // Build evenly-spaced index list (always include first and last)
                    std::vector<size_t> indices;
                    if (total <= kMaxPoints)
                    {
                        indices.resize(total);
                        std::iota(indices.begin(), indices.end(), 0);
                    }
                    else
                    {
                        indices.reserve(kMaxPoints);
                        for (size_t k = 0; k < kMaxPoints; ++k)
                            indices.push_back(static_cast<size_t>(
                                std::round(static_cast<double>(k) * (total - 1) / (kMaxPoints - 1))));
                    }

                    double peak = nav.empty() ? 1.0 : nav.front();
                    for (size_t idx : indices)
                    {
                        double v = nav[idx];
                        if (v > peak) peak = v;
                        double dd = (peak > 0.0) ? -(peak - v) / peak : 0.0;
                        double ret = (idx == 0 || nav[idx - 1] <= 0.0)
                                     ? 0.0
                                     : (v - nav[idx - 1]) / nav[idx - 1];

                        nlohmann::json point;
                        point["nav"] = v;
                        point["return"] = ret;
                        point["drawdown"] = dd;
                        if (idx < res.date_series.size())
                            point["date"] = res.date_series[idx];
                        else
                            point["date"] = static_cast<int>(idx);
                        nav_arr.push_back(point);
                    }
                }
                run["nav_series"] = nav_arr;

                // Trade summary
                const auto& ts = res.trade_summary;
                run["trade_summary"]["rebalance_count"] = ts.rebalance_count;
                run["trade_summary"]["turnover"] = ts.turnover;
                run["trade_summary"]["avg_cost_per_trade"] = ts.avg_cost_per_trade;
                run["trade_summary"]["total_costs"] = ts.total_costs;

                // Per-asset weights and tickers for Brinson-Fachler attribution
                run["tickers"] = res.tickers;
                nlohmann::json weights_arr = nlohmann::json::array();
                for (int i = 0; i < static_cast<int>(res.avg_weights.size()); ++i)
                    weights_arr.push_back(res.avg_weights[i]);
                run["avg_weights"] = weights_arr;

                // Sector weights
                run["sector_weights"] = res.sector_weights;
                run["benchmark_sector_weights"] = res.benchmark_sector_weights;

                j["runs"].push_back(run);
            }

            return j.dump(2);
        }

        void BenchmarkResult::print_summary() const
        {
            std::cout << "============================================================" << std::endl;
            std::cout << "Benchmark Results" << std::endl;
            std::cout << "============================================================" << std::endl;
            std::cout << "Solver                    | Size | Sharpe | Return | vs MV" << std::endl;
            std::cout << "--------------------------|------|--------|--------|-------" << std::endl;

            for (const auto& res : results)
            {
                std::cout << std::left << std::setw(25) << res.solver_name << " | "
                          << std::right << std::setw(4) << res.problem_size << " | "
                          << std::fixed << std::setprecision(3) << std::setw(6) << res.mean_sharpe << " | "
                          << std::fixed << std::setprecision(1) << std::setw(6) << (res.mean_portfolio_return * 100) << "% | "
                          << std::fixed << std::setprecision(2) << std::setw(5) << res.solution_quality_vs_classical << "x"
                          << std::endl;
            }
            std::cout << "============================================================" << std::endl;
        }

        void BenchmarkResult::export_comparison_json(const std::string& path) const
        {
            std::ofstream out(path);
            if (!out)
                throw std::runtime_error("Unable to open file for writing: " + path);
            out << to_json();
        }

    } // namespace quantum
} // namespace portfolio