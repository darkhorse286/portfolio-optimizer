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

            for (int i = 0; i < num_runs; ++i)
            {
                try
                {
                    portfolio::backtest::BacktestEngine engine(backtest_params_);
                    portfolio::backtest::BacktestResult bt_result;

                    if (s.classical)
                    {
                        bt_result = engine.run(data_, *s.classical);
                    }
                    else if (s.quantum)
                    {
                        bt_result = engine.run(data_, *s.quantum);
                    }

                    if (bt_result.success)
                    {
                        auto analytics = bt_result.compute_analytics(backtest_params_.risk_free_rate, 252);

                        sharpes.push_back(analytics.sharpe_ratio());
                        returns.push_back(analytics.total_return());
                        volatilities.push_back(analytics.annualized_volatility());

                        if (s.quantum)
                        {
                            solve_times.push_back(s.quantum->solve_time_ms());
                            circuit_times.push_back(s.quantum->circuit_execution_us());
                        }
                        else
                        {
                            solve_times.push_back(0.0);  // or measure?
                            circuit_times.push_back(-1.0);
                        }

                        result.successful_runs++;

                        // Store nav_series and dates from the last successful run
                        result.nav_series = bt_result.nav_series;
                        result.date_series = bt_result.dates;
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: Run " << i << " for " << s.name << " failed: " << e.what() << std::endl;
                }
            }

            if (result.successful_runs > 0)
            {
                auto mean_std = [](const std::vector<double>& v) {
                    double sum = std::accumulate(v.begin(), v.end(), 0.0);
                    double mean = sum / v.size();
                    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
                    double std = std::sqrt(sq_sum / v.size() - mean * mean);
                    return std::make_pair(mean, std);
                };

                auto [m_sharpe, s_sharpe] = mean_std(sharpes);
                result.mean_sharpe = m_sharpe;
                result.std_sharpe = s_sharpe;

                auto [m_ret, s_ret] = mean_std(returns);
                result.mean_portfolio_return = m_ret;

                auto [m_vol, s_vol] = mean_std(volatilities);
                result.mean_portfolio_volatility = m_vol;

                auto [m_solve, s_solve] = mean_std(solve_times);
                result.mean_solve_time_ms = m_solve;

                auto [m_circuit, s_circuit] = mean_std(circuit_times);
                result.mean_circuit_execution_us = m_circuit;
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

                // nav_series with dates for the viz
                nlohmann::json nav_arr = nlohmann::json::array();
                for (size_t i = 0; i < res.nav_series.size(); ++i)
                {
                    nlohmann::json point;
                    point["nav"] = res.nav_series[i];
                    if (i < res.date_series.size())
                        point["date"] = res.date_series[i];
                    else
                        point["date"] = static_cast<int>(i);
                    nav_arr.push_back(point);
                }
                run["nav_series"] = nav_arr;

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