// SPDX-License-Identifier: MIT

#include "backtest/backtest_engine.hpp"
#include "backtest/transaction_cost_model.hpp"
#include "analytics/performance_metrics.hpp"

#include <fstream>
#include <iostream>
#include <numeric>
#include <cmath>

namespace portfolio
{
    namespace backtest
    {

        using namespace portfolio::risk;
        using namespace portfolio::optimizer;

        // ------------------------- BacktestParams -------------------------------
        BacktestParams BacktestParams::from_config(const portfolio::PortfolioConfig &config)
        {
            BacktestParams p;
            p.initial_capital = config.backtest.initial_capital;
            p.lookback_window = config.backtest.lookback_window;
            p.min_history = config.backtest.min_history;
            p.transaction_costs = TransactionCostConfig::from_json(nlohmann::json::object());
            try
            {
                p.transaction_costs = TransactionCostConfig::from_json(nlohmann::json{
                    {"commission_rate", config.backtest.commission_rate},
                    {"slippage_bps", config.backtest.slippage_bps}});
            }
            catch (...)
            {
            }
            p.rebalance = RebalanceConfig::from_string(config.backtest.rebalance_frequency);
            p.constraints = OptimizationConstraints::from_json(nlohmann::json::object());
            p.risk_model_type = config.risk_model.type;
            p.risk_free_rate = config.backtest.cash_return;
            return p;
        }

        // ------------------------- BacktestEngine -------------------------------
        BacktestEngine::BacktestEngine(const BacktestParams &params) : params_(params) {}

        Eigen::MatrixXd BacktestEngine::extract_window(const MarketData &market_data, int date_index) const
        {
            // returns matrix has rows = dates-1
            Eigen::MatrixXd all_returns = market_data.calculate_returns();
            int max_row = all_returns.rows(); // = num_dates - 1

            int end_row = date_index; // we want rows [date_index - lookback, date_index)
            int look = params_.lookback_window;
            int start_row = std::max(0, end_row - look);
            if (start_row >= end_row)
            {
                // empty window -> return full available
                start_row = 0;
                end_row = std::min(max_row, end_row);
            }

            int rows = end_row - start_row;
            if (rows <= 0)
                return Eigen::MatrixXd(0, market_data.num_assets());

            Eigen::MatrixXd window(rows, static_cast<int>(market_data.num_assets()));
            window = all_returns.block(start_row, 0, rows, static_cast<int>(market_data.num_assets()));
            return window;
        }

        Eigen::MatrixXd BacktestEngine::estimate_risk(const Eigen::MatrixXd &returns) const
        {
            if (returns.rows() == 0)
                return Eigen::MatrixXd();
            portfolio::risk::RiskModelConfig cfg;
            cfg.type = params_.risk_model_type;
            auto model = RiskModelFactory::create(cfg);
            Eigen::MatrixXd cov = model->estimate_covariance(returns);
            // annualize
            return cov * 252.0;
        }

        Eigen::VectorXd BacktestEngine::estimate_returns(const Eigen::MatrixXd &returns) const
        {
            if (returns.rows() == 0)
                return Eigen::VectorXd();
            Eigen::VectorXd mu = returns.colwise().mean();
            // annualize mean
            return mu * 252.0;
        }

        Eigen::VectorXd BacktestEngine::optimize(const Eigen::VectorXd &expected_returns,
                                                 const Eigen::MatrixXd &covariance,
                                                 const Eigen::VectorXd &current_weights,
                                                 MeanVarianceOptimizer &optimizer) const
        {
            OptimizationConstraints constraints = params_.constraints;
            try
            {
                auto res = optimizer.optimize(expected_returns, covariance, constraints, current_weights);
                if (!res.success)
                    throw std::runtime_error(res.message);
                return res.weights;
            }
            catch (const std::exception &e)
            {
                throw;
            }
        }

        Eigen::VectorXd BacktestEngine::optimize_quantum(const Eigen::VectorXd &expected_returns,
                                                         const Eigen::MatrixXd &covariance,
                                                         const Eigen::VectorXd &current_weights,
                                                         quantum::QuantumOptimizer &optimizer) const
        {
            try
            {
                return optimizer.optimize(covariance, expected_returns, params_.constraints);
            }
            catch (const std::exception &e)
            {
                throw;
            }
        }

        Eigen::VectorXd BacktestEngine::project_weights(
            const Eigen::VectorXd &raw_weights,
            const optimizer::OptimizationConstraints &constraints) const
        {
            int n = raw_weights.size();
            if (n == 0) return raw_weights;

            Eigen::VectorXd w = raw_weights;
            double w_min = constraints.min_weight > 0.0 ? constraints.min_weight : 0.0;
            double w_max = constraints.max_weight < 1.0 ? constraints.max_weight : 1.0;

            // Clip to [0, w_max]
            w = w.cwiseMax(0.0).cwiseMin(w_max);

            // Zero out positions below min_weight (treat as no position)
            for (int i = 0; i < n; ++i)
                if (w[i] > 0.0 && w[i] < w_min) w[i] = 0.0;

            // Fallback when everything zeroed: equal weight across the maximum
            // number of slots that can each hold >= w_min (avoids a fallback that
            // itself violates the min-weight constraint).
            double total = w.sum();
            if (total <= 1e-10) {
                int k = (w_min > 1e-10) ? std::min(n, (int)(1.0 / w_min)) : n;
                w = Eigen::VectorXd::Zero(n);
                double eq = 1.0 / k;
                for (int i = 0; i < k; ++i) w[i] = eq;
                return w;
            }

            // Iterative clip-and-redistribute: after renorm a capped asset can push
            // others above w_max, so we clip and spread the deficit to uncapped slots.
            for (int iter = 0; iter < 20; ++iter) {
                w /= w.sum();
                w = w.cwiseMin(w_max);
                double deficit = 1.0 - w.sum();
                if (std::abs(deficit) < 1e-10) break;
                int uncapped = 0;
                for (int i = 0; i < n; ++i)
                    if (w[i] < w_max - 1e-12) ++uncapped;
                if (uncapped == 0) break;
                double share = deficit / uncapped;
                for (int i = 0; i < n; ++i)
                    if (w[i] < w_max - 1e-12) w[i] += share;
            }

            total = w.sum();
            if (total > 1e-10) w /= total;

            return w;
        }

        // ------------------------- BacktestResult helpers -----------------------
        double BacktestResult::total_return() const
        {
            if (nav_series.empty())
                return 0.0;
            double init = nav_series.front();
            double last = nav_series.back();
            if (init == 0.0)
                return 0.0;
            return last / init - 1.0;
        }

        double BacktestResult::annualized_return(int trading_days_per_year) const
        {
            if (nav_series.size() < 2)
                return 0.0;
            double tot = total_return();
            double n_days = static_cast<double>(nav_series.size());
            return std::pow(1.0 + tot, static_cast<double>(trading_days_per_year) / n_days) - 1.0;
        }

        double BacktestResult::annualized_volatility(int trading_days_per_year) const
        {
            if (return_series.empty())
                return 0.0;
            double mean = std::accumulate(return_series.begin(), return_series.end(), 0.0) / return_series.size();
            double sumsq = 0.0;
            for (double r : return_series)
                sumsq += (r - mean) * (r - mean);
            double variance = sumsq / static_cast<double>(return_series.size());
            return std::sqrt(variance) * std::sqrt(static_cast<double>(trading_days_per_year));
        }

        double BacktestResult::sharpe_ratio(double risk_free_rate, int trading_days_per_year) const
        {
            double ar = annualized_return(trading_days_per_year);
            double av = annualized_volatility(trading_days_per_year);
            if (av <= 0.0)
                return 0.0;
            return (ar - risk_free_rate) / av;
        }

        double BacktestResult::max_drawdown() const
        {
            if (nav_series.empty())
                return 0.0;
            double peak = nav_series.front();
            double max_dd = 0.0;
            for (double v : nav_series)
            {
                if (v > peak)
                    peak = v;
                double dd = (peak - v) / peak;
                if (dd > max_dd)
                    max_dd = dd;
            }
            return max_dd;
        }

        void BacktestResult::print_summary() const
        {
            std::cout << "Backtest success=" << success << " message='" << message << "'\n";
            std::cout << "Total return: " << total_return() << " Annualized: " << annualized_return() << "\n";
            std::cout << "Ann vol: " << annualized_volatility() << " Sharpe: " << sharpe_ratio() << "\n";
            std::cout << "Max drawdown: " << max_drawdown() << "\n";
        }

        void BacktestResult::export_nav_to_csv(const std::string &filepath) const
        {
            std::ofstream out(filepath);
            if (!out)
                throw std::runtime_error("unable to open file for writing: " + filepath);
            out << "date,nav,daily_return,cumulative_return\n";
            for (size_t i = 0; i < nav_series.size() && i < dates.size(); ++i)
            {
                double cum = 0.0;
                if (!nav_series.empty())
                    cum = nav_series[i] / nav_series.front() - 1.0;
                out << dates[i] << "," << nav_series[i] << "," << (i < return_series.size() ? return_series[i] : 0.0) << "," << cum << "\n";
            }
        }

        // ------------------------- Run methods ---------------------------------
        BacktestResult BacktestEngine::run(const MarketData &market_data)
        {
            MeanVarianceOptimizer optimizer(optimizer::ObjectiveType::MAX_SHARPE, params_.risk_free_rate);
            return run(market_data, optimizer);
        }

        BacktestResult BacktestEngine::run(const MarketData &market_data, MeanVarianceOptimizer &optimizer)
        {
            BacktestResult result;

            size_t N = market_data.num_dates();
            if (static_cast<int>(N) < params_.min_history)
            {
                throw std::invalid_argument("market_data has fewer dates than min_history");
            }

            std::vector<std::string> dates = market_data.get_dates();
            std::vector<std::string> tickers = market_data.get_tickers();

            const Eigen::MatrixXd &prices_matrix = market_data.get_Prices();

            Portfolio portfolio(params_.initial_capital, tickers);
            TransactionCostModel cost_model(params_.transaction_costs);
            RebalanceScheduler scheduler(params_.rebalance);
            TradeLogger logger;

            // Precompute returns matrix once
            Eigen::MatrixXd all_returns = market_data.calculate_returns();

            for (int date_idx = 0; date_idx < static_cast<int>(N); ++date_idx)
            {
                const std::string &date = dates[static_cast<size_t>(date_idx)];
                Eigen::VectorXd prices = prices_matrix.row(date_idx).transpose();

                portfolio.update_prices(date, prices);

                // Rebalance logic
                if (date_idx >= params_.min_history)
                {
                    // Extract window and estimate
                    Eigen::MatrixXd window = extract_window(market_data, date_idx);
                    if (window.rows() > 0)
                    {
                        if (window.rows() < 2)
                        {
                            // not enough observations for covariance estimation
                            if (params_.verbose)
                                std::cerr << "Skipping optimization on " << date << " - insufficient return observations\n";
                        }
                        else
                        {
                            Eigen::VectorXd expected = estimate_returns(window);
                            Eigen::MatrixXd cov = estimate_risk(window);

                            Eigen::VectorXd pre_weights = portfolio.current_weights();
                            Eigen::VectorXd target_weights;
                            try
                            {
                                auto res = optimizer.optimize(expected, cov, params_.constraints, pre_weights);
                                if (!res.success)
                                {
                                    std::cerr << "Optimizer failed on " << date << ": " << res.message << "\n";
                                }
                                else
                                {
                                    target_weights = res.weights;
                                    bool do_rebalance = scheduler.should_rebalance(date, pre_weights, target_weights);
                                    if (do_rebalance)
                                    {
                                        // Execute trades
                                        Eigen::VectorXd shares_traded = portfolio.set_target_weights(target_weights, prices);

                                        // Build costs vector
                                        std::vector<TradeCost> costs;
                                        costs.reserve(static_cast<size_t>(shares_traded.size()));
                                        for (int i = 0; i < shares_traded.size(); ++i)
                                        {
                                            TradeOrder o{tickers[static_cast<size_t>(i)], shares_traded[i], prices[i]};
                                            costs.push_back(cost_model.calculate_cost(o));
                                        }

                                        double total_cost = cost_model.calculate_total_cost(shares_traded, prices, tickers);
                                        if (total_cost > 0.0)
                                        {
                                            try
                                            {
                                                portfolio.deduct_costs(total_cost);
                                            }
                                            catch (...)
                                            {
                                            }
                                        }

                                        Eigen::VectorXd post_weights = portfolio.current_weights();
                                        logger.log_rebalance(date, tickers, shares_traded, prices, costs, pre_weights, post_weights);
                                        scheduler.record_rebalance(date);
                                    }
                                }
                            }
                            catch (const std::exception &e)
                            {
                                std::cerr << "Optimization exception on " << date << ": " << e.what() << "\n";
                            }
                        }
                    }
                } // <-- closes if (false && ...) — THIS WAS THE MISSING BRACE

                // Record snapshot starting from second date to align with returns
                if (date_idx > 0)
                {
                    portfolio.record_snapshot();
                    const auto &hist = portfolio.history();
                    const PortfolioSnapshot &s = hist.back();
                    result.snapshots.push_back(s);
                    result.nav_series.push_back(s.nav);
                    result.return_series.push_back(s.daily_return);
                    result.dates.push_back(s.date);
                }
            } // <-- closes for loop

            result.success = true;
            result.trade_summary = logger.get_summary();
            result.trades = logger.trades();
            result.message = "completed";
            return result;
        }

        BacktestResult BacktestEngine::run(const MarketData &market_data, quantum::QuantumOptimizer &optimizer)
        {
            BacktestResult result;

            size_t N = market_data.num_dates();
            if (static_cast<int>(N) < params_.min_history)
            {
                throw std::invalid_argument("market_data has fewer dates than min_history");
            }

            std::vector<std::string> dates = market_data.get_dates();
            std::vector<std::string> tickers = market_data.get_tickers();

            const Eigen::MatrixXd &prices_matrix = market_data.get_Prices();

            Portfolio portfolio(params_.initial_capital, tickers);
            TransactionCostModel cost_model(params_.transaction_costs);
            RebalanceScheduler scheduler(params_.rebalance);
            TradeLogger logger;

            // Precompute returns matrix once
            Eigen::MatrixXd all_returns = market_data.calculate_returns();

            for (int date_idx = 0; date_idx < static_cast<int>(N); ++date_idx)
            {
                const std::string &date = dates[static_cast<size_t>(date_idx)];
                Eigen::VectorXd prices = prices_matrix.row(date_idx).transpose();

                portfolio.update_prices(date, prices);

                // Rebalance logic
                if (date_idx >= params_.min_history)
                {
                    // Extract window and estimate
                    Eigen::MatrixXd window = extract_window(market_data, date_idx);
                    if (window.rows() > 0)
                    {
                        if (window.rows() < 2)
                        {
                            // not enough observations for covariance estimation
                            if (params_.verbose)
                                std::cerr << "Skipping optimization on " << date << " - insufficient return observations\n";
                        }
                        else if (scheduler.is_calendar_trigger(date))
                        {
                            // Only invoke the quantum optimizer on scheduled rebalance dates.
                            // Calling it every day would spawn two Python subprocesses per
                            // trading day (~252/year) instead of per rebalancing period.
                            Eigen::VectorXd expected = estimate_returns(window);
                            Eigen::MatrixXd cov = estimate_risk(window);

                            Eigen::VectorXd pre_weights = portfolio.current_weights();
                            try
                            {
                                Eigen::VectorXd raw_weights = optimize_quantum(expected, cov, pre_weights, optimizer);
                                Eigen::VectorXd target_weights = project_weights(raw_weights, params_.constraints);

                                Eigen::VectorXd shares_traded = portfolio.set_target_weights(target_weights, prices);

                                std::vector<TradeCost> costs;
                                costs.reserve(static_cast<size_t>(shares_traded.size()));
                                for (int i = 0; i < shares_traded.size(); ++i)
                                {
                                    TradeOrder o{tickers[static_cast<size_t>(i)], shares_traded[i], prices[i]};
                                    costs.push_back(cost_model.calculate_cost(o));
                                }

                                double total_cost = cost_model.calculate_total_cost(shares_traded, prices, tickers);
                                if (total_cost > 0.0)
                                {
                                    try { portfolio.deduct_costs(total_cost); } catch (...) {}
                                }

                                Eigen::VectorXd post_weights = portfolio.current_weights();
                                logger.log_rebalance(date, tickers, shares_traded, prices, costs, pre_weights, post_weights);
                                scheduler.record_rebalance(date);
                            }
                            catch (const std::exception &e)
                            {
                                std::cerr << "Optimization exception on " << date << ": " << e.what() << "\n";
                            }
                        }
                    }
                }

                // Record snapshot starting from second date to align with returns
                if (date_idx > 0)
                {
                    portfolio.record_snapshot();
                    const auto &hist = portfolio.history();
                    const PortfolioSnapshot &s = hist.back();
                    result.snapshots.push_back(s);
                    result.nav_series.push_back(s.nav);
                    result.return_series.push_back(s.daily_return);
                    result.dates.push_back(s.date);
                }
            } // <-- closes for loop

            result.success = true;
            result.trade_summary = logger.get_summary();
            result.trades = logger.trades();
            result.message = "completed";
            return result;
        }

        analytics::PerformanceMetrics BacktestResult::compute_analytics(
            double risk_free_rate, int trading_days_per_year) const
        {
            return analytics::PerformanceMetrics(*this, risk_free_rate, trading_days_per_year);
        }

    } // namespace backtest
} // namespace portfolio