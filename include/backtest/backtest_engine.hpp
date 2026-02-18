// SPDX-License-Identifier: MIT
#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "data/market_data.hpp"
#include "data/data_loader.hpp"
#include "backtest/portfolio.hpp"
#include "backtest/rebalance_scheduler.hpp"
#include "backtest/trade_logger.hpp"
#include "backtest/transaction_cost_model.hpp"
#include "risk/risk_model_factory.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
#include "analytics/performance_metrics.hpp"

namespace portfolio
{
    namespace backtest
    {

        /**
         * @struct BacktestResult
         * @brief Container for backtest output data and analytics.
         *
         * Holds the NAV series, return series, portfolio snapshots, and trade
         * records produced by a backtest run. Provides both lightweight inline
         * metrics (total_return, sharpe_ratio, max_drawdown) and a bridge to
         * the full analytics layer via compute_analytics().
         */
        struct BacktestResult
        {
            bool success = false;
            std::string message;

            std::vector<PortfolioSnapshot> snapshots;
            std::vector<double> nav_series;
            std::vector<double> return_series;
            std::vector<std::string> dates;

            TradeSummary trade_summary;
            std::vector<TradeRecord> trades;

            // -------------------------------------------------------------------
            // Inline metrics (lightweight, no analytics dependency)
            // -------------------------------------------------------------------

            double total_return() const;
            double annualized_return(int trading_days_per_year = 252) const;
            double annualized_volatility(int trading_days_per_year = 252) const;
            double sharpe_ratio(double risk_free_rate = 0.02,
                                int trading_days_per_year = 252) const;
            double max_drawdown() const;

            // -------------------------------------------------------------------
            // Analytics integration
            // -------------------------------------------------------------------

            /**
             * @brief Construct a full PerformanceMetrics object from this result.
             *
             * Creates a PerformanceMetrics instance initialized from the NAV series,
             * return series, and dates stored in this result. This provides access
             * to the complete analytics suite: VaR/CVaR, drawdown analysis, skewness,
             * kurtosis, Sortino, Calmar, Omega, monthly return tables, and more.
             *
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @return PerformanceMetrics instance for this result.
             * @throws std::invalid_argument If the result contains insufficient data.
             *
             * @note Requires at least 2 NAV observations and a non-empty return series.
             */
            analytics::PerformanceMetrics compute_analytics(
                double risk_free_rate = 0.02,
                int trading_days_per_year = 252) const;

            // -------------------------------------------------------------------
            // Export
            // -------------------------------------------------------------------

            /**
             * @brief Print a summary of performance metrics to stdout.
             *
             * When the result contains sufficient data for analytics (>= 2 NAV
             * observations), prints the full analytics summary including VaR,
             * Sortino, Calmar, Omega, and drawdown detail. Otherwise falls back
             * to the basic inline metrics.
             */
            void print_summary() const;

            /**
             * @brief Export NAV and return series to CSV.
             * @param filepath Path to write the CSV file.
             * @throws std::runtime_error If the file cannot be opened.
             */
            void export_nav_to_csv(const std::string &filepath) const;

            /**
             * @brief Export full analytics to JSON.
             *
             * Computes analytics and returns a JSON string containing all metrics.
             * Includes return metrics, risk metrics, risk-adjusted metrics, and
             * drawdown detail.
             *
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @return JSON string with all analytics metrics.
             * @throws std::invalid_argument If the result contains insufficient data.
             */
            std::string export_analytics_json(
                double risk_free_rate = 0.02,
                int trading_days_per_year = 252) const;

            /**
             * @brief Export full analytics time series to CSV.
             *
             * Writes date, NAV, return, and drawdown columns to CSV via the
             * analytics layer.
             *
             * @param filepath Path to write the CSV file.
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @throws std::runtime_error If the file cannot be opened.
             * @throws std::invalid_argument If the result contains insufficient data.
             */
            void export_analytics_csv(
                const std::string &filepath,
                double risk_free_rate = 0.02,
                int trading_days_per_year = 252) const;
        };

        struct BacktestParams
        {
            double initial_capital = 1000000.0;
            int lookback_window = 252;
            int min_history = 60;
            RebalanceConfig rebalance;
            TransactionCostConfig transaction_costs;
            optimizer::OptimizationConstraints constraints;
            std::string risk_model_type = "ewma";
            double risk_free_rate = 0.02;
            bool verbose = false;

            static BacktestParams from_config(const portfolio::PortfolioConfig &config);
        };

        class BacktestEngine
        {
        public:
            explicit BacktestEngine(const BacktestParams &params);
            ~BacktestEngine() = default;

            BacktestResult run(const MarketData &market_data);
            BacktestResult run(const MarketData &market_data,
                               optimizer::MeanVarianceOptimizer &optimizer);

            const BacktestParams &params() const { return params_; }

        private:
            BacktestParams params_;

            Eigen::MatrixXd extract_window(const MarketData &market_data,
                                           int date_index) const;

            Eigen::MatrixXd estimate_risk(const Eigen::MatrixXd &returns) const;

            Eigen::VectorXd estimate_returns(const Eigen::MatrixXd &returns) const;

            Eigen::VectorXd optimize(const Eigen::VectorXd &expected_returns,
                                     const Eigen::MatrixXd &covariance,
                                     const Eigen::VectorXd &current_weights,
                                     optimizer::MeanVarianceOptimizer &optimizer) const;
        };

    } // namespace backtest
} // namespace portfolio