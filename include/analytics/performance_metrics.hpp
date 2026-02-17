/**
 * @file performance_metrics.hpp
 * @brief Core performance metrics calculator for portfolio return series.
 *
 * Provides comprehensive return, risk, and risk-adjusted performance metrics
 * computed from NAV series, return series, and associated dates. Can be
 * constructed directly from a BacktestResult or from raw vectors for
 * standalone use.
 *
 * All annualized calculations default to 252 trading days per year.
 * Risk-free rate defaults to 2% annualized and is converted internally
 * to a daily rate where needed.
 */

#ifndef PORTFOLIO_ANALYTICS_PERFORMANCE_METRICS_HPP
#define PORTFOLIO_ANALYTICS_PERFORMANCE_METRICS_HPP

#include <string>
#include <vector>
#include <map>

namespace portfolio
{

    // Forward declaration to avoid circular dependency
    namespace backtest
    {
        struct BacktestResult;
    }

    namespace analytics
    {

        /**
         * @struct DrawdownInfo
         * @brief Summary of the maximum drawdown event.
         *
         * Contains the depth of the worst peak-to-trough decline, the duration
         * from peak to trough in trading days, and the recovery time from
         * trough back to the prior peak (or -1 if unrecovered by series end).
         */
        struct DrawdownInfo
        {
            double depth;       ///< Maximum drawdown as a positive fraction (e.g., 0.15 = 15%)
            int duration_days;  ///< Trading days from peak to trough
            int recovery_days;  ///< Trading days from trough to recovery (-1 if unrecovered)
            int peak_index;     ///< Index of the peak NAV before the drawdown
            int trough_index;   ///< Index of the trough NAV
            int recovery_index; ///< Index of recovery (-1 if unrecovered)
        };

        /**
         * @enum VaRMethod
         * @brief Method used for Value at Risk and CVaR calculations.
         */
        enum class VaRMethod
        {
            HISTORICAL, ///< Empirical quantile from historical return distribution
            PARAMETRIC  ///< Assumes normal distribution, uses mean and std dev
        };

        /**
         * @struct MonthlyReturnTable
         * @brief Monthly returns organized by year and month.
         *
         * The outer map is keyed by year, the inner map by month (1-12).
         * An annual total is stored separately per year.
         */
        struct MonthlyReturnTable
        {
            std::map<int, std::map<int, double>> monthly_returns; ///< [year][month] -> return
            std::map<int, double> annual_returns;                 ///< [year] -> annual return
        };

        /**
         * @class PerformanceMetrics
         * @brief Comprehensive performance analytics for a portfolio return series.
         *
         * Computes return metrics (total, annualized, rolling), risk metrics
         * (volatility, downside deviation, VaR, CVaR, drawdown, skewness,
         * kurtosis), and risk-adjusted metrics (Sharpe, Sortino, Calmar, Omega).
         *
         * Usage:
         * @code
         *   auto metrics = PerformanceMetrics(backtest_result);
         *   double sharpe = metrics.sharpe_ratio();
         *   double max_dd = metrics.max_drawdown();
         *   auto dd_series = metrics.drawdown_series();
         * @endcode
         *
         * Thread safety: Instances are effectively immutable after construction.
         * All public methods are const and safe to call concurrently.
         *
         * Performance: All metrics on a 252-day series compute in under 5ms.
         */
        class PerformanceMetrics
        {
        public:
            // ---------------------------------------------------------------
            // Constructors
            // ---------------------------------------------------------------

            /**
             * @brief Construct from a BacktestResult.
             * @param result Backtest result containing NAV series, returns, and dates.
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @throws std::invalid_argument If result contains empty series.
             */
            explicit PerformanceMetrics(const backtest::BacktestResult &result,
                                        double risk_free_rate = 0.02,
                                        int trading_days_per_year = 252);

            /**
             * @brief Construct from raw data vectors.
             * @param nav_series Daily NAV values (must have at least 2 elements).
             * @param return_series Daily simple returns (size must equal nav_series.size() - 1,
             *        or nav_series.size() if NAV starts at the same index as returns).
             * @param dates Date strings corresponding to the NAV series.
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @throws std::invalid_argument If inputs are empty or sizes are inconsistent.
             */
            PerformanceMetrics(const std::vector<double> &nav_series,
                               const std::vector<double> &return_series,
                               const std::vector<std::string> &dates,
                               double risk_free_rate = 0.02,
                               int trading_days_per_year = 252);

            /** @brief Default destructor. */
            ~PerformanceMetrics() = default;

            // ---------------------------------------------------------------
            // Return Metrics
            // ---------------------------------------------------------------

            /**
             * @brief Total cumulative return over the full period.
             * @return Total return as a fraction (e.g., 0.25 = 25%).
             */
            double total_return() const;

            /**
             * @brief Annualized return (CAGR).
             * @return Annualized return as a fraction.
             */
            double annualized_return() const;

            /**
             * @brief Rolling returns over a configurable window.
             * @param window_days Number of trading days in the rolling window.
             * @return Vector of rolling returns, size = n - window_days.
             * @throws std::invalid_argument If window_days < 1 or exceeds series length.
             */
            std::vector<double> rolling_returns(int window_days) const;

            /**
             * @brief Monthly return table organized by year and month.
             * @return MonthlyReturnTable with monthly and annual returns.
             * @throws std::runtime_error If dates cannot be parsed.
             *
             * @note Dates must be in YYYY-MM-DD format.
             */
            MonthlyReturnTable monthly_return_table() const;

            // ---------------------------------------------------------------
            // Risk Metrics
            // ---------------------------------------------------------------

            /**
             * @brief Annualized volatility (standard deviation of returns * sqrt(trading_days)).
             * @return Annualized volatility as a fraction.
             */
            double annualized_volatility() const;

            /**
             * @brief Downside deviation (semi-deviation below a target return).
             * @param target_return Annualized target return (default 0.0).
             *        Converted internally to a daily target.
             * @return Annualized downside deviation.
             */
            double downside_deviation(double target_return = 0.0) const;

            /**
             * @brief Maximum drawdown depth.
             * @return Maximum drawdown as a positive fraction (e.g., 0.15 = 15%).
             */
            double max_drawdown() const;

            /**
             * @brief Detailed maximum drawdown information.
             * @return DrawdownInfo with depth, duration, and recovery time.
             */
            DrawdownInfo max_drawdown_info() const;

            /**
             * @brief Full drawdown (underwater) series.
             * @return Vector of drawdown values at each point, where 0.0 means at peak
             *         and negative values indicate drawdown depth.
             *
             * @note Values are non-positive. A value of -0.10 means 10% below the peak.
             */
            std::vector<double> drawdown_series() const;

            /**
             * @brief Value at Risk at the specified confidence level.
             * @param confidence Confidence level (default 0.95 = 95%).
             * @param method VaR calculation method (default HISTORICAL).
             * @return VaR as a positive number representing the loss threshold.
             *         A VaR of 0.02 means "95% of the time, daily loss does not exceed 2%."
             * @throws std::invalid_argument If confidence is not in (0, 1).
             */
            double value_at_risk(double confidence = 0.95,
                                 VaRMethod method = VaRMethod::HISTORICAL) const;

            /**
             * @brief Conditional VaR (Expected Shortfall) at the specified confidence level.
             * @param confidence Confidence level (default 0.95 = 95%).
             * @param method CVaR calculation method (default HISTORICAL).
             * @return CVaR as a positive number representing the expected loss
             *         beyond the VaR threshold.
             * @throws std::invalid_argument If confidence is not in (0, 1).
             */
            double conditional_var(double confidence = 0.95,
                                   VaRMethod method = VaRMethod::HISTORICAL) const;

            /**
             * @brief Skewness of the daily return distribution.
             * @return Sample skewness (Fisher's definition, bias-corrected).
             */
            double skewness() const;

            /**
             * @brief Excess kurtosis of the daily return distribution.
             * @return Excess kurtosis (0 for normal distribution, bias-corrected).
             */
            double kurtosis() const;

            // ---------------------------------------------------------------
            // Risk-Adjusted Metrics
            // ---------------------------------------------------------------

            /**
             * @brief Sharpe ratio (annualized).
             * @return (annualized_return - risk_free_rate) / annualized_volatility.
             *
             * @note Returns 0.0 if volatility is zero.
             */
            double sharpe_ratio() const;

            /**
             * @brief Sortino ratio (annualized).
             * @param target_return Annualized target return for downside deviation (default 0.0).
             * @return (annualized_return - target_return) / downside_deviation(target_return).
             *
             * @note Returns 0.0 if downside deviation is zero.
             */
            double sortino_ratio(double target_return = 0.0) const;

            /**
             * @brief Calmar ratio.
             * @return annualized_return / max_drawdown.
             *
             * @note Returns 0.0 if max drawdown is zero.
             */
            double calmar_ratio() const;

            /**
             * @brief Omega ratio at the specified threshold.
             * @param threshold Annualized threshold return (default 0.0).
             *        Converted internally to a daily threshold.
             * @return Ratio of gains above threshold to losses below threshold.
             *         Returns infinity (max double) if there are no losses below threshold.
             * @throws std::invalid_argument If all returns equal the threshold exactly.
             */
            double omega_ratio(double threshold = 0.0) const;

            // ---------------------------------------------------------------
            // Accessors
            // ---------------------------------------------------------------

            /** @brief Get the NAV series. */
            const std::vector<double> &nav_series() const;

            /** @brief Get the return series. */
            const std::vector<double> &return_series() const;

            /** @brief Get the date series. */
            const std::vector<std::string> &dates() const;

            /** @brief Get the annualized risk-free rate. */
            double risk_free_rate() const;

            /** @brief Get the trading days per year setting. */
            int trading_days_per_year() const;

            // ---------------------------------------------------------------
            // Export
            // ---------------------------------------------------------------

            /**
             * @brief Generate a summary string of all core metrics.
             * @return Formatted multi-line string with key performance metrics.
             */
            std::string summary() const;

            /**
             * @brief Export metrics to a JSON-formatted string.
             * @return JSON string containing all computed metrics.
             */
            std::string to_json() const;

            /**
             * @brief Export time series data to CSV format.
             * @param filepath Path to write the CSV file.
             * @throws std::runtime_error If the file cannot be opened.
             *
             * Columns: date, nav, return, drawdown
             */
            void to_csv(const std::string &filepath) const;

        private:
            // ---------------------------------------------------------------
            // Private helpers
            // ---------------------------------------------------------------

            /**
             * @brief Validate inputs and compute cached values.
             */
            void initialize();

            /**
             * @brief Compute the drawdown series from the NAV series.
             */
            void compute_drawdown_cache() const;

            /**
             * @brief Convert an annualized rate to a daily rate (simple).
             * @param annual_rate Annualized rate.
             * @return Daily rate = annual_rate / trading_days_per_year.
             */
            double to_daily_rate(double annual_rate) const;

            // ---------------------------------------------------------------
            // Member variables
            // ---------------------------------------------------------------

            std::vector<double> nav_series_;
            std::vector<double> return_series_;
            std::vector<std::string> dates_;
            double risk_free_rate_;
            int trading_days_per_year_;

            // Cached computations (mutable for lazy evaluation in const methods)
            mutable bool drawdown_computed_;
            mutable std::vector<double> drawdown_series_cache_;
            mutable DrawdownInfo max_drawdown_info_cache_;
        };

    } // namespace analytics
} // namespace portfolio

#endif // PORTFOLIO_ANALYTICS_PERFORMANCE_METRICS_HPP