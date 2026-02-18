/**
 * @file rolling_statistics.hpp
 * @brief Rolling window computations for portfolio time series.
 *
 * Computes rolling versions of key performance metrics over configurable
 * windows, producing vectors aligned with the original date series.
 * Supports rolling volatility, Sharpe ratio, beta, drawdown, and
 * generic user-defined rolling functions.
 *
 * Output vectors have size (n - window + 1) where n is the input
 * series length, and are aligned to the end of each window (i.e.,
 * result[0] corresponds to the window ending at index (window - 1)).
 */

#ifndef PORTFOLIO_ANALYTICS_ROLLING_STATISTICS_HPP
#define PORTFOLIO_ANALYTICS_ROLLING_STATISTICS_HPP

#include <functional>
#include <string>
#include <vector>

namespace portfolio
{
    namespace analytics
    {

        /**
         * @struct RollingConfig
         * @brief Configuration for rolling window calculations.
         */
        struct RollingConfig
        {
            int window_days;           ///< Rolling window size in trading days
            int trading_days_per_year; ///< Trading days per year for annualization (default 252)
            double risk_free_rate;     ///< Annualized risk-free rate (default 0.02)
            int min_periods;           ///< Minimum observations required in window (default = window_days)

            /**
             * @brief Construct with defaults.
             * @param window Rolling window size in trading days.
             */
            explicit RollingConfig(int window = 63)
                : window_days(window), trading_days_per_year(252), risk_free_rate(0.02), min_periods(window) {}
        };

        /**
         * @class RollingStatistics
         * @brief Computes rolling metrics over configurable windows.
         *
         * Provides rolling versions of standard performance and risk metrics.
         * All output vectors are aligned so that result[i] contains the metric
         * computed over the window ending at index (i + window - 1) of the
         * input series.
         *
         * Usage:
         * @code
         *   RollingConfig config(63);  // Quarterly rolling window
         *   RollingStatistics rolling(return_series, config);
         *   auto vol = rolling.volatility();
         *   auto sharpe = rolling.sharpe_ratio();
         * @endcode
         *
         * Thread safety: Instances are effectively immutable after construction.
         *
         * Performance: Rolling Sharpe over 1000 days with 252-day window
         * completes in under 10ms.
         */
        class RollingStatistics
        {
        public:
            // ---------------------------------------------------------------
            // Constructors
            // ---------------------------------------------------------------

            /**
             * @brief Construct from a return series and configuration.
             * @param return_series Daily simple returns.
             * @param config Rolling window configuration.
             * @throws std::invalid_argument If return_series is empty or
             *         shorter than config.window_days.
             */
            RollingStatistics(const std::vector<double> &return_series,
                              const RollingConfig &config);

            /**
             * @brief Construct with a default configuration.
             * @param return_series Daily simple returns.
             * @param window_days Rolling window size (default 63).
             * @throws std::invalid_argument If return_series is empty or
             *         shorter than window_days.
             */
            explicit RollingStatistics(const std::vector<double> &return_series,
                                       int window_days = 63);

            /** @brief Default destructor. */
            ~RollingStatistics() = default;

            // ---------------------------------------------------------------
            // Rolling Metrics
            // ---------------------------------------------------------------

            /**
             * @brief Rolling annualized volatility.
             * @return Vector of rolling volatility values.
             *         Size = return_series.size() - window + 1.
             */
            std::vector<double> volatility() const;

            /**
             * @brief Rolling annualized mean return.
             * @return Vector of rolling mean returns.
             */
            std::vector<double> mean_return() const;

            /**
             * @brief Rolling Sharpe ratio (annualized).
             * @return Vector of rolling Sharpe ratios.
             */
            std::vector<double> sharpe_ratio() const;

            /**
             * @brief Rolling Sortino ratio (annualized).
             * @param target_return Annualized target return (default 0.0).
             * @return Vector of rolling Sortino ratios.
             */
            std::vector<double> sortino_ratio(double target_return = 0.0) const;

            /**
             * @brief Rolling maximum drawdown within each window.
             * @return Vector of rolling max drawdown values (positive fractions).
             */
            std::vector<double> max_drawdown() const;

            /**
             * @brief Rolling beta against a benchmark.
             * @param benchmark_returns Daily simple returns of the benchmark
             *        (must have the same length as the portfolio return series).
             * @return Vector of rolling beta values.
             * @throws std::invalid_argument If benchmark_returns size does not match.
             */
            std::vector<double> beta(const std::vector<double> &benchmark_returns) const;

            /**
             * @brief Rolling tracking error against a benchmark.
             * @param benchmark_returns Daily simple returns of the benchmark
             *        (must have the same length as the portfolio return series).
             * @return Vector of rolling annualized tracking error values.
             * @throws std::invalid_argument If benchmark_returns size does not match.
             */
            std::vector<double> tracking_error(const std::vector<double> &benchmark_returns) const;

            /**
             * @brief Rolling skewness.
             * @return Vector of rolling skewness values.
             */
            std::vector<double> skewness() const;

            /**
             * @brief Rolling excess kurtosis.
             * @return Vector of rolling kurtosis values.
             */
            std::vector<double> kurtosis() const;

            // ---------------------------------------------------------------
            // Generic Rolling Application
            // ---------------------------------------------------------------

            /**
             * @brief Apply a custom function to each rolling window.
             * @param func Function that takes a vector slice (window) and returns a double.
             * @return Vector of function results over each rolling window.
             *
             * @note The function receives a vector of size window_days for each window.
             */
            std::vector<double> apply(
                const std::function<double(const std::vector<double> &)> &func) const;

            // ---------------------------------------------------------------
            // Accessors
            // ---------------------------------------------------------------

            /** @brief Get the rolling window configuration. */
            const RollingConfig &config() const;

            /** @brief Get the number of output values (series_length - window + 1). */
            int output_size() const;

            /**
             * @brief Get the dates corresponding to rolling output values.
             * @param dates Full date series aligned with the return series.
             * @return Subset of dates aligned with rolling output (end of each window).
             * @throws std::invalid_argument If dates size does not match return series size.
             */
            std::vector<std::string> aligned_dates(
                const std::vector<std::string> &dates) const;

        private:
            // ---------------------------------------------------------------
            // Private helpers
            // ---------------------------------------------------------------

            /**
             * @brief Compute rolling window statistics using incremental updates
             *        where possible for performance.
             * @param func Function to apply per window.
             * @return Vector of results.
             */
            std::vector<double> rolling_apply_internal(
                const std::function<double(const double *, int)> &func) const;

            // ---------------------------------------------------------------
            // Member variables
            // ---------------------------------------------------------------

            std::vector<double> return_series_;
            RollingConfig config_;
        };

    } // namespace analytics
} // namespace portfolio

#endif // PORTFOLIO_ANALYTICS_ROLLING_STATISTICS_HPP