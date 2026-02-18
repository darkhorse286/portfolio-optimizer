/**
 * @file benchmark_analysis.hpp
 * @brief Relative performance analysis against a benchmark.
 *
 * Computes CAPM regression statistics (alpha, beta, R-squared),
 * tracking error, information ratio, and up/down capture ratios
 * by comparing portfolio returns against a benchmark return series.
 *
 * The CAPM regression model is:
 *   R_p - R_f = alpha + beta * (R_b - R_f) + epsilon
 *
 * where R_p is portfolio return, R_b is benchmark return, and
 * R_f is the risk-free rate.
 */

#ifndef PORTFOLIO_ANALYTICS_BENCHMARK_ANALYSIS_HPP
#define PORTFOLIO_ANALYTICS_BENCHMARK_ANALYSIS_HPP

#include <string>
#include <vector>

namespace portfolio
{
    namespace analytics
    {

        /**
         * @struct RegressionResult
         * @brief Results of the CAPM regression.
         */
        struct RegressionResult
        {
            double alpha;         ///< Jensen's alpha (annualized intercept)
            double beta;          ///< Slope of regression (market sensitivity)
            double r_squared;     ///< Coefficient of determination (0 to 1)
            double alpha_daily;   ///< Daily alpha (raw intercept)
            double residual_std;  ///< Standard deviation of regression residuals
            int num_observations; ///< Number of data points in the regression
        };

        /**
         * @struct CaptureRatios
         * @brief Up-capture and down-capture ratios vs benchmark.
         */
        struct CaptureRatios
        {
            double up_capture;   ///< Ratio of portfolio up-market return to benchmark up-market return
            double down_capture; ///< Ratio of portfolio down-market return to benchmark down-market return
            int up_periods;      ///< Number of periods where benchmark return > 0
            int down_periods;    ///< Number of periods where benchmark return <= 0
        };

        /**
         * @class BenchmarkAnalysis
         * @brief Computes relative performance metrics against a benchmark.
         *
         * Performs CAPM regression and computes a suite of relative metrics
         * for evaluating how a portfolio performs compared to a benchmark
         * index or strategy.
         *
         * Usage:
         * @code
         *   BenchmarkAnalysis bench(portfolio_returns, benchmark_returns);
         *   double alpha = bench.alpha();
         *   double beta = bench.beta();
         *   double info_ratio = bench.information_ratio();
         *   auto capture = bench.capture_ratios();
         * @endcode
         *
         * Thread safety: Instances are effectively immutable after construction.
         *
         * Performance: All metrics on a 252-day series compute in under 1ms.
         */
        class BenchmarkAnalysis
        {
        public:
            // ---------------------------------------------------------------
            // Constructors
            // ---------------------------------------------------------------

            /**
             * @brief Construct from portfolio and benchmark return series.
             * @param portfolio_returns Daily simple returns of the portfolio.
             * @param benchmark_returns Daily simple returns of the benchmark.
             * @param risk_free_rate Annualized risk-free rate (default 0.02 = 2%).
             * @param trading_days_per_year Number of trading days per year (default 252).
             * @throws std::invalid_argument If return series are empty, have different
             *         sizes, or contain fewer than 3 observations.
             */
            BenchmarkAnalysis(const std::vector<double> &portfolio_returns,
                              const std::vector<double> &benchmark_returns,
                              double risk_free_rate = 0.02,
                              int trading_days_per_year = 252);

            /** @brief Default destructor. */
            ~BenchmarkAnalysis() = default;

            // ---------------------------------------------------------------
            // CAPM Regression Metrics
            // ---------------------------------------------------------------

            /**
             * @brief Jensen's alpha (annualized).
             * @return Annualized alpha from the CAPM regression.
             *
             * @note Alpha = daily_alpha * trading_days_per_year.
             */
            double alpha() const;

            /**
             * @brief CAPM beta (market sensitivity).
             * @return Beta coefficient from the regression.
             *
             * A beta of 1.0 indicates market-equivalent sensitivity.
             */
            double beta() const;

            /**
             * @brief R-squared of the CAPM regression.
             * @return Fraction of portfolio variance explained by the benchmark (0 to 1).
             */
            double r_squared() const;

            /**
             * @brief Full regression results.
             * @return RegressionResult with all regression statistics.
             */
            const RegressionResult &regression() const;

            // ---------------------------------------------------------------
            // Tracking and Relative Risk
            // ---------------------------------------------------------------

            /**
             * @brief Tracking error (annualized).
             * @return Annualized standard deviation of excess returns
             *         (portfolio returns minus benchmark returns).
             */
            double tracking_error() const;

            /**
             * @brief Information ratio.
             * @return Annualized excess return / tracking error.
             *
             * @note Returns 0.0 if tracking error is zero.
             */
            double information_ratio() const;

            /**
             * @brief Active return (annualized excess return over benchmark).
             * @return Annualized mean of (portfolio_return - benchmark_return).
             */
            double active_return() const;

            // ---------------------------------------------------------------
            // Capture Ratios
            // ---------------------------------------------------------------

            /**
             * @brief Up-capture ratio.
             * @return Mean portfolio return in up-market periods / mean benchmark
             *         return in up-market periods. A value > 1 means the portfolio
             *         captures more upside than the benchmark.
             *
             * @note Returns 0.0 if there are no up-market periods.
             */
            double up_capture_ratio() const;

            /**
             * @brief Down-capture ratio.
             * @return Mean portfolio return in down-market periods / mean benchmark
             *         return in down-market periods. A value < 1 means the portfolio
             *         loses less in down markets than the benchmark.
             *
             * @note Returns 0.0 if there are no down-market periods.
             */
            double down_capture_ratio() const;

            /**
             * @brief Full capture ratio details.
             * @return CaptureRatios with both ratios and period counts.
             */
            const CaptureRatios &capture_ratios() const;

            // ---------------------------------------------------------------
            // Excess Return Series
            // ---------------------------------------------------------------

            /**
             * @brief Get the excess return series (portfolio - benchmark).
             * @return Vector of daily excess returns.
             */
            const std::vector<double> &excess_returns() const;

            // ---------------------------------------------------------------
            // Export
            // ---------------------------------------------------------------

            /**
             * @brief Generate a formatted summary of benchmark-relative metrics.
             * @return Formatted multi-line string with all relative metrics.
             */
            std::string summary() const;

            /**
             * @brief Export benchmark analysis to a JSON-formatted string.
             * @return JSON string containing all relative metrics.
             */
            std::string to_json() const;

        private:
            // ---------------------------------------------------------------
            // Private helpers
            // ---------------------------------------------------------------

            /**
             * @brief Run the CAPM regression.
             */
            void run_regression();

            /**
             * @brief Compute tracking error and information ratio.
             */
            void compute_tracking_metrics();

            /**
             * @brief Compute up-capture and down-capture ratios.
             */
            void compute_capture_ratios();

            /**
             * @brief Convert an annualized rate to a daily rate.
             * @param annual_rate Annualized rate.
             * @return Daily rate = annual_rate / trading_days_per_year.
             */
            double to_daily_rate(double annual_rate) const;

            // ---------------------------------------------------------------
            // Member variables
            // ---------------------------------------------------------------

            std::vector<double> portfolio_returns_;
            std::vector<double> benchmark_returns_;
            std::vector<double> excess_returns_;
            double risk_free_rate_;
            int trading_days_per_year_;

            RegressionResult regression_result_;
            CaptureRatios capture_ratios_;
            double tracking_error_;
            double information_ratio_;
            double active_return_;
        };

    } // namespace analytics
} // namespace portfolio

#endif // PORTFOLIO_ANALYTICS_BENCHMARK_ANALYSIS_HPP