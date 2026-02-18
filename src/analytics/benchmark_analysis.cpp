/**
 * @file benchmark_analysis.cpp
 * @brief Implementation of the BenchmarkAnalysis class.
 *
 * Performs CAPM regression via ordinary least squares and computes
 * relative performance metrics. The regression model is:
 *
 *   (R_p - R_f) = alpha + beta * (R_b - R_f) + epsilon
 *
 * Alpha is annualized by multiplying the daily intercept by
 * trading_days_per_year. All other annualization uses sqrt(252)
 * scaling where appropriate.
 */

#include "analytics/benchmark_analysis.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace portfolio
{
    namespace analytics
    {

        // ===================================================================
        // Constructors
        // ===================================================================

        BenchmarkAnalysis::BenchmarkAnalysis(const std::vector<double> &portfolio_returns,
                                             const std::vector<double> &benchmark_returns,
                                             double risk_free_rate,
                                             int trading_days_per_year)
            : portfolio_returns_(portfolio_returns), benchmark_returns_(benchmark_returns), risk_free_rate_(risk_free_rate), trading_days_per_year_(trading_days_per_year), tracking_error_(0.0), information_ratio_(0.0), active_return_(0.0)
        {
            if (portfolio_returns_.empty())
            {
                throw std::invalid_argument("Portfolio return series cannot be empty");
            }
            if (benchmark_returns_.empty())
            {
                throw std::invalid_argument("Benchmark return series cannot be empty");
            }
            if (portfolio_returns_.size() != benchmark_returns_.size())
            {
                throw std::invalid_argument(
                    "Portfolio return series size (" + std::to_string(portfolio_returns_.size()) + ") must match benchmark return series size (" + std::to_string(benchmark_returns_.size()) + ")");
            }
            if (portfolio_returns_.size() < 3)
            {
                throw std::invalid_argument(
                    "At least 3 observations are required for regression, got: " + std::to_string(portfolio_returns_.size()));
            }
            if (trading_days_per_year <= 0)
            {
                throw std::invalid_argument(
                    "Expected positive value for parameter 'trading_days_per_year', got: " + std::to_string(trading_days_per_year));
            }

            // Compute excess returns (portfolio - benchmark)
            int n = static_cast<int>(portfolio_returns_.size());
            excess_returns_.resize(n);
            for (int i = 0; i < n; ++i)
            {
                excess_returns_[i] = portfolio_returns_[i] - benchmark_returns_[i];
            }

            run_regression();
            compute_tracking_metrics();
            compute_capture_ratios();
        }

        // ===================================================================
        // CAPM Regression Metrics
        // ===================================================================

        double BenchmarkAnalysis::alpha() const
        {
            return regression_result_.alpha;
        }

        double BenchmarkAnalysis::beta() const
        {
            return regression_result_.beta;
        }

        double BenchmarkAnalysis::r_squared() const
        {
            return regression_result_.r_squared;
        }

        const RegressionResult &BenchmarkAnalysis::regression() const
        {
            return regression_result_;
        }

        // ===================================================================
        // Tracking and Relative Risk
        // ===================================================================

        double BenchmarkAnalysis::tracking_error() const
        {
            return tracking_error_;
        }

        double BenchmarkAnalysis::information_ratio() const
        {
            return information_ratio_;
        }

        double BenchmarkAnalysis::active_return() const
        {
            return active_return_;
        }

        // ===================================================================
        // Capture Ratios
        // ===================================================================

        double BenchmarkAnalysis::up_capture_ratio() const
        {
            return capture_ratios_.up_capture;
        }

        double BenchmarkAnalysis::down_capture_ratio() const
        {
            return capture_ratios_.down_capture;
        }

        const CaptureRatios &BenchmarkAnalysis::capture_ratios() const
        {
            return capture_ratios_;
        }

        // ===================================================================
        // Excess Return Series
        // ===================================================================

        const std::vector<double> &BenchmarkAnalysis::excess_returns() const
        {
            return excess_returns_;
        }

        // ===================================================================
        // Export
        // ===================================================================

        std::string BenchmarkAnalysis::summary() const
        {
            std::ostringstream oss;
            oss << std::fixed;

            oss << "Benchmark Analysis Summary\n";
            oss << "==========================\n\n";

            oss << "CAPM Regression:\n";
            oss << "  Alpha (annualized):  " << std::setprecision(4)
                << regression_result_.alpha * 100.0 << "%\n";
            oss << "  Beta:                " << std::setprecision(4)
                << regression_result_.beta << "\n";
            oss << "  R-squared:           " << std::setprecision(4)
                << regression_result_.r_squared << "\n";
            oss << "  Observations:        " << regression_result_.num_observations << "\n";
            oss << "\n";

            oss << "Relative Risk:\n";
            oss << "  Tracking Error:      " << std::setprecision(4)
                << tracking_error_ * 100.0 << "%\n";
            oss << "  Information Ratio:   " << std::setprecision(4)
                << information_ratio_ << "\n";
            oss << "  Active Return:       " << std::setprecision(4)
                << active_return_ * 100.0 << "%\n";
            oss << "\n";

            oss << "Capture Ratios:\n";
            oss << "  Up Capture:          " << std::setprecision(4)
                << capture_ratios_.up_capture * 100.0 << "%\n";
            oss << "  Down Capture:        " << std::setprecision(4)
                << capture_ratios_.down_capture * 100.0 << "%\n";
            oss << "  Up Periods:          " << capture_ratios_.up_periods << "\n";
            oss << "  Down Periods:        " << capture_ratios_.down_periods << "\n";

            return oss.str();
        }

        std::string BenchmarkAnalysis::to_json() const
        {
            nlohmann::json j;

            j["regression"]["alpha"] = regression_result_.alpha;
            j["regression"]["alpha_daily"] = regression_result_.alpha_daily;
            j["regression"]["beta"] = regression_result_.beta;
            j["regression"]["r_squared"] = regression_result_.r_squared;
            j["regression"]["residual_std"] = regression_result_.residual_std;
            j["regression"]["num_observations"] = regression_result_.num_observations;

            j["relative_risk"]["tracking_error"] = tracking_error_;
            j["relative_risk"]["information_ratio"] = information_ratio_;
            j["relative_risk"]["active_return"] = active_return_;

            j["capture"]["up_capture"] = capture_ratios_.up_capture;
            j["capture"]["down_capture"] = capture_ratios_.down_capture;
            j["capture"]["up_periods"] = capture_ratios_.up_periods;
            j["capture"]["down_periods"] = capture_ratios_.down_periods;

            j["settings"]["risk_free_rate"] = risk_free_rate_;
            j["settings"]["trading_days_per_year"] = trading_days_per_year_;

            return j.dump(2);
        }

        // ===================================================================
        // Private Helpers
        // ===================================================================

        void BenchmarkAnalysis::run_regression()
        {
            // OLS regression: y = alpha + beta * x
            // where y = (R_p - R_f) and x = (R_b - R_f)

            int n = static_cast<int>(portfolio_returns_.size());
            double daily_rf = to_daily_rate(risk_free_rate_);

            // Compute excess returns over risk-free rate
            double sum_x = 0.0;
            double sum_y = 0.0;
            double sum_xx = 0.0;
            double sum_xy = 0.0;

            for (int i = 0; i < n; ++i)
            {
                double x = benchmark_returns_[i] - daily_rf;
                double y = portfolio_returns_[i] - daily_rf;
                sum_x += x;
                sum_y += y;
                sum_xx += x * x;
                sum_xy += x * y;
            }

            double nd = static_cast<double>(n);
            double mean_x = sum_x / nd;
            double mean_y = sum_y / nd;

            // Beta = Cov(x,y) / Var(x)
            double cov_xy = sum_xy / nd - mean_x * mean_y;
            double var_x = sum_xx / nd - mean_x * mean_x;

            if (std::abs(var_x) < 1e-18)
            {
                // Benchmark has zero variance - beta is undefined
                regression_result_.beta = 0.0;
                regression_result_.alpha_daily = mean_y;
            }
            else
            {
                regression_result_.beta = cov_xy / var_x;
                regression_result_.alpha_daily = mean_y - regression_result_.beta * mean_x;
            }

            // Annualize alpha (simple scaling)
            regression_result_.alpha = regression_result_.alpha_daily * static_cast<double>(trading_days_per_year_);

            // Compute R-squared and residual standard deviation
            double ss_res = 0.0;
            double ss_tot = 0.0;

            for (int i = 0; i < n; ++i)
            {
                double x = benchmark_returns_[i] - daily_rf;
                double y = portfolio_returns_[i] - daily_rf;
                double y_hat = regression_result_.alpha_daily + regression_result_.beta * x;
                double residual = y - y_hat;
                ss_res += residual * residual;
                ss_tot += (y - mean_y) * (y - mean_y);
            }

            regression_result_.r_squared = (std::abs(ss_tot) < 1e-18)
                                               ? 0.0
                                               : (1.0 - ss_res / ss_tot);

            // Residual standard deviation (using n-2 degrees of freedom)
            regression_result_.residual_std = (n > 2)
                                                  ? std::sqrt(ss_res / static_cast<double>(n - 2))
                                                  : 0.0;

            regression_result_.num_observations = n;
        }

        void BenchmarkAnalysis::compute_tracking_metrics()
        {
            int n = static_cast<int>(excess_returns_.size());

            // Active return: annualized mean of excess returns
            double mean_excess = std::accumulate(excess_returns_.begin(),
                                                 excess_returns_.end(), 0.0) /
                                 static_cast<double>(n);
            active_return_ = mean_excess * static_cast<double>(trading_days_per_year_);

            // Tracking error: annualized std dev of excess returns
            double sum_sq = 0.0;
            for (double er : excess_returns_)
            {
                double diff = er - mean_excess;
                sum_sq += diff * diff;
            }
            double daily_te = (n > 1)
                                  ? std::sqrt(sum_sq / static_cast<double>(n - 1))
                                  : 0.0;
            tracking_error_ = daily_te * std::sqrt(static_cast<double>(trading_days_per_year_));

            // Information ratio: active return / tracking error
            information_ratio_ = (tracking_error_ > 1e-18)
                                     ? (active_return_ / tracking_error_)
                                     : 0.0;
        }

        void BenchmarkAnalysis::compute_capture_ratios()
        {
            int n = static_cast<int>(benchmark_returns_.size());

            double up_port_sum = 0.0;
            double up_bench_sum = 0.0;
            int up_count = 0;

            double down_port_sum = 0.0;
            double down_bench_sum = 0.0;
            int down_count = 0;

            for (int i = 0; i < n; ++i)
            {
                if (benchmark_returns_[i] > 0.0)
                {
                    up_port_sum += portfolio_returns_[i];
                    up_bench_sum += benchmark_returns_[i];
                    ++up_count;
                }
                else
                {
                    down_port_sum += portfolio_returns_[i];
                    down_bench_sum += benchmark_returns_[i];
                    ++down_count;
                }
            }

            capture_ratios_.up_periods = up_count;
            capture_ratios_.down_periods = down_count;

            // Up capture: mean portfolio return in up markets / mean benchmark return in up markets
            if (up_count > 0 && std::abs(up_bench_sum) > 1e-18)
            {
                capture_ratios_.up_capture = (up_port_sum / static_cast<double>(up_count)) / (up_bench_sum / static_cast<double>(up_count));
            }
            else
            {
                capture_ratios_.up_capture = 0.0;
            }

            // Down capture: mean portfolio return in down markets / mean benchmark return in down markets
            if (down_count > 0 && std::abs(down_bench_sum) > 1e-18)
            {
                capture_ratios_.down_capture = (down_port_sum / static_cast<double>(down_count)) / (down_bench_sum / static_cast<double>(down_count));
            }
            else
            {
                capture_ratios_.down_capture = 0.0;
            }
        }

        double BenchmarkAnalysis::to_daily_rate(double annual_rate) const
        {
            return annual_rate / static_cast<double>(trading_days_per_year_);
        }

    } // namespace analytics
} // namespace portfolio