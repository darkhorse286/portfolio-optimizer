/**
 * @file performance_metrics.cpp
 * @brief Implementation of the PerformanceMetrics class.
 *
 * Computes return, risk, and risk-adjusted performance metrics from
 * NAV and return time series. All annualized values use simple scaling
 * (multiply/divide by trading_days_per_year or its square root) unless
 * otherwise noted (CAGR uses geometric compounding).
 */

#include "analytics/performance_metrics.hpp"
#include "backtest/backtest_result.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace portfolio
{
    namespace analytics
    {

        // ===================================================================
        // Anonymous namespace: math helpers for parametric VaR/CVaR
        // ===================================================================

        namespace
        {

            /**
             * @brief Standard normal PDF: phi(x) = (1/sqrt(2*pi)) * exp(-x^2/2).
             */
            double normal_pdf(double x)
            {
                static const double INV_SQRT_2PI = 0.3989422804014327;
                return INV_SQRT_2PI * std::exp(-0.5 * x * x);
            }

            /**
             * @brief Rational approximation of the inverse normal CDF (probit).
             *
             * Uses the Beasley-Springer-Moro algorithm, accurate to approximately
             * 1e-9 for p in [1e-8, 1 - 1e-8].
             *
             * @param p Probability in (0, 1).
             * @return z such that Phi(z) = p.
             */
            double inverse_normal_cdf(double p)
            {
                static const double a[] = {
                    -3.969683028665376e+01, 2.209460984245205e+02,
                    -2.759285104469687e+02, 1.383577518672690e+02,
                    -3.066479806614716e+01, 2.506628277459239e+00};
                static const double b[] = {
                    -5.447609879822406e+01, 1.615858368580409e+02,
                    -1.556989798598866e+02, 6.680131188771972e+01,
                    -1.328068155288572e+01};
                static const double c[] = {
                    -7.784894002430293e-03, -3.223964580411365e-01,
                    -2.400758277161838e+00, -2.549732539343734e+00,
                    4.374664141464968e+00, 2.938163982698783e+00};
                static const double d[] = {
                    7.784695709041462e-03, 3.224671290700398e-01,
                    2.445134137142996e+00, 3.754408661907416e+00};

                static const double P_LOW = 0.02425;
                static const double P_HIGH = 1.0 - P_LOW;

                double q;
                double r;

                if (p < P_LOW)
                {
                    q = std::sqrt(-2.0 * std::log(p));
                    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
                }
                else if (p <= P_HIGH)
                {
                    q = p - 0.5;
                    r = q * q;
                    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
                }
                else
                {
                    q = std::sqrt(-2.0 * std::log(1.0 - p));
                    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
                }
            }

        } // anonymous namespace

        // ===================================================================
        // Constructors
        // ===================================================================

        PerformanceMetrics::PerformanceMetrics(const backtest::BacktestResult &result,
                                               double risk_free_rate,
                                               int trading_days_per_year)
            : nav_series_(result.nav_series), return_series_(result.return_series), dates_(result.dates), risk_free_rate_(risk_free_rate), trading_days_per_year_(trading_days_per_year), drawdown_computed_(false)
        {
            initialize();
        }

        PerformanceMetrics::PerformanceMetrics(const std::vector<double> &nav_series,
                                               const std::vector<double> &return_series,
                                               const std::vector<std::string> &dates,
                                               double risk_free_rate,
                                               int trading_days_per_year)
            : nav_series_(nav_series), return_series_(return_series), dates_(dates), risk_free_rate_(risk_free_rate), trading_days_per_year_(trading_days_per_year), drawdown_computed_(false)
        {
            initialize();
        }

        // ===================================================================
        // Return Metrics
        // ===================================================================

        double PerformanceMetrics::total_return() const
        {
            return (nav_series_.back() / nav_series_.front()) - 1.0;
        }

        double PerformanceMetrics::annualized_return() const
        {
            double total = total_return();
            double years = static_cast<double>(nav_series_.size() - 1) / static_cast<double>(trading_days_per_year_);
            if (years <= 0.0)
            {
                return 0.0;
            }
            // CAGR = (1 + total_return)^(1/years) - 1
            return std::pow(1.0 + total, 1.0 / years) - 1.0;
        }

        std::vector<double> PerformanceMetrics::rolling_returns(int window_days) const
        {
            if (window_days < 1)
            {
                throw std::invalid_argument(
                    "Expected positive value for parameter 'window_days', got: " + std::to_string(window_days));
            }
            int n = static_cast<int>(nav_series_.size());
            if (window_days >= n)
            {
                throw std::invalid_argument(
                    "Rolling window (" + std::to_string(window_days) + ") must be less than series length (" + std::to_string(n) + ")");
            }

            std::vector<double> result(n - window_days);
            for (int i = 0; i < static_cast<int>(result.size()); ++i)
            {
                result[i] = (nav_series_[i + window_days] / nav_series_[i]) - 1.0;
            }
            return result;
        }

        MonthlyReturnTable PerformanceMetrics::monthly_return_table() const
        {
            MonthlyReturnTable table;

            struct MonthKey
            {
                int year;
                int month;
            };

            std::vector<MonthKey> keys;
            keys.reserve(dates_.size());

            for (const auto &date_str : dates_)
            {
                if (date_str.size() < 10)
                {
                    throw std::runtime_error(
                        "Date string too short for YYYY-MM-DD format: '" + date_str + "'");
                }
                int year = std::stoi(date_str.substr(0, 4));
                int month = std::stoi(date_str.substr(5, 2));
                keys.push_back({year, month});
            }

            // For each month, compute return from first NAV to last NAV in that month
            int i = 0;
            int n = static_cast<int>(nav_series_.size());

            while (i < n)
            {
                int year = keys[i].year;
                int month = keys[i].month;
                double start_nav = nav_series_[i];

                // Find last index in this month
                int j = i;
                while (j + 1 < n && keys[j + 1].year == year && keys[j + 1].month == month)
                {
                    ++j;
                }
                double end_nav = nav_series_[j];
                double monthly_ret = (end_nav / start_nav) - 1.0;
                table.monthly_returns[year][month] = monthly_ret;

                i = j + 1;
            }

            // Compute annual returns by compounding monthly returns
            for (const auto &[year, months] : table.monthly_returns)
            {
                double compounded = 1.0;
                for (const auto &[mon, ret] : months)
                {
                    compounded *= (1.0 + ret);
                }
                table.annual_returns[year] = compounded - 1.0;
            }

            return table;
        }

        // ===================================================================
        // Risk Metrics
        // ===================================================================

        double PerformanceMetrics::annualized_volatility() const
        {
            int n = static_cast<int>(return_series_.size());
            if (n < 2)
            {
                return 0.0;
            }

            double mean = std::accumulate(return_series_.begin(), return_series_.end(), 0.0) / static_cast<double>(n);

            double sum_sq = 0.0;
            for (double r : return_series_)
            {
                double diff = r - mean;
                sum_sq += diff * diff;
            }
            double daily_vol = std::sqrt(sum_sq / static_cast<double>(n - 1));
            return daily_vol * std::sqrt(static_cast<double>(trading_days_per_year_));
        }

        double PerformanceMetrics::downside_deviation(double target_return) const
        {
            int n = static_cast<int>(return_series_.size());
            if (n < 2)
            {
                return 0.0;
            }

            double daily_target = to_daily_rate(target_return);
            double sum_sq = 0.0;

            for (double r : return_series_)
            {
                double diff = r - daily_target;
                if (diff < 0.0)
                {
                    sum_sq += diff * diff;
                }
            }

            // Use full count in denominator (standard semi-deviation formula)
            double daily_dd = std::sqrt(sum_sq / static_cast<double>(n - 1));
            return daily_dd * std::sqrt(static_cast<double>(trading_days_per_year_));
        }

        double PerformanceMetrics::max_drawdown() const
        {
            compute_drawdown_cache();
            return max_drawdown_info_cache_.depth;
        }

        DrawdownInfo PerformanceMetrics::max_drawdown_info() const
        {
            compute_drawdown_cache();
            return max_drawdown_info_cache_;
        }

        std::vector<double> PerformanceMetrics::drawdown_series() const
        {
            compute_drawdown_cache();
            return drawdown_series_cache_;
        }

        double PerformanceMetrics::value_at_risk(double confidence,
                                                 VaRMethod method) const
        {
            if (confidence <= 0.0 || confidence >= 1.0)
            {
                throw std::invalid_argument(
                    "Confidence level must be in (0, 1), got: " + std::to_string(confidence));
            }

            int n = static_cast<int>(return_series_.size());

            if (method == VaRMethod::PARAMETRIC)
            {
                double mean = std::accumulate(return_series_.begin(),
                                              return_series_.end(), 0.0) /
                              static_cast<double>(n);
                double sum_sq = 0.0;
                for (double r : return_series_)
                {
                    double diff = r - mean;
                    sum_sq += diff * diff;
                }
                double std_dev = std::sqrt(sum_sq / static_cast<double>(n - 1));

                double z = inverse_normal_cdf(1.0 - confidence);
                return -(mean + z * std_dev);
            }

            // Historical VaR: sort returns, find the (1-confidence) quantile
            std::vector<double> sorted_returns(return_series_);
            std::sort(sorted_returns.begin(), sorted_returns.end());

            double index = (1.0 - confidence) * static_cast<double>(n - 1);
            int lower = static_cast<int>(std::floor(index));
            int upper = static_cast<int>(std::ceil(index));

            double var_value;
            if (lower == upper || upper >= n)
            {
                var_value = sorted_returns[lower];
            }
            else
            {
                double frac = index - static_cast<double>(lower);
                var_value = sorted_returns[lower] * (1.0 - frac) + sorted_returns[upper] * frac;
            }

            return -var_value;
        }

        double PerformanceMetrics::conditional_var(double confidence,
                                                   VaRMethod method) const
        {
            if (confidence <= 0.0 || confidence >= 1.0)
            {
                throw std::invalid_argument(
                    "Confidence level must be in (0, 1), got: " + std::to_string(confidence));
            }

            int n = static_cast<int>(return_series_.size());

            if (method == VaRMethod::PARAMETRIC)
            {
                double mean = std::accumulate(return_series_.begin(),
                                              return_series_.end(), 0.0) /
                              static_cast<double>(n);
                double sum_sq = 0.0;
                for (double r : return_series_)
                {
                    double diff = r - mean;
                    sum_sq += diff * diff;
                }
                double std_dev = std::sqrt(sum_sq / static_cast<double>(n - 1));

                double z = inverse_normal_cdf(1.0 - confidence);
                double phi_z = normal_pdf(z);
                return -(mean + std_dev * phi_z / (1.0 - confidence));
            }

            // Historical CVaR: average of returns at or below the VaR threshold
            std::vector<double> sorted_returns(return_series_);
            std::sort(sorted_returns.begin(), sorted_returns.end());

            int cutoff = static_cast<int>(std::floor((1.0 - confidence) * static_cast<double>(n)));
            if (cutoff < 1)
            {
                cutoff = 1;
            }

            double sum = 0.0;
            for (int i = 0; i < cutoff; ++i)
            {
                sum += sorted_returns[i];
            }
            return -(sum / static_cast<double>(cutoff));
        }

        double PerformanceMetrics::skewness() const
        {
            int n = static_cast<int>(return_series_.size());
            if (n < 3)
            {
                return 0.0;
            }

            double mean = std::accumulate(return_series_.begin(),
                                          return_series_.end(), 0.0) /
                          static_cast<double>(n);

            double m2 = 0.0;
            double m3 = 0.0;
            for (double r : return_series_)
            {
                double diff = r - mean;
                double d2 = diff * diff;
                m2 += d2;
                m3 += d2 * diff;
            }
            m2 /= static_cast<double>(n);
            m3 /= static_cast<double>(n);

            if (m2 < 1e-18)
            {
                return 0.0;
            }

            double raw_skew = m3 / std::pow(m2, 1.5);

            // Bias correction (Fisher's definition, matches scipy.stats.skew)
            double nd = static_cast<double>(n);
            double adjust = std::sqrt(nd * (nd - 1.0)) / (nd - 2.0);
            return adjust * raw_skew;
        }

        double PerformanceMetrics::kurtosis() const
        {
            int n = static_cast<int>(return_series_.size());
            if (n < 4)
            {
                return 0.0;
            }

            double mean = std::accumulate(return_series_.begin(),
                                          return_series_.end(), 0.0) /
                          static_cast<double>(n);

            double m2 = 0.0;
            double m4 = 0.0;
            for (double r : return_series_)
            {
                double diff = r - mean;
                double d2 = diff * diff;
                m2 += d2;
                m4 += d2 * d2;
            }
            m2 /= static_cast<double>(n);
            m4 /= static_cast<double>(n);

            if (m2 < 1e-18)
            {
                return 0.0;
            }

            double raw_kurt = m4 / (m2 * m2);

            // Bias-corrected excess kurtosis (matches scipy.stats.kurtosis)
            double nd = static_cast<double>(n);
            double excess = ((nd + 1.0) * raw_kurt - 3.0 * (nd - 1.0)) * (nd - 1.0) / ((nd - 2.0) * (nd - 3.0));
            return excess;
        }

        // ===================================================================
        // Risk-Adjusted Metrics
        // ===================================================================

        double PerformanceMetrics::sharpe_ratio() const
        {
            double vol = annualized_volatility();
            if (vol < 1e-18)
            {
                return 0.0;
            }
            return (annualized_return() - risk_free_rate_) / vol;
        }

        double PerformanceMetrics::sortino_ratio(double target_return) const
        {
            double dd = downside_deviation(target_return);
            if (dd < 1e-18)
            {
                return 0.0;
            }
            return (annualized_return() - target_return) / dd;
        }

        double PerformanceMetrics::calmar_ratio() const
        {
            double mdd = max_drawdown();
            if (mdd < 1e-18)
            {
                return 0.0;
            }
            return annualized_return() / mdd;
        }

        double PerformanceMetrics::omega_ratio(double threshold) const
        {
            double daily_threshold = to_daily_rate(threshold);

            double gains = 0.0;
            double losses = 0.0;

            for (double r : return_series_)
            {
                double excess = r - daily_threshold;
                if (excess > 0.0)
                {
                    gains += excess;
                }
                else
                {
                    losses -= excess;
                }
            }

            if (losses < 1e-18 && gains < 1e-18)
            {
                throw std::invalid_argument(
                    "All returns equal the threshold; Omega ratio is undefined");
            }

            if (losses < 1e-18)
            {
                return std::numeric_limits<double>::max();
            }

            return gains / losses;
        }

        // ===================================================================
        // Accessors
        // ===================================================================

        const std::vector<double> &PerformanceMetrics::nav_series() const
        {
            return nav_series_;
        }

        const std::vector<double> &PerformanceMetrics::return_series() const
        {
            return return_series_;
        }

        const std::vector<std::string> &PerformanceMetrics::dates() const
        {
            return dates_;
        }

        double PerformanceMetrics::risk_free_rate() const
        {
            return risk_free_rate_;
        }

        int PerformanceMetrics::trading_days_per_year() const
        {
            return trading_days_per_year_;
        }

        // ===================================================================
        // Export
        // ===================================================================

        std::string PerformanceMetrics::summary() const
        {
            std::ostringstream oss;
            oss << std::fixed;

            oss << "Performance Summary\n";
            oss << "===================\n";
            oss << "\n";

            oss << "Return Metrics:\n";
            oss << "  Total Return:        " << std::setprecision(4)
                << total_return() * 100.0 << "%\n";
            oss << "  Annualized Return:   " << std::setprecision(4)
                << annualized_return() * 100.0 << "%\n";
            oss << "\n";

            oss << "Risk Metrics:\n";
            oss << "  Annualized Vol:      " << std::setprecision(4)
                << annualized_volatility() * 100.0 << "%\n";
            oss << "  Downside Deviation:  " << std::setprecision(4)
                << downside_deviation() * 100.0 << "%\n";
            oss << "  Max Drawdown:        " << std::setprecision(4)
                << max_drawdown() * 100.0 << "%\n";
            oss << "  VaR (95%, Hist):     " << std::setprecision(4)
                << value_at_risk(0.95) * 100.0 << "%\n";
            oss << "  CVaR (95%, Hist):    " << std::setprecision(4)
                << conditional_var(0.95) * 100.0 << "%\n";
            oss << "  Skewness:            " << std::setprecision(4)
                << skewness() << "\n";
            oss << "  Kurtosis:            " << std::setprecision(4)
                << kurtosis() << "\n";
            oss << "\n";

            oss << "Risk-Adjusted Metrics:\n";
            oss << "  Sharpe Ratio:        " << std::setprecision(4)
                << sharpe_ratio() << "\n";
            oss << "  Sortino Ratio:       " << std::setprecision(4)
                << sortino_ratio() << "\n";
            oss << "  Calmar Ratio:        " << std::setprecision(4)
                << calmar_ratio() << "\n";
            oss << "  Omega Ratio:         " << std::setprecision(4)
                << omega_ratio() << "\n";
            oss << "\n";

            auto dd = max_drawdown_info();
            oss << "Max Drawdown Detail:\n";
            oss << "  Depth:               " << std::setprecision(4)
                << dd.depth * 100.0 << "%\n";
            oss << "  Duration (days):     " << dd.duration_days << "\n";
            oss << "  Recovery (days):     ";
            if (dd.recovery_days < 0)
            {
                oss << "Unrecovered\n";
            }
            else
            {
                oss << dd.recovery_days << "\n";
            }

            return oss.str();
        }

        std::string PerformanceMetrics::to_json() const
        {
            nlohmann::json j;

            j["return_metrics"]["total_return"] = total_return();
            j["return_metrics"]["annualized_return"] = annualized_return();

            j["risk_metrics"]["annualized_volatility"] = annualized_volatility();
            j["risk_metrics"]["downside_deviation"] = downside_deviation();
            j["risk_metrics"]["max_drawdown"] = max_drawdown();
            j["risk_metrics"]["var_95_historical"] = value_at_risk(0.95, VaRMethod::HISTORICAL);
            j["risk_metrics"]["var_95_parametric"] = value_at_risk(0.95, VaRMethod::PARAMETRIC);
            j["risk_metrics"]["cvar_95_historical"] = conditional_var(0.95, VaRMethod::HISTORICAL);
            j["risk_metrics"]["cvar_95_parametric"] = conditional_var(0.95, VaRMethod::PARAMETRIC);
            j["risk_metrics"]["skewness"] = skewness();
            j["risk_metrics"]["kurtosis"] = kurtosis();

            j["risk_adjusted"]["sharpe_ratio"] = sharpe_ratio();
            j["risk_adjusted"]["sortino_ratio"] = sortino_ratio();
            j["risk_adjusted"]["calmar_ratio"] = calmar_ratio();
            j["risk_adjusted"]["omega_ratio"] = omega_ratio();

            auto dd = max_drawdown_info();
            j["max_drawdown_detail"]["depth"] = dd.depth;
            j["max_drawdown_detail"]["duration_days"] = dd.duration_days;
            j["max_drawdown_detail"]["recovery_days"] = dd.recovery_days;
            j["max_drawdown_detail"]["peak_index"] = dd.peak_index;
            j["max_drawdown_detail"]["trough_index"] = dd.trough_index;

            j["settings"]["risk_free_rate"] = risk_free_rate_;
            j["settings"]["trading_days_per_year"] = trading_days_per_year_;
            j["settings"]["num_observations"] = static_cast<int>(return_series_.size());

            return j.dump(2);
        }

        void PerformanceMetrics::to_csv(const std::string &filepath) const
        {
            std::ofstream file(filepath);
            if (!file.is_open())
            {
                throw std::runtime_error("Cannot open file for writing: " + filepath);
            }

            compute_drawdown_cache();

            file << "date,nav,return,drawdown\n";
            file << std::fixed << std::setprecision(8);

            // First row: NAV and drawdown only (no return for day 0)
            file << dates_[0] << "," << nav_series_[0] << ",,"
                 << drawdown_series_cache_[0] << "\n";

            for (size_t i = 1; i < nav_series_.size(); ++i)
            {
                file << dates_[i] << ","
                     << nav_series_[i] << ",";
                if (i - 1 < return_series_.size())
                {
                    file << return_series_[i - 1];
                }
                file << "," << drawdown_series_cache_[i] << "\n";
            }
        }

        // ===================================================================
        // Private Helpers
        // ===================================================================

        void PerformanceMetrics::initialize()
        {
            if (nav_series_.size() < 2)
            {
                throw std::invalid_argument(
                    "NAV series must have at least 2 elements, got: " + std::to_string(nav_series_.size()));
            }
            if (return_series_.empty())
            {
                throw std::invalid_argument("Return series cannot be empty");
            }

            // Return series can be nav_series.size()-1 or nav_series.size()
            size_t expected_min = nav_series_.size() - 1;
            if (return_series_.size() != expected_min && return_series_.size() != nav_series_.size())
            {
                throw std::invalid_argument(
                    "Return series size (" + std::to_string(return_series_.size()) + ") is inconsistent with NAV series size (" + std::to_string(nav_series_.size()) + ")");
            }

            if (dates_.size() != nav_series_.size())
            {
                throw std::invalid_argument(
                    "Dates size (" + std::to_string(dates_.size()) + ") must match NAV series size (" + std::to_string(nav_series_.size()) + ")");
            }

            if (trading_days_per_year_ <= 0)
            {
                throw std::invalid_argument(
                    "Expected positive value for parameter 'trading_days_per_year', got: " + std::to_string(trading_days_per_year_));
            }
        }

        void PerformanceMetrics::compute_drawdown_cache() const
        {
            if (drawdown_computed_)
            {
                return;
            }

            int n = static_cast<int>(nav_series_.size());
            drawdown_series_cache_.resize(n);

            double peak = nav_series_[0];
            double max_dd = 0.0;
            int peak_idx = 0;
            int best_peak_idx = 0;
            int best_trough_idx = 0;

            for (int i = 0; i < n; ++i)
            {
                if (nav_series_[i] > peak)
                {
                    peak = nav_series_[i];
                    peak_idx = i;
                }
                double dd = (nav_series_[i] - peak) / peak; // Non-positive
                drawdown_series_cache_[i] = dd;

                if (-dd > max_dd)
                {
                    max_dd = -dd;
                    best_peak_idx = peak_idx;
                    best_trough_idx = i;
                }
            }

            // Find recovery point after the worst trough
            int recovery_idx = -1;
            double peak_nav = nav_series_[best_peak_idx];
            for (int i = best_trough_idx + 1; i < n; ++i)
            {
                if (nav_series_[i] >= peak_nav)
                {
                    recovery_idx = i;
                    break;
                }
            }

            max_drawdown_info_cache_.depth = max_dd;
            max_drawdown_info_cache_.duration_days = best_trough_idx - best_peak_idx;
            max_drawdown_info_cache_.recovery_days =
                (recovery_idx >= 0) ? (recovery_idx - best_trough_idx) : -1;
            max_drawdown_info_cache_.peak_index = best_peak_idx;
            max_drawdown_info_cache_.trough_index = best_trough_idx;
            max_drawdown_info_cache_.recovery_index = recovery_idx;

            drawdown_computed_ = true;
        }

        double PerformanceMetrics::to_daily_rate(double annual_rate) const
        {
            return annual_rate / static_cast<double>(trading_days_per_year_);
        }

    } // namespace analytics
} // namespace portfolio
