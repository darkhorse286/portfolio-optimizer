/**
 * @file rolling_statistics.cpp
 * @brief Implementation of the RollingStatistics class.
 *
 * Computes rolling window metrics over return series. Uses direct
 * per-window computation for correctness and clarity. Output vectors
 * are aligned to the end of each window: result[i] corresponds to
 * the window [i, i + window - 1] of the input series.
 */

#include "analytics/rolling_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace portfolio
{
    namespace analytics
    {

        // ===================================================================
        // Constructors
        // ===================================================================

        RollingStatistics::RollingStatistics(const std::vector<double> &return_series,
                                             const RollingConfig &config)
            : return_series_(return_series), config_(config)
        {
            if (return_series_.empty())
            {
                throw std::invalid_argument("Return series cannot be empty");
            }
            if (config_.window_days < 2)
            {
                throw std::invalid_argument(
                    "Expected window_days >= 2 for rolling statistics, got: " + std::to_string(config_.window_days));
            }
            if (static_cast<int>(return_series_.size()) < config_.window_days)
            {
                throw std::invalid_argument(
                    "Return series length (" + std::to_string(return_series_.size()) + ") must be >= window_days (" + std::to_string(config_.window_days) + ")");
            }
            if (config_.min_periods <= 0)
            {
                config_.min_periods = config_.window_days;
            }
        }

        RollingStatistics::RollingStatistics(const std::vector<double> &return_series,
                                             int window_days)
            : return_series_(return_series), config_(window_days)
        {
            if (return_series_.empty())
            {
                throw std::invalid_argument("Return series cannot be empty");
            }
            if (window_days < 2)
            {
                throw std::invalid_argument(
                    "Expected window_days >= 2 for rolling statistics, got: " + std::to_string(window_days));
            }
            if (static_cast<int>(return_series_.size()) < window_days)
            {
                throw std::invalid_argument(
                    "Return series length (" + std::to_string(return_series_.size()) + ") must be >= window_days (" + std::to_string(window_days) + ")");
            }
        }

        // ===================================================================
        // Rolling Metrics
        // ===================================================================

        std::vector<double> RollingStatistics::volatility() const
        {
            int w = config_.window_days;
            int tdy = config_.trading_days_per_year;

            return rolling_apply_internal(
                [w, tdy](const double *data, int size) -> double
                {
                    double nd = static_cast<double>(size);
                    double sum = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                    }
                    double mean = sum / nd;

                    double sum_sq = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        double diff = data[i] - mean;
                        sum_sq += diff * diff;
                    }
                    double daily_vol = std::sqrt(sum_sq / (nd - 1.0));
                    return daily_vol * std::sqrt(static_cast<double>(tdy));
                });
        }

        std::vector<double> RollingStatistics::mean_return() const
        {
            int tdy = config_.trading_days_per_year;

            return rolling_apply_internal(
                [tdy](const double *data, int size) -> double
                {
                    double sum = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                    }
                    return (sum / static_cast<double>(size)) * static_cast<double>(tdy);
                });
        }

        std::vector<double> RollingStatistics::sharpe_ratio() const
        {
            int tdy = config_.trading_days_per_year;
            double rf = config_.risk_free_rate;

            return rolling_apply_internal(
                [tdy, rf](const double *data, int size) -> double
                {
                    double nd = static_cast<double>(size);
                    double sum = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                    }
                    double mean = sum / nd;

                    double sum_sq = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        double diff = data[i] - mean;
                        sum_sq += diff * diff;
                    }
                    double daily_vol = std::sqrt(sum_sq / (nd - 1.0));

                    if (daily_vol < 1e-18)
                    {
                        return 0.0;
                    }

                    double ann_ret = mean * static_cast<double>(tdy);
                    double ann_vol = daily_vol * std::sqrt(static_cast<double>(tdy));
                    return (ann_ret - rf) / ann_vol;
                });
        }

        std::vector<double> RollingStatistics::sortino_ratio(double target_return) const
        {
            int tdy = config_.trading_days_per_year;
            double daily_target = target_return / static_cast<double>(tdy);

            return rolling_apply_internal(
                [tdy, target_return, daily_target](const double *data, int size) -> double
                {
                    double nd = static_cast<double>(size);
                    double sum = 0.0;
                    double sum_sq_down = 0.0;

                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                        double diff = data[i] - daily_target;
                        if (diff < 0.0)
                        {
                            sum_sq_down += diff * diff;
                        }
                    }

                    double daily_dd = std::sqrt(sum_sq_down / (nd - 1.0));
                    double ann_dd = daily_dd * std::sqrt(static_cast<double>(tdy));

                    if (ann_dd < 1e-18)
                    {
                        return 0.0;
                    }

                    double ann_ret = (sum / nd) * static_cast<double>(tdy);
                    return (ann_ret - target_return) / ann_dd;
                });
        }

        std::vector<double> RollingStatistics::max_drawdown() const
        {
            // For rolling max drawdown, we need to compute the max drawdown
            // within each window. This requires building a NAV-like curve
            // from the returns in each window.
            return rolling_apply_internal(
                [](const double *data, int size) -> double
                {
                    // Build cumulative wealth curve within the window
                    double nav = 1.0;
                    double peak = 1.0;
                    double max_dd = 0.0;

                    for (int i = 0; i < size; ++i)
                    {
                        nav *= (1.0 + data[i]);
                        if (nav > peak)
                        {
                            peak = nav;
                        }
                        double dd = (peak - nav) / peak;
                        if (dd > max_dd)
                        {
                            max_dd = dd;
                        }
                    }
                    return max_dd;
                });
        }

        std::vector<double> RollingStatistics::beta(
            const std::vector<double> &benchmark_returns) const
        {
            if (benchmark_returns.size() != return_series_.size())
            {
                throw std::invalid_argument(
                    "Benchmark return series size (" + std::to_string(benchmark_returns.size()) + ") must match portfolio return series size (" + std::to_string(return_series_.size()) + ")");
            }

            int n = static_cast<int>(return_series_.size());
            int w = config_.window_days;
            int out_size = n - w + 1;

            std::vector<double> result(out_size);

            for (int i = 0; i < out_size; ++i)
            {
                const double *port = return_series_.data() + i;
                const double *bench = benchmark_returns.data() + i;

                double sum_x = 0.0;
                double sum_y = 0.0;
                double sum_xx = 0.0;
                double sum_xy = 0.0;

                for (int j = 0; j < w; ++j)
                {
                    sum_x += bench[j];
                    sum_y += port[j];
                    sum_xx += bench[j] * bench[j];
                    sum_xy += bench[j] * port[j];
                }

                double wd = static_cast<double>(w);
                double mean_x = sum_x / wd;
                double var_x = sum_xx / wd - mean_x * mean_x;

                if (std::abs(var_x) < 1e-18)
                {
                    result[i] = 0.0;
                }
                else
                {
                    double cov_xy = sum_xy / wd - mean_x * (sum_y / wd);
                    result[i] = cov_xy / var_x;
                }
            }

            return result;
        }

        std::vector<double> RollingStatistics::tracking_error(
            const std::vector<double> &benchmark_returns) const
        {
            if (benchmark_returns.size() != return_series_.size())
            {
                throw std::invalid_argument(
                    "Benchmark return series size (" + std::to_string(benchmark_returns.size()) + ") must match portfolio return series size (" + std::to_string(return_series_.size()) + ")");
            }

            int n = static_cast<int>(return_series_.size());
            int w = config_.window_days;
            int tdy = config_.trading_days_per_year;
            int out_size = n - w + 1;

            std::vector<double> result(out_size);

            for (int i = 0; i < out_size; ++i)
            {
                double sum = 0.0;
                double sum_sq = 0.0;

                for (int j = 0; j < w; ++j)
                {
                    double excess = return_series_[i + j] - benchmark_returns[i + j];
                    sum += excess;
                    sum_sq += excess * excess;
                }

                double wd = static_cast<double>(w);
                double mean = sum / wd;
                double var = (sum_sq - wd * mean * mean) / (wd - 1.0);
                double daily_te = std::sqrt(std::max(0.0, var));
                result[i] = daily_te * std::sqrt(static_cast<double>(tdy));
            }

            return result;
        }

        std::vector<double> RollingStatistics::skewness() const
        {
            return rolling_apply_internal(
                [](const double *data, int size) -> double
                {
                    if (size < 3)
                    {
                        return 0.0;
                    }

                    double nd = static_cast<double>(size);
                    double sum = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                    }
                    double mean = sum / nd;

                    double m2 = 0.0;
                    double m3 = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        double diff = data[i] - mean;
                        double d2 = diff * diff;
                        m2 += d2;
                        m3 += d2 * diff;
                    }
                    m2 /= nd;
                    m3 /= nd;

                    if (m2 < 1e-18)
                    {
                        return 0.0;
                    }

                    double raw = m3 / std::pow(m2, 1.5);
                    double adjust = std::sqrt(nd * (nd - 1.0)) / (nd - 2.0);
                    return adjust * raw;
                });
        }

        std::vector<double> RollingStatistics::kurtosis() const
        {
            return rolling_apply_internal(
                [](const double *data, int size) -> double
                {
                    if (size < 4)
                    {
                        return 0.0;
                    }

                    double nd = static_cast<double>(size);
                    double sum = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        sum += data[i];
                    }
                    double mean = sum / nd;

                    double m2 = 0.0;
                    double m4 = 0.0;
                    for (int i = 0; i < size; ++i)
                    {
                        double diff = data[i] - mean;
                        double d2 = diff * diff;
                        m2 += d2;
                        m4 += d2 * d2;
                    }
                    m2 /= nd;
                    m4 /= nd;

                    if (m2 < 1e-18)
                    {
                        return 0.0;
                    }

                    double raw = m4 / (m2 * m2);
                    double excess = ((nd + 1.0) * raw - 3.0 * (nd - 1.0)) * (nd - 1.0) / ((nd - 2.0) * (nd - 3.0));
                    return excess;
                });
        }

        // ===================================================================
        // Generic Rolling Application
        // ===================================================================

        std::vector<double> RollingStatistics::apply(
            const std::function<double(const std::vector<double> &)> &func) const
        {
            int n = static_cast<int>(return_series_.size());
            int w = config_.window_days;
            int out_size = n - w + 1;

            std::vector<double> result(out_size);

            for (int i = 0; i < out_size; ++i)
            {
                std::vector<double> window(return_series_.begin() + i,
                                           return_series_.begin() + i + w);
                result[i] = func(window);
            }

            return result;
        }

        // ===================================================================
        // Accessors
        // ===================================================================

        const RollingConfig &RollingStatistics::config() const
        {
            return config_;
        }

        int RollingStatistics::output_size() const
        {
            return static_cast<int>(return_series_.size()) - config_.window_days + 1;
        }

        std::vector<std::string> RollingStatistics::aligned_dates(
            const std::vector<std::string> &dates) const
        {
            if (dates.size() != return_series_.size())
            {
                throw std::invalid_argument(
                    "Dates size (" + std::to_string(dates.size()) + ") must match return series size (" + std::to_string(return_series_.size()) + ")");
            }

            int w = config_.window_days;
            // Return dates aligned to the end of each window
            return std::vector<std::string>(dates.begin() + w - 1, dates.end());
        }

        // ===================================================================
        // Private Helpers
        // ===================================================================

        std::vector<double> RollingStatistics::rolling_apply_internal(
            const std::function<double(const double *, int)> &func) const
        {
            int n = static_cast<int>(return_series_.size());
            int w = config_.window_days;
            int out_size = n - w + 1;

            std::vector<double> result(out_size);

            for (int i = 0; i < out_size; ++i)
            {
                result[i] = func(return_series_.data() + i, w);
            }

            return result;
        }

    } // namespace analytics
} // namespace portfolio
