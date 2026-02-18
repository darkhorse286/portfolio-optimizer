/**
 * @file attribution.cpp
 * @brief Implementation of the Attribution class.
 *
 * Implements Brinson-Fachler single-period performance attribution and
 * multi-period linking via the Carino smoothing method. The Carino method
 * applies logarithmic smoothing factors so that single-period additive
 * effects aggregate to match the total geometric active return.
 *
 * Brinson-Fachler decomposition:
 *   Allocation  = (w_p - w_b) * (R_b_sector - R_b_total)
 *   Selection   = w_b * (R_p_sector - R_b_sector)
 *   Interaction = (w_p - w_b) * (R_p_sector - R_b_sector)
 */

#include "analytics/attribution.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace portfolio
{
    namespace analytics
    {

        // ===================================================================
        // Single-Period Attribution
        // ===================================================================

        SinglePeriodAttribution Attribution::single_period(
            const std::vector<SectorAllocation> &portfolio_sectors,
            const std::vector<SectorAllocation> &benchmark_sectors)
        {
            if (portfolio_sectors.empty())
            {
                throw std::invalid_argument(
                    "Portfolio sector allocations cannot be empty");
            }
            if (benchmark_sectors.empty())
            {
                throw std::invalid_argument(
                    "Benchmark sector allocations cannot be empty");
            }

            // Build lookup maps for benchmark sectors
            std::unordered_map<std::string, const SectorAllocation *> bench_map;
            for (const auto &sector : benchmark_sectors)
            {
                bench_map[sector.sector_name] = &sector;
            }

            // Verify all portfolio sectors exist in benchmark
            for (const auto &sector : portfolio_sectors)
            {
                if (bench_map.find(sector.sector_name) == bench_map.end())
                {
                    throw std::invalid_argument(
                        "Portfolio sector '" + sector.sector_name + "' not found in benchmark sectors");
                }
            }

            // Compute total benchmark return (weighted sum)
            double total_bench_return = 0.0;
            for (const auto &sector : benchmark_sectors)
            {
                total_bench_return += sector.weight * sector.return_value;
            }

            // Compute total portfolio return (weighted sum)
            double total_port_return = 0.0;
            for (const auto &sector : portfolio_sectors)
            {
                total_port_return += sector.weight * sector.return_value;
            }

            // Compute per-sector effects
            SinglePeriodAttribution result;
            result.portfolio_return = total_port_return;
            result.benchmark_return = total_bench_return;
            result.total_allocation = 0.0;
            result.total_selection = 0.0;
            result.total_interaction = 0.0;

            for (const auto &port_sector : portfolio_sectors)
            {
                const auto &bench_sector = *bench_map[port_sector.sector_name];

                double w_p = port_sector.weight;
                double w_b = bench_sector.weight;
                double r_p = port_sector.return_value;
                double r_b = bench_sector.return_value;

                AttributionEffect effect;
                effect.sector_name = port_sector.sector_name;
                effect.allocation = (w_p - w_b) * (r_b - total_bench_return);
                effect.selection = w_b * (r_p - r_b);
                effect.interaction = (w_p - w_b) * (r_p - r_b);
                effect.total = effect.allocation + effect.selection + effect.interaction;

                result.total_allocation += effect.allocation;
                result.total_selection += effect.selection;
                result.total_interaction += effect.interaction;

                result.sector_effects.push_back(effect);
            }

            result.total_active = result.total_allocation + result.total_selection + result.total_interaction;

            return result;
        }

        // ===================================================================
        // Multi-Period Attribution
        // ===================================================================

        void Attribution::add_period(
            const std::vector<SectorAllocation> &portfolio_sectors,
            const std::vector<SectorAllocation> &benchmark_sectors)
        {
            // Compute single-period attribution and store it
            periods_.push_back(single_period(portfolio_sectors, benchmark_sectors));
        }

        MultiPeriodAttribution Attribution::linked_attribution() const
        {
            if (periods_.empty())
            {
                throw std::runtime_error(
                    "No attribution periods have been added; "
                    "call add_period() before linked_attribution()");
            }

            MultiPeriodAttribution result;
            result.period_results = periods_;
            result.num_periods = static_cast<int>(periods_.size());

            // Compute total compounded returns
            double port_compound = 1.0;
            double bench_compound = 1.0;
            std::vector<double> port_returns;
            std::vector<double> bench_returns;

            for (const auto &period : periods_)
            {
                port_compound *= (1.0 + period.portfolio_return);
                bench_compound *= (1.0 + period.benchmark_return);
                port_returns.push_back(period.portfolio_return);
                bench_returns.push_back(period.benchmark_return);
            }

            result.total_portfolio_return = port_compound - 1.0;
            result.total_benchmark_return = bench_compound - 1.0;

            // Compute Carino smoothing factors
            auto factors = carino_factors(port_returns, bench_returns);

            // Aggregate linked effects per sector
            result.cumulative_allocation = 0.0;
            result.cumulative_selection = 0.0;
            result.cumulative_interaction = 0.0;

            // Collect all sector names
            std::unordered_map<std::string, AttributionEffect> sector_accum;

            for (int t = 0; t < result.num_periods; ++t)
            {
                double factor = factors[t];
                const auto &period = periods_[t];

                for (const auto &effect : period.sector_effects)
                {
                    auto &accum = sector_accum[effect.sector_name];
                    accum.sector_name = effect.sector_name;
                    accum.allocation += effect.allocation * factor;
                    accum.selection += effect.selection * factor;
                    accum.interaction += effect.interaction * factor;
                    accum.total += effect.total * factor;
                }

                result.cumulative_allocation += period.total_allocation * factor;
                result.cumulative_selection += period.total_selection * factor;
                result.cumulative_interaction += period.total_interaction * factor;
            }

            result.cumulative_active = result.cumulative_allocation + result.cumulative_selection + result.cumulative_interaction;

            // Copy to the map
            for (const auto &[name, effect] : sector_accum)
            {
                result.cumulative_sector_effects[name] = effect;
            }

            return result;
        }

        int Attribution::period_count() const
        {
            return static_cast<int>(periods_.size());
        }

        void Attribution::clear()
        {
            periods_.clear();
        }

        // ===================================================================
        // Export
        // ===================================================================

        std::string Attribution::report() const
        {
            if (periods_.empty())
            {
                throw std::runtime_error(
                    "No attribution periods have been added");
            }

            std::ostringstream oss;
            oss << std::fixed;

            oss << "Performance Attribution Report\n";
            oss << "==============================\n\n";

            // Per-period summary
            for (int t = 0; t < static_cast<int>(periods_.size()); ++t)
            {
                const auto &period = periods_[t];
                oss << "Period " << (t + 1) << ":\n";
                oss << "  Portfolio Return:  " << std::setprecision(4)
                    << period.portfolio_return * 100.0 << "%\n";
                oss << "  Benchmark Return:  " << std::setprecision(4)
                    << period.benchmark_return * 100.0 << "%\n";
                oss << "  Active Return:     " << std::setprecision(4)
                    << period.total_active * 100.0 << "%\n";
                oss << "  Allocation:        " << std::setprecision(4)
                    << period.total_allocation * 100.0 << "%\n";
                oss << "  Selection:         " << std::setprecision(4)
                    << period.total_selection * 100.0 << "%\n";
                oss << "  Interaction:       " << std::setprecision(4)
                    << period.total_interaction * 100.0 << "%\n";

                oss << "  Sector Detail:\n";
                for (const auto &effect : period.sector_effects)
                {
                    oss << "    " << std::left << std::setw(16) << effect.sector_name
                        << " Alloc: " << std::setprecision(4)
                        << std::setw(8) << std::right << effect.allocation * 100.0 << "%"
                        << " Sel: " << std::setw(8) << effect.selection * 100.0 << "%"
                        << " Inter: " << std::setw(8) << effect.interaction * 100.0 << "%"
                        << " Total: " << std::setw(8) << effect.total * 100.0 << "%\n";
                }
                oss << "\n";
            }

            // Linked totals (if multi-period)
            if (periods_.size() > 1)
            {
                auto linked = linked_attribution();
                oss << "Linked Cumulative Attribution (Carino Method):\n";
                oss << "  Total Portfolio Return:  " << std::setprecision(4)
                    << linked.total_portfolio_return * 100.0 << "%\n";
                oss << "  Total Benchmark Return:  " << std::setprecision(4)
                    << linked.total_benchmark_return * 100.0 << "%\n";
                oss << "  Cumulative Active:       " << std::setprecision(4)
                    << linked.cumulative_active * 100.0 << "%\n";
                oss << "  Cumulative Allocation:   " << std::setprecision(4)
                    << linked.cumulative_allocation * 100.0 << "%\n";
                oss << "  Cumulative Selection:    " << std::setprecision(4)
                    << linked.cumulative_selection * 100.0 << "%\n";
                oss << "  Cumulative Interaction:  " << std::setprecision(4)
                    << linked.cumulative_interaction * 100.0 << "%\n";
            }

            return oss.str();
        }

        // ===================================================================
        // Private Helpers
        // ===================================================================

        std::vector<double> Attribution::carino_factors(
            const std::vector<double> &portfolio_returns,
            const std::vector<double> &benchmark_returns)
        {
            int n = static_cast<int>(portfolio_returns.size());
            std::vector<double> factors(n, 1.0);

            if (n <= 1)
            {
                return factors;
            }

            // Compute total compounded returns
            double port_compound = 1.0;
            double bench_compound = 1.0;
            for (int t = 0; t < n; ++t)
            {
                port_compound *= (1.0 + portfolio_returns[t]);
                bench_compound *= (1.0 + benchmark_returns[t]);
            }

            double total_port = port_compound - 1.0;
            double total_bench = bench_compound - 1.0;

            // Carino logarithmic linking factor:
            //
            // k_t = [ln(1 + r_p,t) - ln(1 + r_b,t)] / (r_p,t - r_b,t)
            //        * (r_P - r_B) / [ln(1 + r_P) - ln(1 + r_B)]
            //
            // where r_p,t and r_b,t are single-period returns,
            //       r_P and r_B are total compounded returns.
            //
            // Special cases: if any denominator is near zero, use L'Hopital
            // limit (the factor approaches 1).

            double total_active = total_port - total_bench;
            double total_log_ratio;

            if (std::abs(total_active) < 1e-14)
            {
                // Total active return near zero; factors are all 1.0
                return factors;
            }

            double log_port = std::log(1.0 + total_port);
            double log_bench = std::log(1.0 + total_bench);
            double log_diff = log_port - log_bench;

            if (std::abs(log_diff) < 1e-14)
            {
                total_log_ratio = 1.0;
            }
            else
            {
                total_log_ratio = total_active / log_diff;
            }

            for (int t = 0; t < n; ++t)
            {
                double r_p = portfolio_returns[t];
                double r_b = benchmark_returns[t];
                double active = r_p - r_b;

                double period_factor;
                if (std::abs(active) < 1e-14)
                {
                    // Single-period active near zero: use limit = 1 / (1+r_p)
                    // (from L'Hopital on the log ratio)
                    period_factor = 1.0 / (1.0 + r_p);
                }
                else
                {
                    double log_p = std::log(1.0 + r_p);
                    double log_b = std::log(1.0 + r_b);
                    period_factor = (log_p - log_b) / active;
                }

                factors[t] = period_factor / (1.0 / total_log_ratio);
            }

            return factors;
        }

    } // namespace analytics
} // namespace portfolio
