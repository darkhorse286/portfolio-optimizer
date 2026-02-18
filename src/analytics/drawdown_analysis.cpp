/**
 * @file drawdown_analysis.cpp
 * @brief Implementation of the DrawdownAnalysis class.
 *
 * Decomposes a NAV series into discrete drawdown events, computing
 * per-event statistics (depth, decline duration, recovery time) and
 * aggregate summary metrics. Events are identified by tracking the
 * running peak and detecting transitions into and out of drawdown.
 */

#include "analytics/drawdown_analysis.hpp"

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

        DrawdownAnalysis::DrawdownAnalysis(const std::vector<double> &nav_series,
                                           const std::vector<std::string> &dates)
            : nav_series_(nav_series), dates_(dates), worst_event_index_(-1)
        {
            if (nav_series_.size() < 2)
            {
                throw std::invalid_argument(
                    "NAV series must have at least 2 elements, got: " + std::to_string(nav_series_.size()));
            }
            if (nav_series_.size() != dates_.size())
            {
                throw std::invalid_argument(
                    "NAV series size (" + std::to_string(nav_series_.size()) + ") must match dates size (" + std::to_string(dates_.size()) + ")");
            }

            compute_underwater_curve();
            identify_events();
        }

        // ===================================================================
        // Event Access
        // ===================================================================

        const std::vector<DrawdownEvent> &DrawdownAnalysis::all_events() const
        {
            return events_;
        }

        std::vector<DrawdownEvent> DrawdownAnalysis::top_drawdowns(int n) const
        {
            if (n < 1)
            {
                throw std::invalid_argument(
                    "Expected positive value for parameter 'n', got: " + std::to_string(n));
            }

            std::vector<DrawdownEvent> sorted_events(events_);
            std::sort(sorted_events.begin(), sorted_events.end(),
                      [](const DrawdownEvent &a, const DrawdownEvent &b)
                      {
                          return a.depth > b.depth; // Deepest first
                      });

            int count = std::min(n, static_cast<int>(sorted_events.size()));
            sorted_events.resize(count);
            return sorted_events;
        }

        const DrawdownEvent &DrawdownAnalysis::worst_drawdown() const
        {
            if (events_.empty())
            {
                throw std::runtime_error(
                    "No drawdown events found in the NAV series");
            }
            return events_[worst_event_index_];
        }

        int DrawdownAnalysis::event_count() const
        {
            return static_cast<int>(events_.size());
        }

        // ===================================================================
        // Aggregate Statistics
        // ===================================================================

        DrawdownSummary DrawdownAnalysis::summary() const
        {
            DrawdownSummary result{};
            result.total_events = static_cast<int>(events_.size());

            if (events_.empty())
            {
                result.time_in_drawdown_pct = 0.0;
                return result;
            }

            double sum_depth = 0.0;
            double sum_decline = 0.0;
            double sum_recovery = 0.0;
            int recovery_count = 0;
            int total_underwater_days = 0;
            double max_depth = 0.0;
            int longest_decline = 0;
            int longest_recovery = 0;

            for (const auto &event : events_)
            {
                sum_depth += event.depth;
                sum_decline += static_cast<double>(event.decline_days);

                if (event.depth > max_depth)
                {
                    max_depth = event.depth;
                }
                if (event.decline_days > longest_decline)
                {
                    longest_decline = event.decline_days;
                }

                if (event.recovery_days >= 0)
                {
                    sum_recovery += static_cast<double>(event.recovery_days);
                    ++recovery_count;
                    total_underwater_days += event.total_days;
                    if (event.recovery_days > longest_recovery)
                    {
                        longest_recovery = event.recovery_days;
                    }
                }
                else
                {
                    // Unrecovered: count days from peak to end of series
                    int days_to_end = static_cast<int>(nav_series_.size()) - 1 - event.peak_index;
                    total_underwater_days += days_to_end;
                }
            }

            int n_events = static_cast<int>(events_.size());
            result.average_depth = sum_depth / static_cast<double>(n_events);
            result.average_decline_days = sum_decline / static_cast<double>(n_events);
            result.average_recovery_days = (recovery_count > 0)
                                               ? (sum_recovery / static_cast<double>(recovery_count))
                                               : -1.0;
            result.unrecovered_count = n_events - recovery_count;
            result.time_in_drawdown_pct = static_cast<double>(total_underwater_days) / static_cast<double>(nav_series_.size() - 1);
            result.max_depth = max_depth;
            result.longest_decline_days = static_cast<double>(longest_decline);
            result.longest_recovery_days = static_cast<double>(longest_recovery);

            return result;
        }

        const std::vector<double> &DrawdownAnalysis::underwater_curve() const
        {
            return underwater_curve_;
        }

        // ===================================================================
        // Export
        // ===================================================================

        std::string DrawdownAnalysis::report(int max_events) const
        {
            std::ostringstream oss;
            oss << std::fixed;

            auto sum = summary();

            oss << "Drawdown Analysis Report\n";
            oss << "========================\n\n";

            oss << "Summary:\n";
            oss << "  Total Events:           " << sum.total_events << "\n";
            oss << "  Max Depth:              " << std::setprecision(4)
                << sum.max_depth * 100.0 << "%\n";
            oss << "  Average Depth:          " << std::setprecision(4)
                << sum.average_depth * 100.0 << "%\n";
            oss << "  Avg Decline Duration:   " << std::setprecision(1)
                << sum.average_decline_days << " days\n";
            oss << "  Avg Recovery Duration:  ";
            if (sum.average_recovery_days < 0)
            {
                oss << "N/A\n";
            }
            else
            {
                oss << std::setprecision(1) << sum.average_recovery_days << " days\n";
            }
            oss << "  Unrecovered Events:     " << sum.unrecovered_count << "\n";
            oss << "  Time in Drawdown:       " << std::setprecision(2)
                << sum.time_in_drawdown_pct * 100.0 << "%\n";
            oss << "\n";

            int events_to_show = static_cast<int>(events_.size());
            if (max_events >= 0 && max_events < events_to_show)
            {
                events_to_show = max_events;
            }

            // Show events sorted by depth
            auto sorted = top_drawdowns(events_to_show);

            oss << "Top " << events_to_show << " Drawdowns:\n";
            oss << "  " << std::left
                << std::setw(6) << "Rank"
                << std::setw(10) << "Depth"
                << std::setw(14) << "Peak Date"
                << std::setw(14) << "Trough Date"
                << std::setw(14) << "Recovery"
                << std::setw(10) << "Decline"
                << std::setw(10) << "Recovery"
                << "\n";
            oss << "  " << std::setw(6) << "" << std::setw(10) << ""
                << std::setw(14) << "" << std::setw(14) << "" << std::setw(14) << "Date"
                << std::setw(10) << "Days" << std::setw(10) << "Days" << "\n";
            oss << "  " << std::string(78, '-') << "\n";

            for (int i = 0; i < static_cast<int>(sorted.size()); ++i)
            {
                const auto &e = sorted[i];
                oss << "  " << std::left
                    << std::setw(6) << (i + 1)
                    << std::setw(10) << (std::to_string(static_cast<int>(e.depth * 10000.0 + 0.5) / 100.0) + "%")
                    << std::setw(14) << e.peak_date
                    << std::setw(14) << e.trough_date
                    << std::setw(14) << (e.recovery_date.empty() ? "Unrecovered" : e.recovery_date)
                    << std::setw(10) << e.decline_days;
                if (e.recovery_days >= 0)
                {
                    oss << std::setw(10) << e.recovery_days;
                }
                else
                {
                    oss << std::setw(10) << "N/A";
                }
                oss << "\n";
            }

            return oss.str();
        }

        // ===================================================================
        // Private Helpers
        // ===================================================================

        void DrawdownAnalysis::compute_underwater_curve()
        {
            int n = static_cast<int>(nav_series_.size());
            underwater_curve_.resize(n);

            double peak = nav_series_[0];
            for (int i = 0; i < n; ++i)
            {
                if (nav_series_[i] > peak)
                {
                    peak = nav_series_[i];
                }
                underwater_curve_[i] = (nav_series_[i] - peak) / peak;
            }
        }

        void DrawdownAnalysis::identify_events()
        {
            int n = static_cast<int>(nav_series_.size());

            // State machine: track whether we are in a drawdown event.
            // An event starts when NAV drops below the running peak and ends
            // when NAV recovers to or exceeds the peak.

            double peak = nav_series_[0];
            int peak_idx = 0;
            bool in_drawdown = false;

            // Current event being built
            int event_peak_idx = 0;
            double event_peak_nav = nav_series_[0];
            int event_trough_idx = 0;
            double event_trough_nav = nav_series_[0];

            double worst_depth = 0.0;

            for (int i = 1; i < n; ++i)
            {
                if (nav_series_[i] >= peak)
                {
                    // At or above peak
                    if (in_drawdown)
                    {
                        // Event ends: NAV recovered to peak
                        DrawdownEvent event;
                        event.peak_index = event_peak_idx;
                        event.trough_index = event_trough_idx;
                        event.recovery_index = i;
                        event.peak_date = dates_[event_peak_idx];
                        event.trough_date = dates_[event_trough_idx];
                        event.recovery_date = dates_[i];
                        event.peak_nav = event_peak_nav;
                        event.trough_nav = event_trough_nav;
                        event.depth = (event_peak_nav - event_trough_nav) / event_peak_nav;
                        event.decline_days = event_trough_idx - event_peak_idx;
                        event.recovery_days = i - event_trough_idx;
                        event.total_days = i - event_peak_idx;

                        if (event.depth > worst_depth)
                        {
                            worst_depth = event.depth;
                            worst_event_index_ = static_cast<int>(events_.size());
                        }

                        events_.push_back(event);
                        in_drawdown = false;
                    }
                    peak = nav_series_[i];
                    peak_idx = i;
                }
                else
                {
                    // Below peak
                    if (!in_drawdown)
                    {
                        // New drawdown event starts
                        in_drawdown = true;
                        event_peak_idx = peak_idx;
                        event_peak_nav = peak;
                        event_trough_idx = i;
                        event_trough_nav = nav_series_[i];
                    }
                    else
                    {
                        // Continue existing drawdown; update trough if deeper
                        if (nav_series_[i] < event_trough_nav)
                        {
                            event_trough_idx = i;
                            event_trough_nav = nav_series_[i];
                        }
                    }
                }
            }

            // Handle unrecovered drawdown at end of series
            if (in_drawdown)
            {
                DrawdownEvent event;
                event.peak_index = event_peak_idx;
                event.trough_index = event_trough_idx;
                event.recovery_index = -1;
                event.peak_date = dates_[event_peak_idx];
                event.trough_date = dates_[event_trough_idx];
                event.recovery_date = "";
                event.peak_nav = event_peak_nav;
                event.trough_nav = event_trough_nav;
                event.depth = (event_peak_nav - event_trough_nav) / event_peak_nav;
                event.decline_days = event_trough_idx - event_peak_idx;
                event.recovery_days = -1;
                event.total_days = -1;

                if (event.depth > worst_depth)
                {
                    worst_depth = event.depth;
                    worst_event_index_ = static_cast<int>(events_.size());
                }

                events_.push_back(event);
            }
        }

    } // namespace analytics
} // namespace portfolio