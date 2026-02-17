/**
 * @file drawdown_analysis.hpp
 * @brief Detailed drawdown analysis for portfolio NAV series.
 *
 * Provides identification and characterization of all drawdown events
 * in a NAV series, including the top-N deepest drawdowns, duration
 * analysis, and overlap detection. Complements the basic drawdown
 * metrics in PerformanceMetrics with deeper structural analysis.
 *
 * A drawdown event is defined as a contiguous period where NAV is
 * below a prior peak. Each event has a peak, trough, and optional
 * recovery point.
 */

#ifndef PORTFOLIO_ANALYTICS_DRAWDOWN_ANALYSIS_HPP
#define PORTFOLIO_ANALYTICS_DRAWDOWN_ANALYSIS_HPP

#include <string>
#include <vector>

namespace portfolio
{
    namespace analytics
    {

        /**
         * @struct DrawdownEvent
         * @brief A single drawdown event with full characterization.
         *
         * Represents one peak-to-trough-to-recovery cycle. If recovery has not
         * occurred by the end of the series, recovery_index is -1 and
         * recovery_date is empty.
         */
        struct DrawdownEvent
        {
            int peak_index;     ///< Index of the peak NAV before decline
            int trough_index;   ///< Index of the lowest NAV during the event
            int recovery_index; ///< Index where NAV recovers to peak (-1 if unrecovered)

            std::string peak_date;     ///< Date at peak
            std::string trough_date;   ///< Date at trough
            std::string recovery_date; ///< Date at recovery (empty if unrecovered)

            double peak_nav;   ///< NAV at peak
            double trough_nav; ///< NAV at trough
            double depth;      ///< Drawdown depth as positive fraction (e.g., 0.15 = 15%)

            int decline_days;  ///< Trading days from peak to trough
            int recovery_days; ///< Trading days from trough to recovery (-1 if unrecovered)
            int total_days;    ///< Trading days from peak to recovery (-1 if unrecovered)
        };

        /**
         * @struct DrawdownSummary
         * @brief Aggregate statistics across all drawdown events.
         */
        struct DrawdownSummary
        {
            int total_events;             ///< Number of distinct drawdown events
            double average_depth;         ///< Mean drawdown depth across all events
            double average_decline_days;  ///< Mean decline duration in trading days
            double average_recovery_days; ///< Mean recovery duration (excluding unrecovered)
            int unrecovered_count;        ///< Number of drawdowns not recovered by end
            double time_in_drawdown_pct;  ///< Fraction of total period spent in drawdown
            double max_depth;             ///< Deepest drawdown across all events
            double longest_decline_days;  ///< Longest peak-to-trough duration
            double longest_recovery_days; ///< Longest trough-to-recovery duration
        };

        /**
         * @class DrawdownAnalysis
         * @brief Identifies and characterizes all drawdown events in a NAV series.
         *
         * Decomposes the NAV series into discrete drawdown events, computes
         * per-event statistics, and provides aggregate summaries including
         * the top-N deepest drawdowns.
         *
         * Usage:
         * @code
         *   DrawdownAnalysis analysis(nav_series, dates);
         *   auto top5 = analysis.top_drawdowns(5);
         *   auto summary = analysis.summary();
         *   double pct_underwater = summary.time_in_drawdown_pct;
         * @endcode
         *
         * Thread safety: Instances are effectively immutable after construction.
         *
         * Performance: Analysis of 2500-day series completes in under 5ms.
         */
        class DrawdownAnalysis
        {
        public:
            // ---------------------------------------------------------------
            // Constructors
            // ---------------------------------------------------------------

            /**
             * @brief Construct from NAV series and dates.
             * @param nav_series Daily NAV values (at least 2 elements).
             * @param dates Date strings aligned with NAV series (YYYY-MM-DD format).
             * @throws std::invalid_argument If nav_series has fewer than 2 elements,
             *         or sizes of nav_series and dates do not match.
             */
            DrawdownAnalysis(const std::vector<double> &nav_series,
                             const std::vector<std::string> &dates);

            /** @brief Default destructor. */
            ~DrawdownAnalysis() = default;

            // ---------------------------------------------------------------
            // Event Access
            // ---------------------------------------------------------------

            /**
             * @brief Get all drawdown events in chronological order.
             * @return Vector of DrawdownEvent sorted by peak_index ascending.
             */
            const std::vector<DrawdownEvent> &all_events() const;

            /**
             * @brief Get the top-N deepest drawdown events.
             * @param n Number of events to return (clamped to total events).
             * @return Vector of DrawdownEvent sorted by depth descending.
             * @throws std::invalid_argument If n < 1.
             */
            std::vector<DrawdownEvent> top_drawdowns(int n) const;

            /**
             * @brief Get the single worst drawdown event.
             * @return DrawdownEvent with the deepest drawdown.
             * @throws std::runtime_error If no drawdown events exist.
             */
            const DrawdownEvent &worst_drawdown() const;

            /**
             * @brief Get the number of drawdown events.
             * @return Total count of distinct drawdown events.
             */
            int event_count() const;

            // ---------------------------------------------------------------
            // Aggregate Statistics
            // ---------------------------------------------------------------

            /**
             * @brief Compute aggregate drawdown summary statistics.
             * @return DrawdownSummary with averages and extremes across all events.
             */
            DrawdownSummary summary() const;

            /**
             * @brief Get the full underwater (drawdown) curve.
             * @return Vector of drawdown values at each point. Values are non-positive;
             *         0.0 means at or above previous peak, -0.10 means 10% below peak.
             */
            const std::vector<double> &underwater_curve() const;

            // ---------------------------------------------------------------
            // Export
            // ---------------------------------------------------------------

            /**
             * @brief Generate a formatted string report of all drawdown events.
             * @param max_events Maximum number of events to include (default all).
             * @return Formatted multi-line string with drawdown details.
             */
            std::string report(int max_events = -1) const;

        private:
            // ---------------------------------------------------------------
            // Private helpers
            // ---------------------------------------------------------------

            /**
             * @brief Identify all drawdown events from the NAV series.
             */
            void identify_events();

            /**
             * @brief Compute the underwater curve from the NAV series.
             */
            void compute_underwater_curve();

            // ---------------------------------------------------------------
            // Member variables
            // ---------------------------------------------------------------

            std::vector<double> nav_series_;
            std::vector<std::string> dates_;
            std::vector<DrawdownEvent> events_;
            std::vector<double> underwater_curve_;
            int worst_event_index_; ///< Index into events_ of the deepest drawdown
        };

    } // namespace analytics
} // namespace portfolio

#endif // PORTFOLIO_ANALYTICS_DRAWDOWN_ANALYSIS_HPP