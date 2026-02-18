/**
 * @file attribution.hpp
 * @brief Brinson-Fachler performance attribution analysis.
 *
 * Implements single-period Brinson-Fachler attribution decomposing
 * active return into allocation, selection, and interaction effects.
 * Multi-period linking uses the Carino smoothing method to ensure
 * single-period effects compound to match the total active return.
 *
 * This is a stretch goal component -- implement if time permits after
 * the core analytics classes are complete and tested.
 *
 * Attribution model:
 *   Allocation  = (w_p - w_b) * (R_b_sector - R_b_total)
 *   Selection   = w_b * (R_p_sector - R_b_sector)
 *   Interaction = (w_p - w_b) * (R_p_sector - R_b_sector)
 *   Active      = Allocation + Selection + Interaction
 */

#ifndef PORTFOLIO_ANALYTICS_ATTRIBUTION_HPP
#define PORTFOLIO_ANALYTICS_ATTRIBUTION_HPP

#include <map>
#include <string>
#include <vector>

namespace portfolio
{
    namespace analytics
    {

        /**
         * @struct SectorAllocation
         * @brief Portfolio or benchmark allocation to a single sector.
         */
        struct SectorAllocation
        {
            std::string sector_name; ///< Sector identifier
            double weight;           ///< Weight in the portfolio or benchmark (0 to 1)
            double return_value;     ///< Return of that sector for the period
        };

        /**
         * @struct AttributionEffect
         * @brief Single-period attribution effects for one sector.
         */
        struct AttributionEffect
        {
            std::string sector_name; ///< Sector identifier
            double allocation;       ///< Allocation effect
            double selection;        ///< Selection effect
            double interaction;      ///< Interaction effect
            double total;            ///< Total active effect (sum of above)
        };

        /**
         * @struct SinglePeriodAttribution
         * @brief Complete single-period attribution result.
         */
        struct SinglePeriodAttribution
        {
            std::vector<AttributionEffect> sector_effects; ///< Per-sector effects

            double total_allocation;  ///< Sum of allocation effects across sectors
            double total_selection;   ///< Sum of selection effects across sectors
            double total_interaction; ///< Sum of interaction effects across sectors
            double total_active;      ///< Total active return (should match R_p - R_b)

            double portfolio_return; ///< Total portfolio return for the period
            double benchmark_return; ///< Total benchmark return for the period
        };

        /**
         * @struct MultiPeriodAttribution
         * @brief Multi-period linked attribution using Carino method.
         */
        struct MultiPeriodAttribution
        {
            std::vector<SinglePeriodAttribution> period_results; ///< Per-period results

            /// Linked cumulative effects per sector over all periods
            std::map<std::string, AttributionEffect> cumulative_sector_effects;

            double cumulative_allocation;  ///< Linked total allocation effect
            double cumulative_selection;   ///< Linked total selection effect
            double cumulative_interaction; ///< Linked total interaction effect
            double cumulative_active;      ///< Linked total active return

            double total_portfolio_return; ///< Compounded portfolio return
            double total_benchmark_return; ///< Compounded benchmark return
            int num_periods;               ///< Number of attribution periods
        };

        /**
         * @class Attribution
         * @brief Brinson-Fachler performance attribution with Carino linking.
         *
         * Decomposes the difference between portfolio and benchmark returns
         * into allocation, selection, and interaction effects at the sector
         * level. Supports both single-period and multi-period attribution
         * with arithmetic linking (Carino smoothing).
         *
         * Usage:
         * @code
         *   // Single period
         *   auto result = Attribution::single_period(
         *       portfolio_sectors, benchmark_sectors);
         *
         *   // Multi-period
         *   Attribution attr;
         *   attr.add_period(portfolio_sectors_t1, benchmark_sectors_t1);
         *   attr.add_period(portfolio_sectors_t2, benchmark_sectors_t2);
         *   auto linked = attr.linked_attribution();
         * @endcode
         *
         * Thread safety: single_period is a static pure function. Instances
         * are not thread-safe during add_period calls but are safe to read
         * after all periods are added.
         */
        class Attribution
        {
        public:
            // ---------------------------------------------------------------
            // Constructors
            // ---------------------------------------------------------------

            /** @brief Default constructor for multi-period accumulation. */
            Attribution() = default;

            /** @brief Default destructor. */
            ~Attribution() = default;

            // ---------------------------------------------------------------
            // Single-Period Attribution
            // ---------------------------------------------------------------

            /**
             * @brief Compute single-period Brinson-Fachler attribution.
             * @param portfolio_sectors Sector allocations and returns for the portfolio.
             * @param benchmark_sectors Sector allocations and returns for the benchmark.
             * @return SinglePeriodAttribution with per-sector and aggregate effects.
             * @throws std::invalid_argument If sector lists are empty or sector names
             *         do not match between portfolio and benchmark.
             *
             * @note Weights in both portfolio and benchmark should sum to 1.0.
             *       A warning is logged if they deviate by more than 1e-6.
             */
            static SinglePeriodAttribution single_period(
                const std::vector<SectorAllocation> &portfolio_sectors,
                const std::vector<SectorAllocation> &benchmark_sectors);

            // ---------------------------------------------------------------
            // Multi-Period Attribution
            // ---------------------------------------------------------------

            /**
             * @brief Add a period for multi-period attribution.
             * @param portfolio_sectors Sector allocations and returns for this period.
             * @param benchmark_sectors Sector allocations and returns for this period.
             * @throws std::invalid_argument If sector lists are empty or inconsistent.
             */
            void add_period(const std::vector<SectorAllocation> &portfolio_sectors,
                            const std::vector<SectorAllocation> &benchmark_sectors);

            /**
             * @brief Compute multi-period linked attribution using Carino method.
             * @return MultiPeriodAttribution with linked cumulative effects.
             * @throws std::runtime_error If no periods have been added.
             *
             * The Carino method applies smoothing factors so that single-period
             * effects compound geometrically to match the total active return.
             */
            MultiPeriodAttribution linked_attribution() const;

            /**
             * @brief Get the number of periods added.
             * @return Period count.
             */
            int period_count() const;

            /**
             * @brief Reset all accumulated periods.
             */
            void clear();

            // ---------------------------------------------------------------
            // Export
            // ---------------------------------------------------------------

            /**
             * @brief Generate a formatted attribution report.
             * @return Multi-line string with per-period and cumulative effects.
             * @throws std::runtime_error If no periods have been added.
             */
            std::string report() const;

        private:
            // ---------------------------------------------------------------
            // Private helpers
            // ---------------------------------------------------------------

            /**
             * @brief Compute Carino smoothing factors for multi-period linking.
             * @param portfolio_returns Per-period total portfolio returns.
             * @param benchmark_returns Per-period total benchmark returns.
             * @return Vector of smoothing factors, one per period.
             */
            static std::vector<double> carino_factors(
                const std::vector<double> &portfolio_returns,
                const std::vector<double> &benchmark_returns);

            // ---------------------------------------------------------------
            // Member variables
            // ---------------------------------------------------------------

            std::vector<SinglePeriodAttribution> periods_;
        };

    } // namespace analytics
} // namespace portfolio

#endif // PORTFOLIO_ANALYTICS_ATTRIBUTION_HPP