/**
 * @file optimizer_interface.hpp
 * @brief Abstract interface for portfolio optimization methods
 *
 * Provides a common interface for different portfolio optimization
 * techniques. All optimizers must implement the optimize method.
 *
 * Design follows the RiskModel pattern from the risk layer.
 *
 * Thread Safety: Implementations are expected to be thread-safe for
 * read-only operations.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include <memory>
#include <vector>
#include <nlohmann/json.hpp>

namespace portfolio
{
    namespace optimizer
    {

        /**
         * @struct OptimizationConstraints
         * @brief Container for portfolio constraints
         */
        struct OptimizationConstraints
        {
            double min_weight = 0.0;   ///< Minimum asset weight
            double max_weight = 1.0;   ///< Maximum asset weight
            bool long_only = true;     ///< No short positions
            bool sum_to_one = true;    ///< Weights sum to 1
            double max_turnover = 1.0; ///< Maximum turnover from current
            std::vector<struct GroupConstraint> group_constraints; ///< Group exposure constraints

            /**
             * @brief Validate constraints
             * @throws std::invalid_argument if constraints are inconsistent
             */
            void validate() const;

            /**
             * @brief Create from JSON configuration
             */
            static OptimizationConstraints from_json(const nlohmann::json &j);
        };

        /**
         * @struct GroupConstraint
         * @brief Constraint on a group of assets (sector, industry, factor, etc.)
         *
         * Enforces: min_weight <= sum(w_i : i in asset_indices) <= max_weight
         *
         * Example: GroupConstraint{"Technology", {0,1,2,3}, 0.20, 0.35}
         *          ensures 20-35% allocation to assets 0,1,2,3 (Tech sector)
         */
        struct GroupConstraint {
            std::string name;                ///< Group name (e.g., "Technology Sector")
            std::vector<int> asset_indices;  ///< Asset indices in this group
            double min_weight;               ///< Minimum total weight (e.g., 0.20 = 20%)
            double max_weight;               ///< Maximum total weight (e.g., 0.35 = 35%)

            /**
             * @brief Validate constraint
             * @throws std::invalid_argument if invalid
             */
            void validate() const;

            /**
             * @brief Convert to JSON
             */
            nlohmann::json to_json() const;

            /**
             * @brief Create from JSON
             */
            static GroupConstraint from_json(const nlohmann::json& j);
        };

        /**
         * @struct OptimizationResult
         * @brief Container for optimization results
         */
        struct OptimizationResult
        {
            Eigen::VectorXd weights; ///< Optimal portfolio weights
            double expected_return;  ///< Portfolio expected return
            double volatility;       ///< Portfolio volatility
            double sharpe_ratio;     ///< Sharpe ratio
            bool success;            ///< Optimization succeeded
            std::string message;     ///< Status message
            int iterations;          ///< Number of iterations
            double objective_value;  ///< Final objective value

            /**
             * @brief Default constructor
             */
            OptimizationResult();

            /**
             * @brief Check if result is valid
             */
            bool is_valid() const;

            /**
             * @brief Print summary statistics
             */
            void print_summary() const;
        };

        /**
         * @class OptimizerInterface
         * @brief Abstract base class for portfolio optimizers
         *
         * Defines the interface for portfolio optimization. Concrete
         * implementations include mean-variance, QUBO, and risk parity.
         *
         * Usage Example:
         * @code
         * auto optimizer = std::make_unique<MeanVarianceOptimizer>();
         * auto result = optimizer->optimize(returns, covariance, constraints);
         * result.print_summary();
         * @endcode
         */
        class OptimizerInterface
        {
        public:
            virtual ~OptimizerInterface() = default;

            /**
             * @brief Optimize portfolio weights
             * @param expected_returns Expected returns for each asset (N x 1)
             * @param covariance Covariance matrix (N x N)
             * @param constraints Portfolio constraints
             * @param current_weights Current portfolio weights (for turnover)
             * @return OptimizationResult structure
             * @throws std::invalid_argument if inputs are invalid
             * @throws std::runtime_error if optimization fails
             */
            virtual OptimizationResult optimize(
                const Eigen::VectorXd &expected_returns,
                const Eigen::MatrixXd &covariance,
                const OptimizationConstraints &constraints,
                const Eigen::VectorXd &current_weights = Eigen::VectorXd()) const = 0;

            /**
             * @brief Get optimizer name
             * @return String identifier for the optimizer type
             */
            virtual std::string get_name() const = 0;

            /**
             * @brief Get optimizer parameters as JSON
             * @return JSON object with optimizer configuration
             */
            virtual nlohmann::json get_parameters() const = 0;

            /**
             * @brief Validate input data
             * @param expected_returns Expected returns vector
             * @param covariance Covariance matrix
             * @throws std::invalid_argument if validation fails
             */
            static void validate_inputs(
                const Eigen::VectorXd &expected_returns,
                const Eigen::MatrixXd &covariance);

        protected:
            /**
             * @brief Check if weights satisfy constraints
             * @param weights Portfolio weights
             * @param constraints Constraints to check
             * @param tolerance Numerical tolerance
             * @return true if constraints satisfied
             */
            static bool check_constraints(
                const Eigen::VectorXd &weights,
                const OptimizationConstraints &constraints,
                double tolerance = 1e-6);

            /**
             * @brief Calculate portfolio statistics
             * @param weights Portfolio weights
             * @param expected_returns Expected returns
             * @param covariance Covariance matrix
             * @param risk_free_rate Risk-free rate for Sharpe
             * @return Result structure with calculated statistics
             */
            static OptimizationResult calculate_statistics(
                const Eigen::VectorXd &weights,
                const Eigen::VectorXd &expected_returns,
                const Eigen::MatrixXd &covariance,
                double risk_free_rate = 0.0);
        };

    } // namespace optimizer
} // namespace portfolio