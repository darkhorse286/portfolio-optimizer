/**
 * @file mean_variance_optimizer.hpp
 * @brief Mean-variance portfolio optimizer (Markowitz)
 *
 * Implements classical Markowitz mean-variance optimization using
 * custom quadratic programming solver. Supports multiple objective functions:
 * - Minimum variance
 * - Maximum Sharpe ratio
 * - Target return with minimum variance
 * - Risk aversion utility maximization
 *
 * Mathematical Formulation:
 *
 * Minimize:     (1/2) * w^T * Sigma * w - lambda * mu^T * w
 * Subject to:   sum(w_i) = 1
 *              w_min <= w_i <= w_max
 *              |w_i - w_i^old| <= turnover
 *
 * where:
 * - w: portfolio weights
 * - Sigma: covariance matrix
 * - mu: expected returns
 * - lambda: risk aversion parameter
 *
 * Performance: approximately 50ms for 50 assets with constraints
 */

#pragma once

#include "optimizer/optimizer_interface.hpp"
#include "optimizer/quadratic_solver.hpp"

namespace portfolio
{
    namespace optimizer
    {

        /**
         * @enum ObjectiveType
         * @brief Type of optimization objective
         */
        enum class ObjectiveType
        {
            MIN_VARIANCE,      ///< Minimize variance only
            MAX_SHARPE,        ///< Maximize Sharpe ratio
            TARGET_RETURN,     ///< Target return with min variance
            RISK_AVERSION      ///< Utility: return - lambda*variance
        };

        /**
         * @class MeanVarianceOptimizer
         * @brief Markowitz mean-variance portfolio optimizer
         *
         * Implements classical mean-variance optimization with support
         * for various constraints and objective functions.
         *
         * Solver: Custom quadratic programming implementation
         *
         * Key Features:
         * - Multiple objective functions
         * - Box constraints on weights
         * - Turnover constraints
         * - Long-only or long-short
         *
         * Usage Example:
         * @code
         * MeanVarianceOptimizer optimizer(ObjectiveType::MAX_SHARPE, 0.02);
         * optimizer.set_risk_aversion(1.0);
         * 
         * OptimizationConstraints constraints;
         * constraints.min_weight = 0.0;
         * constraints.max_weight = 0.3;
         * 
         * auto result = optimizer.optimize(returns, covariance, constraints);
         * std::cout << "Sharpe: " << result.sharpe_ratio << "\n";
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class MeanVarianceOptimizer : public OptimizerInterface
        {
        public:
            /**
             * @brief Construct mean-variance optimizer
             * @param objective Optimization objective type
             * @param risk_free_rate Risk-free rate for Sharpe ratio
             * @throws std::invalid_argument if parameters invalid
             */
            explicit MeanVarianceOptimizer(
                ObjectiveType objective = ObjectiveType::MAX_SHARPE,
                double risk_free_rate = 0.0);

            /**
             * @brief Destructor
             */
            ~MeanVarianceOptimizer() override = default;

            /**
             * @brief Optimize portfolio weights
             * @param expected_returns Expected returns (N x 1)
             * @param covariance Covariance matrix (N x N)
             * @param constraints Portfolio constraints
             * @param current_weights Current weights (for turnover)
             * @return Optimization result
             * @throws std::invalid_argument if inputs invalid
             * @throws std::runtime_error if optimization fails
             */
            OptimizationResult optimize(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                const Eigen::VectorXd& current_weights = Eigen::VectorXd()) const override;

            /**
             * @brief Get optimizer name
             * @return "MeanVarianceOptimizer"
             */
            std::string get_name() const override;
            
            /**
             * @brief Get parameters as JSON
             */
            nlohmann::json get_parameters() const override;

            /**
             * @brief Set risk aversion parameter
             * @param lambda Risk aversion (lambda > 0)
             * @throws std::invalid_argument if lambda <= 0
             *
             * Used for RISK_AVERSION objective.
             * Higher values = more risk averse.
             */
            void set_risk_aversion(double lambda);

            /**
             * @brief Set target return
             * @param target_return Target portfolio return
             *
             * Used for TARGET_RETURN objective.
             */
            void set_target_return(double target_return);
            
            /**
             * @brief Get risk aversion parameter
             */
            double get_risk_aversion() const { return risk_aversion_; }
            
            /**
             * @brief Get target return
             */
            double get_target_return() const { return target_return_; }
            
            /**
             * @brief Get objective type
             */
            ObjectiveType get_objective() const { return objective_; }

        private:
            ObjectiveType objective_;       ///< Optimization objective
            double risk_free_rate_;         ///< Risk-free rate
            double risk_aversion_;          ///< Risk aversion parameter
            double target_return_;          ///< Target return

            /**
             * @brief Optimize for minimum variance
             */
            OptimizationResult optimize_min_variance(
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                const Eigen::VectorXd& current_weights) const;

            /**
             * @brief Optimize for maximum Sharpe ratio
             */
            OptimizationResult optimize_max_sharpe(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                const Eigen::VectorXd& current_weights) const;

            /**
             * @brief Optimize for target return
             */
            OptimizationResult optimize_target_return(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                const Eigen::VectorXd& current_weights) const;

            /**
             * @brief Optimize with risk aversion utility
             */
            OptimizationResult optimize_risk_aversion(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                const Eigen::VectorXd& current_weights) const;

            /**
             * @brief Validate optimizer parameters
             */
            void validate_parameters() const;
        };

    } // namespace optimizer
} // namespace portfolio