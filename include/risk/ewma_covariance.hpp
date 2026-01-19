/**
 * @file ewma_covariance.hpp
 * @brief Exponentially Weighted Moving Average covariance estimator
 *
 * Implements EWMA covariance estimation where recent observations are
 * weighted more heavily than older ones. Commonly used in risk management
 * to adapt quickly to changing market conditions.
 *
 * Formula:
 *     Cov_t = λ * Cov_{t-1} + (1-λ) * r_t * r_t^T
 *
 * RiskMetrics (J.P. Morgan) recommends λ = 0.94 for daily data.
 *
 * Performance: ~3ms for 100x100 matrix on modern CPU
 */

#pragma once

#include "risk_model.hpp"

namespace portfolio
{
    namespace risk
    {

        /**
         * @class EWMACovariance
         * @brief Exponentially weighted covariance estimator
         *
         * Computes covariance using exponentially decaying weights, giving more
         * importance to recent observations. Particularly useful for modeling
         * time-varying volatility and correlation.
         *
         * Properties:
         * - Adapts to regime changes faster than sample covariance
         * - Smooths out noise while tracking trends
         * - Always positive semi-definite
         * - Less sensitive to the choice of estimation window
         *
         * Common decay parameters:
         * - Daily data: λ = 0.94 (RiskMetrics standard)
         * - Weekly data: λ = 0.97
         * - Monthly data: λ = 0.98
         *
         * Effective window: ~ 1/(1-λ) observations
         * - λ = 0.94 → ~17 days effective window
         * - λ = 0.97 → ~33 days effective window
         *
         * Usage Example:
         * @code
         * EWMACovariance ewma(0.94);  // RiskMetrics standard
         * Eigen::MatrixXd cov = ewma.estimate_covariance(returns);
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class EWMACovariance : public RiskModel
        {
        public:
            /**
             * @brief Construct EWMA covariance estimator
             * @param lambda Decay factor (0 < λ < 1), typically 0.94-0.99
             * @throws std::invalid_argument if lambda not in (0, 1)
             *
             * Recommended values:
             * - 0.94: RiskMetrics standard for daily data
             * - 0.97: For slower-moving estimates
             * - 0.99: For very stable markets
             */
            explicit EWMACovariance(double lambda = 0.94);

            /**
             * @brief Destructor
             */
            ~EWMACovariance() override = default;

            /**
             * @brief Estimate EWMA covariance matrix
             * @param returns Matrix of returns (T x N: observations x assets)
             * @return EWMA covariance matrix (N x N)
             * @throws std::invalid_argument if returns is empty
             *
             * Implementation uses recursive formula for efficiency.
             * Time complexity: O(N^2 * T)
             * Space complexity: O(N^2)
             */
            Eigen::MatrixXd estimate_covariance(
                const Eigen::MatrixXd &returns) const override;

            /**
             * @brief Get model name
             * @return "EWMACovariance"
             */
            std::string get_name() const override;

            /**
             * @brief Get decay parameter
             * @return Current lambda value
             */
            double get_lambda() const { return lambda_; }

            /**
             * @brief Get effective window size
             * @return Approximate number of observations with significant weight
             *
             * Formula: 1 / (1 - λ)
             */
            double get_effective_window() const { return 1.0 / (1.0 - lambda_); }

        private:
            double lambda_; ///< Decay factor (0 < λ < 1)

            /**
             * @brief Validate lambda parameter
             * @param lambda Value to validate
             * @throws std::invalid_argument if lambda not in (0, 1)
             */
            static void validate_lambda(double lambda);
        };

    } // namespace risk
} // namespace portfolio