/**
 * @file ledoit_wolf_shrinkage.hpp
 * @brief Ledoit-Wolf shrinkage covariance estimator
 *
 * Implements the Ledoit-Wolf (2004) covariance shrinkage method, which
 * combines the sample covariance with a structured target (identity matrix
 * or constant correlation). This reduces estimation error when the number
 * of observations is small relative to the number of assets.
 *
 * Formula:
 *     Cov_shrunk = δ * F + (1-δ) * S
 * where S is sample covariance, F is target, δ is shrinkage intensity
 *
 * Reference: Ledoit & Wolf (2004), "Honey, I Shrunk the Sample Covariance Matrix"
 * Performance: ~5ms for 100x100 matrix on modern CPU
 */

#pragma once

#include "risk_model.hpp"

namespace portfolio
{
    namespace risk
    {

        /**
         * @enum ShrinkageTarget
         * @brief Target matrix for shrinkage
         */
        enum class ShrinkageTarget
        {
            IDENTITY,             ///< Identity matrix (uncorrelated assets)
            CONSTANT_VARIANCE,    ///< Diagonal with constant variance
            CONSTANT_CORRELATION, ///< Single-factor model
            MARKET_MODEL          ///< Sharpe's market model
        };

        /**
         * @class LedoitWolfShrinkage
         * @brief Ledoit-Wolf shrinkage estimator
         *
         * Shrinks the sample covariance matrix towards a structured target to
         * improve estimation accuracy, especially when T < N or T/N is small.
         * The optimal shrinkage intensity is estimated analytically.
         *
         * Properties:
         * - Reduces estimation error for small samples
         * - Always positive definite (if target is positive definite)
         * - Optimal shrinkage intensity computed automatically
         * - Better out-of-sample performance than sample covariance
         *
         * When to use:
         * - T/N < 10 (observations/assets ratio is small)
         * - Need stable, well-conditioned covariance matrix
         * - Out-of-sample performance is critical
         *
         * Usage Example:
         * @code
         * LedoitWolfShrinkage lw(ShrinkageTarget::CONSTANT_CORRELATION);
         * Eigen::MatrixXd cov = lw.estimate_covariance(returns);
         * double shrinkage = lw.get_shrinkage_intensity();
         * std::cout << "Optimal shrinkage: " << shrinkage << std::endl;
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class LedoitWolfShrinkage : public RiskModel
        {
        public:
            /**
             * @brief Construct Ledoit-Wolf estimator
             * @param target Shrinkage target type
             * @param shrinkage_override Optional fixed shrinkage (0 ≤ δ ≤ 1)
             *
             * If shrinkage_override < 0 (default), optimal shrinkage is computed.
             * Otherwise, the specified shrinkage is used (useful for testing).
             */
            explicit LedoitWolfShrinkage(
                ShrinkageTarget target = ShrinkageTarget::CONSTANT_CORRELATION,
                double shrinkage_override = -1.0);

            /**
             * @brief Destructor
             */
            ~LedoitWolfShrinkage() override = default;

            /**
             * @brief Estimate shrinkage covariance matrix
             * @param returns Matrix of returns (T x N: observations x assets)
             * @return Shrunk covariance matrix (N x N)
             * @throws std::invalid_argument if returns is empty or has < 2 observations
             *
             * Computes optimal shrinkage intensity and shrinks sample covariance.
             * Time complexity: O(N^2 * T)
             * Space complexity: O(N^2)
             */
            Eigen::MatrixXd estimate_covariance(
                const Eigen::MatrixXd &returns) const override;

            /**
             * @brief Get model name
             * @return "LedoitWolfShrinkage"
             */
            std::string get_name() const override;

            /**
             * @brief Get shrinkage intensity from last estimation
             * @return Shrinkage coefficient δ ∈ [0, 1]
             *
             * @note Only valid after calling estimate_covariance()
             * @note Returns 0.0 if not yet estimated
             */
            double get_shrinkage_intensity() const { return last_shrinkage_; }

            /**
             * @brief Get shrinkage target type
             * @return Current shrinkage target
             */
            ShrinkageTarget get_target() const { return target_; }

        private:
            ShrinkageTarget target_;        ///< Shrinkage target type
            double shrinkage_override_;     ///< Fixed shrinkage (-1 for auto)
            mutable double last_shrinkage_; ///< Last computed shrinkage intensity

            /**
             * @brief Compute shrinkage target matrix
             * @param returns Return matrix
             * @param sample_cov Sample covariance matrix
             * @return Target matrix F
             */
            Eigen::MatrixXd compute_target(
                const Eigen::MatrixXd &returns,
                const Eigen::MatrixXd &sample_cov) const;

            /**
             * @brief Compute optimal shrinkage intensity
             * @param returns Return matrix
             * @param sample_cov Sample covariance
             * @param target Target matrix
             * @return Optimal shrinkage δ ∈ [0, 1]
             */
            double compute_shrinkage_intensity(
                const Eigen::MatrixXd &returns,
                const Eigen::MatrixXd &sample_cov,
                const Eigen::MatrixXd &target) const;
        };

    } // namespace risk
} // namespace portfolio