/**
 * @file ledoit_wolf_shrinkage.hpp
 * @brief Ledoit-Wolf shrinkage covariance estimator
 *
 * Implements the Ledoit-Wolf (2004) covariance shrinkage method, which
 * combines the sample covariance matrix with a structured target matrix
 * to reduce estimation error. This is particularly valuable when the
 * number of observations is small relative to the number of assets.
 *
 * Mathematical Background:
 * The shrinkage estimator is a convex combination:
 *
 *     Cov_shrunk = δ * F + (1-δ) * S
 *
 * where:
 * - S is the sample covariance matrix
 * - F is the target matrix (structured, low-variance estimate)
 * - δ is the shrinkage intensity (0 ≤ δ ≤ 1)
 *
 * Optimal Shrinkage:
 * The optimal δ minimizes the expected squared Frobenius norm of the
 * estimation error. Ledoit & Wolf provide an analytical formula for
 * δ that can be estimated from the data.
 *
 * Why Shrinkage Works:
 *
 * Problem: Sample covariance has high estimation error when T/N is small
 * - Large variance (sampling error)
 * - Extreme eigenvalues (overestimated largest, underestimated smallest)
 * - Poor out-of-sample performance
 *
 * Solution: Shrink towards a simple, structured target
 * - Reduces variance at the cost of introducing small bias
 * - Stabilizes eigenvalues
 * - Improves out-of-sample forecasting
 * - Better condition number for optimization
 *
 * Common Targets:
 *
 * 1. IDENTITY: F = I (assumes uncorrelated assets)
 *    - Simplest target
 *    - Maximum shrinkage effect
 *    - Use when you have little prior information
 *
 * 2. CONSTANT_VARIANCE: F = diag(σ²_avg) (common variance, no correlation)
 *    - Preserves information about average variance
 *    - Still assumes zero correlation
 *
 * 3. CONSTANT_CORRELATION: F based on constant correlation model
 *    - Most popular choice (Ledoit-Wolf 2004)
 *    - Assumes all pairs have same correlation
 *    - Preserves individual variances
 *    - Best theoretical properties
 *
 * 4. MARKET_MODEL: F based on single-factor (market) model
 *    - Assumes returns driven by one common factor
 *    - Appropriate for equity portfolios
 *    - Requires market index data
 *
 * When to Use Shrinkage:
 * - T/N < 10 (strong shrinkage)
 * - T/N < 5 (very strong shrinkage)
 * - Need stable, well-conditioned matrix
 * - Out-of-sample performance critical
 * - Portfolio optimization (mean-variance, risk parity)
 *
 * Performance: ~5ms for 100x100 matrix on modern CPU
 *
 * References:
 * - Ledoit & Wolf (2004), "Honey, I Shrunk the Sample Covariance Matrix"
 * - Ledoit & Wolf (2004), "A Well-Conditioned Estimator for Large-Dimensional
 *   Covariance Matrices"
 * - Schäfer & Strimmer (2005), "A Shrinkage Approach to Large-Scale Covariance
 *   Matrix Estimation"
 */

#pragma once

#include "risk/risk_model.hpp"

namespace portfolio
{
    namespace risk
    {

        /**
         * @enum ShrinkageTarget
         * @brief Target matrix structure for shrinkage estimation
         *
         * Defines the structured target matrix F towards which the sample
         * covariance is shrunk. Different targets make different assumptions
         * about the structure of the true covariance matrix.
         */
        enum class ShrinkageTarget
        {
            /**
             * @brief Identity matrix: F = I
             *
             * Assumes:
             * - All assets have unit variance
             * - All assets are uncorrelated
             *
             * Most aggressive shrinkage. Use when you have minimal
             * prior information about asset relationships.
             */
            IDENTITY,

            /**
             * @brief Diagonal with constant variance: F = σ²_avg * I
             *
             * Assumes:
             * - All assets have the same variance (average)
             * - All assets are uncorrelated
             *
             * Preserves information about overall variance level
             * while assuming no correlation structure.
             */
            CONSTANT_VARIANCE,

            /**
             * @brief Constant correlation model (recommended)
             *
             * Assumes:
             * - Assets have their individual variances (from sample)
             * - All pairs have the same correlation ρ_avg
             *
             * Formula: F_ij = σ_i * σ_j * ρ_avg for i≠j, F_ii = σ_i²
             *
             * This is the most popular choice and has the best
             * theoretical properties (Ledoit-Wolf 2004).
             */
            CONSTANT_CORRELATION,

            /**
             * @brief Single-factor (market) model
             *
             * Assumes:
             * - Returns driven by one common factor (market)
             * - Assets have factor loadings β_i
             * - Idiosyncratic risks are uncorrelated
             *
             * Formula: F_ij = β_i * β_j * σ²_market + δ_ij * σ²_idiosyncratic
             *
             * Appropriate for equity portfolios. Requires market index.
             * (Note: Current implementation uses simplified version)
             */
            MARKET_MODEL
        };

        /**
         * @class LedoitWolfShrinkage
         * @brief Ledoit-Wolf shrinkage estimator for covariance matrices
         *
         * Implements optimal shrinkage of the sample covariance matrix towards
         * a structured target. The shrinkage intensity is estimated analytically
         * to minimize expected estimation error.
         *
         * Key Benefits:
         * - Reduced estimation error (especially for small samples)
         * - Always positive definite (if target is positive definite)
         * - Better conditioned (smaller condition number)
         * - Improved out-of-sample performance
         * - Theoretical optimality guarantees
         *
         * Typical Shrinkage Intensities:
         * - T/N = 1.0 → δ ≈ 0.8-0.9 (heavy shrinkage)
         * - T/N = 2.0 → δ ≈ 0.5-0.7 (moderate shrinkage)
         * - T/N = 5.0 → δ ≈ 0.2-0.4 (mild shrinkage)
         * - T/N = 10+ → δ ≈ 0.0-0.2 (minimal shrinkage)
         *
         * Impact on Eigenvalues:
         * - Shrinkage pulls extreme eigenvalues towards the average
         * - Reduces eigenvalue spread (condition number)
         * - Makes matrix better suited for portfolio optimization
         *
         * Usage Example:
         * @code
         * // Standard Ledoit-Wolf with constant correlation target
         * LedoitWolfShrinkage lw(ShrinkageTarget::CONSTANT_CORRELATION);
         * Eigen::MatrixXd cov = lw.estimate_covariance(returns);
         *
         * // Check shrinkage intensity
         * double delta = lw.get_shrinkage_intensity();
         * std::cout << "Shrinkage: " << delta << " (0=no shrinkage, 1=full shrinkage)\n";
         *
         * // Compare condition numbers
         * SampleCovariance sample;
         * auto sample_cov = sample.estimate_covariance(returns);
         * std::cout << "Condition number improvement: "
         *           << cond(sample_cov) / cond(cov) << "x\n";
         *
         * // Use fixed shrinkage for testing/comparison
         * LedoitWolfShrinkage lw_fixed(ShrinkageTarget::IDENTITY, 0.5);
         * auto cov_fixed = lw_fixed.estimate_covariance(returns);
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations.
         * Note: Last shrinkage intensity is mutable state, updated on each call.
         */
        class LedoitWolfShrinkage : public RiskModel
        {
        public:
            /**
             * @brief Construct Ledoit-Wolf shrinkage estimator
             * @param target Shrinkage target type (default: constant correlation)
             * @param shrinkage_override Optional fixed shrinkage (0 ≤ δ ≤ 1)
             * @throws std::invalid_argument if shrinkage_override not in [0, 1]
             *
             * If shrinkage_override < 0 (default -1.0):
             * - Optimal shrinkage is computed automatically from data
             * - Recommended for production use
             *
             * If shrinkage_override in [0, 1]:
             * - Use the specified fixed shrinkage intensity
             * - Useful for testing, sensitivity analysis, or when you
             *   have external knowledge of the optimal shrinkage
             *
             * Target selection guide:
             * - CONSTANT_CORRELATION (default): Best all-around choice
             * - IDENTITY: Maximum shrinkage, minimal assumptions
             * - CONSTANT_VARIANCE: Middle ground between IDENTITY and CONSTANT_CORRELATION
             * - MARKET_MODEL: For equity portfolios with market factor
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
             * Algorithm:
             * 1. Validate input
             * 2. Compute sample covariance matrix S
             * 3. Compute target matrix F based on target type
             * 4. Estimate optimal shrinkage intensity δ (if not overridden)
             * 5. Compute shrunk covariance: Cov = δ*F + (1-δ)*S
             * 6. Ensure symmetry
             *
             * Time complexity: O(N^2 * T) for computation
             * Space complexity: O(N^2) for matrices
             *
             * The optimal shrinkage intensity minimizes:
             *     E[||Cov_shrunk - Cov_true||²_F]
             * where ||·||_F is the Frobenius norm.
             *
             * @note The shrinkage intensity is stored internally and can be
             *       retrieved using get_shrinkage_intensity()
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
             * Interpretation:
             * - δ = 0: No shrinkage, pure sample covariance
             * - δ = 0.5: Equal weight to sample and target
             * - δ = 1: Full shrinkage, pure target matrix
             *
             * Typical values:
             * - δ < 0.2: Sample covariance is quite good (large T/N)
             * - 0.2 ≤ δ < 0.5: Moderate shrinkage beneficial
             * - 0.5 ≤ δ < 0.8: Strong shrinkage needed (small T/N)
             * - δ ≥ 0.8: Very strong shrinkage (T ≈ N)
             *
             * @note Only valid after calling estimate_covariance()
             * @note Returns 0.0 if not yet estimated
             * @note Thread safety: Reading this value is not thread-safe
             *       if another thread is calling estimate_covariance()
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
             * @param returns Return matrix (T x N)
             * @param sample_cov Sample covariance matrix (N x N)
             * @return Target matrix F (N x N)
             *
             * Target computation depends on target type:
             * - IDENTITY: F = I
             * - CONSTANT_VARIANCE: F = (1/N) * tr(S) * I
             * - CONSTANT_CORRELATION: F_ij = σ_i * σ_j * ρ_avg
             * - MARKET_MODEL: Simplified single-factor model
             */
            Eigen::MatrixXd compute_target(
                const Eigen::MatrixXd &returns,
                const Eigen::MatrixXd &sample_cov) const;

            /**
             * @brief Compute optimal shrinkage intensity
             * @param returns Return matrix (T x N)
             * @param sample_cov Sample covariance (N x N)
             * @param target Target matrix (N x N)
             * @return Optimal shrinkage δ ∈ [0, 1]
             *
             * Uses Ledoit-Wolf analytical formula for optimal shrinkage.
             * The formula involves computing the variance of the sample
             * covariance estimator.
             *
             * If the computed value is outside [0,1], it is clamped to
             * that range for numerical stability.
             */
            double compute_shrinkage_intensity(
                const Eigen::MatrixXd &returns,
                const Eigen::MatrixXd &sample_cov,
                const Eigen::MatrixXd &target) const;

            /**
             * @brief Validate shrinkage override parameter
             * @param shrinkage Value to validate
             * @throws std::invalid_argument if shrinkage not in [-1, 1] or not in [0, 1] when >= 0
             */
            static void validate_shrinkage(double shrinkage);
        };

    } // namespace risk
} // namespace portfolio