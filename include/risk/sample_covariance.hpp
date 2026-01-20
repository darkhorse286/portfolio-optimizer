/**
 * @file sample_covariance.hpp
 * @brief Classical sample covariance estimator
 *
 * Implements the standard sample covariance matrix estimation with
 * optional bias correction (Bessel's correction). This is the most
 * basic covariance estimator and serves as a baseline.
 *
 * Formula (with bias correction):
 *     Cov = (1/(n-1)) * (X - mean(X))^T * (X - mean(X))
 *
 * Performance: ~2ms for 100x100 matrix on modern CPU
 */

#pragma once

#include "risk_model.hpp"

namespace portfolio
{
    namespace risk
    {

        /**
         * @class SampleCovariance
         * @brief Sample covariance matrix estimator
         *
         * Computes the classical sample covariance matrix from historical returns.
         * Supports bias correction via Bessel's correction (dividing by n-1 instead of n).
         *
         * Properties:
         * - Unbiased estimator (with bias_correction = true)
         * - Positive semi-definite (guaranteed)
         * - Sensitive to outliers
         * - Requires n > p (observations > assets) for non-singular matrix
         *
         * Usage Example:
         * @code
         * SampleCovariance estimator(true);  // With bias correction
         * Eigen::MatrixXd cov = estimator.estimate_covariance(returns);
         * std::cout << "Determinant: " << cov.determinant() << std::endl;
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class SampleCovariance : public RiskModel
        {
        public:
            /**
             * @brief Construct sample covariance estimator
             * @param bias_correction Apply Bessel's correction (divide by n-1 vs n)
             *
             * Default is true to match pandas/NumPy default behavior
             */
            explicit SampleCovariance(bool bias_correction = true);

            /**
             * @brief Destructor
             */
            ~SampleCovariance() override = default;

            /**
             * @brief Estimate covariance matrix
             * @param returns Matrix of returns (T x N: observations x assets)
             * @return Covariance matrix (N x N)
             * @throws std::invalid_argument if returns is empty or has < 2 observations
             *
             * Time complexity: O(N^2 * T)
             * Space complexity: O(N^2)
             */
            Eigen::MatrixXd estimate_covariance(
                const Eigen::MatrixXd &returns) const override;

            /**
             * @brief Get model name
             * @return "SampleCovariance"
             */
            std::string get_name() const override;

            /**
             * @brief Check if using bias correction
             * @return true if using Bessel's correction
             */
            bool uses_bias_correction() const { return bias_correction_; }

        private:
            bool bias_correction_; ///< Whether to apply Bessel's correction
        };

    } // namespace risk
} // namespace portfolio