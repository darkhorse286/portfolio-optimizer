/**
 * @file risk_model.hpp
 * @brief Abstract interface for covariance estimation methods
 *
 * Provides a common interface for different risk model estimation
 * techniques used in portfolio optimization. All risk models must
 * implement the estimate_covariance method.
 *
 * Thread Safety: Implementations are expected to be thread-safe for
 * read-only operations.
 */

#pragma once

#include <Eigen/Dense>
#include <string>
#include <memory>

namespace portfolio
{
    namespace risk
    {

        /**
         * @class RiskModel
         * @brief Abstract base class for risk model estimation
         *
         * Defines the interface for covariance matrix estimation. Concrete
         * implementations include sample covariance, EWMA, and shrinkage methods.
         *
         * Usage Example:
         * @code
         * auto risk_model = std::make_unique<SampleCovariance>(true);
         * Eigen::MatrixXd cov = risk_model->estimate_covariance(returns);
         * @endcode
         */
        class RiskModel
        {
        public:
            virtual ~RiskModel() = default;

            /**
             * @brief Estimate covariance matrix from return data
             * @param returns Matrix of returns (rows = observations, cols = assets)
             * @return Covariance matrix (n_assets x n_assets)
             * @throws std::invalid_argument if returns matrix is empty
             * @throws std::runtime_error if computation fails
             *
             * @note The returned matrix is guaranteed to be symmetric
             * @note Matrix may not be positive definite for all estimators
             */
            virtual Eigen::MatrixXd estimate_covariance(
                const Eigen::MatrixXd &returns) const = 0;

            /**
             * @brief Estimate correlation matrix from return data
             * @param returns Matrix of returns (rows = observations, cols = assets)
             * @return Correlation matrix (n_assets x n_assets)
             * @throws std::invalid_argument if returns matrix is empty
             *
             * Default implementation: Convert covariance to correlation
             */
            virtual Eigen::MatrixXd estimate_correlation(
                const Eigen::MatrixXd &returns) const;

            /**
             * @brief Get the name of the risk model
             * @return String identifier for the model type
             */
            virtual std::string get_name() const = 0;

        protected:
            /**
             * @brief Validate input returns matrix
             * @param returns Matrix to validate
             * @throws std::invalid_argument if validation fails
             */
            static void validate_returns(const Eigen::MatrixXd &returns);

            /**
             * @brief Convert covariance matrix to correlation matrix
             * @param covariance Input covariance matrix
             * @return Correlation matrix
             * @throws std::runtime_error if diagonal elements are non-positive
             */
            static Eigen::MatrixXd covariance_to_correlation(
                const Eigen::MatrixXd &covariance);
        };

    } // namespace risk
} // namespace portfolio