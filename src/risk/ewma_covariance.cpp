/**
 * @file ewma_covariance.cpp
 * @brief Implementation of EWMA covariance estimator
 */

#include "risk/ewma_covariance.hpp"
#include <cmath>

namespace portfolio
{
    namespace risk
    {

        EWMACovariance::EWMACovariance(double lambda)
            : lambda_(lambda)
        {
            validate_lambda(lambda);
        }

        void EWMACovariance::validate_lambda(double lambda)
        {
            if (lambda <= 0.0 || lambda >= 1.0)
            {
                throw std::invalid_argument(
                    "Lambda must be in the range (0, 1), got: " +
                    std::to_string(lambda));
            }
        }

        Eigen::MatrixXd EWMACovariance::estimate_covariance(const Eigen::MatrixXd &returns) const
        {

            // Validate input
            validate_returns(returns);

            const int n_obs = returns.rows();    // Number of observations (T)
            const int n_assets = returns.cols(); // Number of assets (N)

            // Demean returns for better numerical stability
            // This is standard practice in EWMA estimation
            Eigen::RowVectorXd means = returns.colwise().mean();
            Eigen::MatrixXd centered = returns.rowwise() - means;

            // Initialize EWMA covariance with the first observation's outer product
            // This provides a reasonable starting point
            Eigen::VectorXd first_return = centered.row(0).transpose();
            Eigen::MatrixXd ewma_cov = first_return * first_return.transpose();

            // Recursive EWMA update: Cov_t = λ * Cov_{t-1} + (1-λ) * r_t * r_t^T
            // Start from the second observation
            for (int t = 1; t < n_obs; ++t)
            {
                Eigen::VectorXd return_t = centered.row(t).transpose();
                Eigen::MatrixXd outer_product = return_t * return_t.transpose();

                // Update EWMA covariance
                ewma_cov = lambda_ * ewma_cov + (1.0 - lambda_) * outer_product;
            }

            // Ensure exact symmetry
            ewma_cov = ensure_symmetric(ewma_cov);

            return ewma_cov;
        }

        std::string EWMACovariance::get_name() const
        {
            return "EWMACovariance";
        }

        double EWMACovariance::get_weight(int i) const
        {
            if (i < 0)
            {
                throw std::invalid_argument(
                    "Lag must be non-negative, got: " + std::to_string(i));
            }

            // Weight for observation i periods ago: w_i = (1-λ) * λ^i
            return (1.0 - lambda_) * std::pow(lambda_, static_cast<double>(i));
        }

    } // namespace risk
} // namespace portfolio