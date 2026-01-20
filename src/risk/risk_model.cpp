/**
 * @file risk_model.cpp
 * @brief Implementation of RiskModel base class utilities
 */

#include "risk/risk_model.hpp"
#include <cmath>
#include <limits>

namespace portfolio
{
    namespace risk
    {
        Eigen::MatrixXd RiskModel::estimate_correlation(const Eigen::MatrixXd &returns)
            const
        {
            // Default implementation: compute covariance then convert
            Eigen::MatrixXd covariance = estimate_covariance(returns);
            return covariance_to_correlation(covariance);
        }

        void RiskModel::validate_returns(const Eigen::MatrixXd &returns)
        {
            // Check if matrix is empty
            if (returns.rows() == 0 || returns.cols() == 0)
            {
                throw std::invalid_argument("Returns matrix cannot be empty.");
            }

            // Check minimum number of observations
            if (returns.rows() < 2)
            {
                throw std::invalid_argument("Not enough observations to compute risk model. Returns matrix must have at least 2 observations for covariance estimation. Received: " + std::to_string(returns.rows()));
            }

            // Check for NaN or Inf values
            if (!returns.allFinite())
            {
                throw std::invalid_argument("Returns matrix contains NaN or Inf values.");
            }
        }

        Eigen::MatrixXd RiskModel::covariance_to_correlation(const Eigen::MatrixXd &covariance)
        {
            const int n = covariance.rows();

            // Validate input
            if (n == 0 || covariance.cols() != n)
            {
                throw std::runtime_error("Covariance matrix must be square and non-empty");
            }

            // Extract standard deviations from diagonal
            Eigen::VectorXd std_devs = covariance.diagonal().array().sqrt();

            // Check for non-positive variances
            for (int i = 0; i < n; ++i)
            {
                if (std_devs(i) <= 0.0)
                {
                    throw std::runtime_error("Covariance matrix has non-positive diagonal element at index " + std::to_string(i) + std::to_string(covariance(i, i)));
                }
            }

            // Compute correlation matrix
            // corr(i,j) = cov(i,j) / (std(i) * std(j))
            Eigen::MatrixXd correlation = Eigen::MatrixXd::Zero(n, n);

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == j)
                    {
                        // Diagonal elements are exactly 1.0
                        correlation(i, j) = 1.0;
                    }

                    else
                    {
                        correlation(i, j) = covariance(i, j) / (std_devs(i) * std_devs(j));
                    }

                    // Clamp to [-1, 1] to handle numerical errors
                    if (correlation(i, j) > 1.0)
                    {
                        correlation(i, j) = 1.0;
                    }
                    else
                    {
                        if (correlation(i, j) < -1.0)
                        {
                            correlation(i, j) = -1.0;
                        }
                    }
                }
            }

            return correlation;
        }

        Eigen::MatrixXd RiskModel::ensure_symmetric(const Eigen::MatrixXd &matrix)
        {
            // Enforce exact symmetry: (M + M^T) / 2
            return 0.5 * (matrix + matrix.transpose());
        }
    } // namespace risk
} // namespace portfolio