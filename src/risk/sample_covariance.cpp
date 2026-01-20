/**
 * @file sample_covariance.cpp
 * @brief Implementation of sample covariance estimator
 */

#include "risk/sample_covariance.hpp"

namespace portfolio
{
    namespace risk
    {

        SampleCovariance::SampleCovariance(bool bias_correction) : bias_correction_(bias_correction)
        {
            // No additional validation needed for boolean parameter
        }

        Eigen::MatrixXd SampleCovariance::estimate_covariance(const Eigen::MatrixXd &returns) const
        {

            // Validate input using base class method
            validate_returns(returns);

            const int n_obs = returns.rows();    // Number of observations (T)
            const int n_assets = returns.cols(); // Number of assets (N)

            // Compute column-wise mean (mean return for each asset)
            Eigen::RowVectorXd means = returns.colwise().mean();

            // Center the data: subtract mean from each column
            // This creates a temporary matrix but is more numerically stable
            Eigen::MatrixXd centered = returns.rowwise() - means;

            // Compute covariance: (1/k) * X^T * X where X is centered data
            // Using matrix multiplication for efficiency
            Eigen::MatrixXd covariance = (centered.transpose() * centered);

            // Apply normalization factor
            double normalization;
            if (bias_correction_)
            {
                // Bessel's correction: divide by (n-1) for unbiased estimate
                normalization = static_cast<double>(n_obs - 1);
            }
            else
            {
                // Maximum likelihood estimate: divide by n
                normalization = static_cast<double>(n_obs);
            }

            covariance /= normalization;

            // Ensure exact symmetry to eliminate numerical errors
            // This is important for downstream eigenvalue computations
            covariance = ensure_symmetric(covariance);

            return covariance;
        }

        std::string SampleCovariance::get_name() const
        {
            return "SampleCovariance";
        }

    } // namespace risk
} // namespace portfolio