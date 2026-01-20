/**
 * @file ledoit_wolf_shrinkage.cpp
 * @brief Implementation of Ledoit-Wolf shrinkage estimator
 */

#include "risk/ledoit_wolf_shrinkage.hpp"
#include "risk/sample_covariance.hpp"
#include <cmath>
#include <algorithm>

namespace portfolio
{
    namespace risk
    {

        LedoitWolfShrinkage::LedoitWolfShrinkage(ShrinkageTarget target, double shrinkage_override): target_(target), shrinkage_override_(shrinkage_override), last_shrinkage_(0.0)
        {

            validate_shrinkage(shrinkage_override);
        }

        void LedoitWolfShrinkage::validate_shrinkage(double shrinkage)
        {
            if (shrinkage < -1.0 || shrinkage > 1.0)
            {
                throw std::invalid_argument(
                    "Shrinkage override must be in [0, 1] or -1 for auto, got: " +
                    std::to_string(shrinkage));
            }

            if (shrinkage >= 0.0 && shrinkage > 1.0)
            {
                throw std::invalid_argument(
                    "Shrinkage override must be in [0, 1], got: " +
                    std::to_string(shrinkage));
            }
        }

        Eigen::MatrixXd LedoitWolfShrinkage::estimate_covariance(const Eigen::MatrixXd &returns) const
        {

            // Validate input
            validate_returns(returns);

            const int n_obs = returns.rows();    // T
            const int n_assets = returns.cols(); // N

            // Step 1: Compute sample covariance matrix
            SampleCovariance sample_estimator(true); // With bias correction
            Eigen::MatrixXd sample_cov = sample_estimator.estimate_covariance(returns);

            // Step 2: Compute target matrix
            Eigen::MatrixXd target = compute_target(returns, sample_cov);

            // Step 3: Determine shrinkage intensity
            double shrinkage;
            if (shrinkage_override_ >= 0.0)
            {
                // Use fixed shrinkage
                shrinkage = shrinkage_override_;
            }
            else
            {
                // Compute optimal shrinkage
                shrinkage = compute_shrinkage_intensity(returns, sample_cov, target);
            }

            // Store for later retrieval
            last_shrinkage_ = shrinkage;

            // Step 4: Compute shrunk covariance: δ*F + (1-δ)*S
            Eigen::MatrixXd shrunk_cov = shrinkage * target + (1.0 - shrinkage) * sample_cov;

            // Ensure exact symmetry
            shrunk_cov = ensure_symmetric(shrunk_cov);

            return shrunk_cov;
        }

        std::string LedoitWolfShrinkage::get_name() const
        {
            return "LedoitWolfShrinkage";
        }

        Eigen::MatrixXd LedoitWolfShrinkage::compute_target(const Eigen::MatrixXd &returns, const Eigen::MatrixXd &sample_cov) const
        {

            const int n = sample_cov.rows();
            Eigen::MatrixXd target(n, n);

            switch (target_)
            {
            case ShrinkageTarget::IDENTITY:
            {
                // F = I (identity matrix)
                target = Eigen::MatrixXd::Identity(n, n);
                break;
            }

            case ShrinkageTarget::CONSTANT_VARIANCE:
            {
                // F = (tr(S)/N) * I (average variance on diagonal)
                double avg_variance = sample_cov.trace() / n;
                target = avg_variance * Eigen::MatrixXd::Identity(n, n);
                break;
            }

            case ShrinkageTarget::CONSTANT_CORRELATION:
            {
                // F_ij = σ_i * σ_j * ρ_avg for i≠j, F_ii = σ_i²

                // Extract standard deviations
                Eigen::VectorXd std_devs = sample_cov.diagonal().array().sqrt();

                // Compute average correlation
                double sum_corr = 0.0;
                int count = 0;
                for (int i = 0; i < n; ++i)
                {
                    for (int j = i + 1; j < n; ++j)
                    {
                        double corr = sample_cov(i, j) / (std_devs(i) * std_devs(j));
                        sum_corr += corr;
                        count++;
                    }
                }
                double avg_corr = (count > 0) ? sum_corr / count : 0.0;

                // Construct target with constant correlation
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        if (i == j)
                        {
                            target(i, j) = sample_cov(i, i); // Keep individual variances
                        }
                        else
                        {
                            target(i, j) = std_devs(i) * std_devs(j) * avg_corr;
                        }
                    }
                }
                break;
            }

            case ShrinkageTarget::MARKET_MODEL:
            {
                // Simplified single-factor model
                // More sophisticated implementation would use actual market returns

                // Compute average pairwise covariance
                double sum_cov = 0.0;
                int count = 0;
                for (int i = 0; i < n; ++i)
                {
                    for (int j = i + 1; j < n; ++j)
                    {
                        sum_cov += sample_cov(i, j);
                        count++;
                    }
                }
                double avg_cov = (count > 0) ? sum_cov / count : 0.0;

                // Construct target with constant covariance
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        if (i == j)
                        {
                            target(i, j) = sample_cov(i, i); // Keep individual variances
                        }
                        else
                        {
                            target(i, j) = avg_cov;
                        }
                    }
                }
                break;
            }

            default:
                throw std::runtime_error("Unknown shrinkage target type");
            }

            return target;
        }

        double LedoitWolfShrinkage::compute_shrinkage_intensity(const Eigen::MatrixXd &returns, const Eigen::MatrixXd &sample_cov, const Eigen::MatrixXd &target) const
        {

            const int n_obs = returns.rows();    // T
            const int n_assets = returns.cols(); // N

            // Demean returns
            Eigen::RowVectorXd means = returns.colwise().mean();
            Eigen::MatrixXd centered = returns.rowwise() - means;

            // Compute variance of sample covariance estimator
            // This is the key step in the Ledoit-Wolf formula

            // π-hat: Sum of asymptotic variances of sample covariances
            double pi_hat = 0.0;
            for (int t = 0; t < n_obs; ++t)
            {
                Eigen::VectorXd r_t = centered.row(t).transpose();
                Eigen::MatrixXd outer = r_t * r_t.transpose();
                Eigen::MatrixXd diff = outer - sample_cov;
                pi_hat += diff.squaredNorm();
            }
            pi_hat /= n_obs;

            // ρ-hat: Misspecification term (difference between sample and target)
            double rho_hat = (sample_cov - target).squaredNorm();

            // γ-hat: Normalization term
            double gamma_hat = rho_hat;

            // Compute optimal shrinkage intensity
            double shrinkage;
            if (gamma_hat > 1e-10)
            { // Avoid division by zero
                shrinkage = std::max(0.0, std::min(1.0, pi_hat / gamma_hat));
            }
            else
            {
                // If target equals sample covariance, no shrinkage needed
                shrinkage = 0.0;
            }

            return shrinkage;
        }

    } // namespace risk
} // namespace portfolio