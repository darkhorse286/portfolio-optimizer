/**
 * @file optimizer_interface.cpp
 * @brief Implementation of optimizer interface and common structures
 */

#include "optimizer/optimizer_interface.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace portfolio
{
    namespace optimizer
    {

        // ============================================================================
        // OptimizationConstraints Implementation
        // ============================================================================

        void OptimizationConstraints::validate() const
        {
            if (min_weight < 0.0 && long_only)
            {
                throw std::invalid_argument(
                    "min_weight must be non-negative for long_only portfolios, got: " +
                    std::to_string(min_weight));
            }

            if (max_weight <= 0.0)
            {
                throw std::invalid_argument(
                    "max_weight must be positive, got: " +
                    std::to_string(max_weight));
            }

            if (min_weight > max_weight)
            {
                throw std::invalid_argument(
                    "min_weight (" + std::to_string(min_weight) +
                    ") cannot exceed max_weight (" + std::to_string(max_weight) + ")");
            }

            if (max_turnover < 0.0)
            {
                throw std::invalid_argument(
                    "max_turnover must be non-negative, got: " +
                    std::to_string(max_turnover));
            }

            // Check if constraints are feasible
            if (sum_to_one && min_weight > 1.0)
            {
                throw std::invalid_argument(
                    "Infeasible constraints: min_weight > 1.0 with sum_to_one constraint");
            }
        }

        OptimizationConstraints OptimizationConstraints::from_json(const nlohmann::json &j)
        {
            OptimizationConstraints constraints;

            constraints.min_weight = j.value("min_weight", 0.0);
            constraints.max_weight = j.value("max_weight", 1.0);
            constraints.long_only = j.value("long_only", true);
            constraints.sum_to_one = j.value("sum_to_one", true);
            constraints.max_turnover = j.value("max_turnover", 1.0);

            constraints.validate();
            return constraints;
        }

        // ============================================================================
        // OptimizationResult Implementation
        // ============================================================================

        OptimizationResult::OptimizationResult()
            : expected_return(0.0),
              volatility(0.0),
              sharpe_ratio(0.0),
              success(false),
              iterations(0),
              objective_value(0.0)
        {
        }

        bool OptimizationResult::is_valid() const
        {
            if (!success)
                return false;
            if (weights.size() == 0)
                return false;
            if (!weights.allFinite())
                return false;
            if (volatility < 0.0)
                return false;

            return true;
        }

        void OptimizationResult::print_summary() const
        {
            std::cout << "\n=== Optimization Result ===\n";
            std::cout << "Status: " << (success ? "SUCCESS" : "FAILED") << "\n";
            std::cout << "Message: " << message << "\n";
            std::cout << "Iterations: " << iterations << "\n";
            std::cout << std::string(50, '-') << "\n";

            if (success && weights.size() > 0)
            {
                std::cout << "Portfolio Statistics:\n";
                std::cout << "  Expected Return:  " << std::fixed << std::setprecision(4)
                          << expected_return * 100 << "%\n";
                std::cout << "  Volatility:       " << volatility * 100 << "%\n";
                std::cout << "  Sharpe Ratio:     " << std::setprecision(3)
                          << sharpe_ratio << "\n";
                std::cout << "  Objective Value:  " << std::setprecision(6)
                          << objective_value << "\n";

                std::cout << "\nWeight Statistics:\n";
                std::cout << "  Number of assets: " << weights.size() << "\n";
                std::cout << "  Sum of weights:   " << std::setprecision(4)
                          << weights.sum() << "\n";
                std::cout << "  Min weight:       " << weights.minCoeff() << "\n";
                std::cout << "  Max weight:       " << weights.maxCoeff() << "\n";

                // Count non-zero positions
                int non_zero = 0;
                for (int i = 0; i < weights.size(); ++i)
                {
                    if (std::abs(weights(i)) > 1e-6)
                    {
                        ++non_zero;
                    }
                }
                std::cout << "  Non-zero positions: " << non_zero << "\n";
            }

            std::cout << "===========================\n"
                      << std::endl;
        }

        // ============================================================================
        // OptimizerInterface Protected Methods
        // ============================================================================

        void OptimizerInterface::validate_inputs(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance)
        {
            // Check dimensions
            if (expected_returns.size() == 0)
            {
                throw std::invalid_argument("Expected returns vector is empty");
            }

            if (covariance.rows() == 0 || covariance.cols() == 0)
            {
                throw std::invalid_argument("Covariance matrix is empty");
            }

            // Check consistency
            if (expected_returns.size() != covariance.rows() ||
                expected_returns.size() != covariance.cols())
            {
                throw std::invalid_argument(
                    "Dimension mismatch: expected returns size (" +
                    std::to_string(expected_returns.size()) +
                    ") does not match covariance dimensions (" +
                    std::to_string(covariance.rows()) + "x" +
                    std::to_string(covariance.cols()) + ")");
            }

            // Check for NaN or Inf
            if (!expected_returns.allFinite())
            {
                throw std::invalid_argument(
                    "Expected returns contain NaN or Inf values");
            }

            if (!covariance.allFinite())
            {
                throw std::invalid_argument(
                    "Covariance matrix contains NaN or Inf values");
            }

            // Check symmetry of covariance
            double asymmetry = (covariance - covariance.transpose()).cwiseAbs().maxCoeff();
            if (asymmetry > 1e-8)
            {
                throw std::invalid_argument(
                    "Covariance matrix is not symmetric (max asymmetry: " +
                    std::to_string(asymmetry) + ")");
            }

            // Check positive semi-definiteness (via eigenvalues)
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
            double min_eigenvalue = solver.eigenvalues().minCoeff();
            if (min_eigenvalue < -1e-8)
            {
                throw std::invalid_argument(
                    "Covariance matrix is not positive semi-definite (min eigenvalue: " +
                    std::to_string(min_eigenvalue) + ")");
            }
        }

        bool OptimizerInterface::check_constraints(
            const Eigen::VectorXd &weights,
            const OptimizationConstraints &constraints,
            double tolerance)
        {
            const int n = weights.size();

            // Check box constraints
            for (int i = 0; i < n; ++i)
            {
                if (weights(i) < constraints.min_weight - tolerance)
                {
                    return false;
                }
                if (weights(i) > constraints.max_weight + tolerance)
                {
                    return false;
                }
                if (constraints.long_only && weights(i) < -tolerance)
                {
                    return false;
                }
            }

            // Check sum constraint
            if (constraints.sum_to_one)
            {
                double sum = weights.sum();
                if (std::abs(sum - 1.0) > tolerance)
                {
                    return false;
                }
            }

            return true;
        }

        OptimizationResult OptimizerInterface::calculate_statistics(
            const Eigen::VectorXd &weights,
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            double risk_free_rate)
        {
            OptimizationResult result;
            result.weights = weights;
            result.success = true;

            // Calculate expected return
            result.expected_return = weights.dot(expected_returns);

            // Calculate volatility
            double variance = weights.transpose() * covariance * weights;
            result.volatility = std::sqrt(std::max(0.0, variance));

            // Calculate Sharpe ratio
            if (result.volatility > 1e-10)
            {
                result.sharpe_ratio = (result.expected_return - risk_free_rate) / result.volatility;
            }
            else
            {
                result.sharpe_ratio = 0.0;
            }

            // Objective value (variance)
            result.objective_value = variance;

            return result;
        }

    } // namespace optimizer
} // namespace portfolio