/**
 * @file mean_variance_optimizer.cpp
 * @brief Implementation of mean-variance portfolio optimizer
 */

#include "optimizer/mean_variance_optimizer.hpp"
#include "optimizer/osqp_solver.hpp"
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>

namespace portfolio
{
    namespace optimizer
    {

        MeanVarianceOptimizer::MeanVarianceOptimizer(
            ObjectiveType objective,
            double risk_free_rate)
            : objective_(objective),
              risk_free_rate_(risk_free_rate),
              risk_aversion_(1.0),
              target_return_(0.0)
        {
            if (risk_free_rate < 0.0)
            {
                throw std::invalid_argument(
                    "Risk-free rate must be non-negative, got: " +
                    std::to_string(risk_free_rate));
            }
        }

        std::string MeanVarianceOptimizer::get_name() const
        {
            return "MeanVarianceOptimizer";
        }

        nlohmann::json MeanVarianceOptimizer::get_parameters() const
        {
            nlohmann::json params;
            params["optimizer_type"] = "MeanVariance";
            params["risk_free_rate"] = risk_free_rate_;
            params["risk_aversion"] = risk_aversion_;
            params["target_return"] = target_return_;

            std::string obj_str;
            switch (objective_)
            {
            case ObjectiveType::MIN_VARIANCE:
                obj_str = "MIN_VARIANCE";
                break;
            case ObjectiveType::MAX_SHARPE:
                obj_str = "MAX_SHARPE";
                break;
            case ObjectiveType::TARGET_RETURN:
                obj_str = "TARGET_RETURN";
                break;
            case ObjectiveType::RISK_AVERSION:
                obj_str = "RISK_AVERSION";
                break;
            }
            params["objective"] = obj_str;

            return params;
        }

        void MeanVarianceOptimizer::set_risk_aversion(double lambda)
        {
            if (lambda <= 0.0)
            {
                throw std::invalid_argument(
                    "Risk aversion must be positive, got: " +
                    std::to_string(lambda));
            }
            risk_aversion_ = lambda;
        }

        void MeanVarianceOptimizer::set_target_return(double target_return)
        {
            target_return_ = target_return;
        }

        void MeanVarianceOptimizer::validate_parameters() const
        {
            if (objective_ == ObjectiveType::RISK_AVERSION && risk_aversion_ <= 0.0)
            {
                throw std::runtime_error(
                    "Risk aversion must be set for RISK_AVERSION objective");
            }

            if (objective_ == ObjectiveType::TARGET_RETURN && target_return_ == 0.0)
            {
                throw std::runtime_error(
                    "Target return must be set for TARGET_RETURN objective");
            }
        }

        OptimizationResult MeanVarianceOptimizer::optimize(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            // Validate inputs
            validate_inputs(expected_returns, covariance);
            constraints.validate();
            validate_parameters();

            // Validate current_weights dimensions if provided
            if (current_weights.size() > 0 &&
                current_weights.size() != static_cast<Eigen::Index>(expected_returns.size()))
            {
                throw std::invalid_argument(
                    "current_weights size (" + std::to_string(current_weights.size()) +
                    ") does not match number of assets (" +
                    std::to_string(expected_returns.size()) + ")");
            }

            // Route to appropriate optimizer based on objective
            switch (objective_)
            {
            case ObjectiveType::MIN_VARIANCE:
                return optimize_min_variance(covariance, constraints, current_weights);

            case ObjectiveType::MAX_SHARPE:
                return optimize_max_sharpe(expected_returns, covariance, constraints, current_weights);

            case ObjectiveType::TARGET_RETURN:
                return optimize_target_return(expected_returns, covariance, constraints, current_weights);

            case ObjectiveType::RISK_AVERSION:
                return optimize_risk_aversion(expected_returns, covariance, constraints, current_weights);

            default:
                throw std::runtime_error("Unknown objective type");
            }
        }

        OptimizationResult MeanVarianceOptimizer::optimize_min_variance(
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            const int n = covariance.rows();

            // Set up quadratic problem
            // Minimize: (1/2) * w^T * Sigma * w
            // Subject to: sum(w) = 1, w_min <= w <= w_max
            QuadraticProblem problem;
            problem.P = covariance;
            problem.q = Eigen::VectorXd::Zero(n);

            // Tracking error (soft penalty) integration
            if (constraints.benchmark_weights.size() > 0)
            {
                if (constraints.benchmark_weights.size() != n)
                {
                    throw std::invalid_argument(
                        "benchmark_weights size mismatch: " + std::to_string(constraints.benchmark_weights.size()) +
                        " vs " + std::to_string(n));
                }

                double penalty = constraints.tracking_error_penalty;
                if (penalty <= 0.0)
                {
                    penalty = constraints.compute_tracking_error_penalty(covariance);
                }

                // QuadraticProblem uses objective (1/2) x^T P x + q^T x
                // We need to add lambda_te * w^T Sigma w -> P += 2 * lambda_te * Sigma
                problem.P += 2.0 * penalty * covariance;
                problem.q -= 2.0 * penalty * covariance * constraints.benchmark_weights;
            }

            // Equality constraint: sum(w) = 1
            if (constraints.sum_to_one)
            {
                problem.A_eq = Eigen::MatrixXd::Ones(1, n);
                problem.b_eq = Eigen::VectorXd::Ones(1);
            }

            // Box constraints intersected with per-asset turnover
            compute_turnover_bounds(n, constraints, current_weights,
                                    problem.lower_bounds, problem.upper_bounds);

            // Solve
            OSQPSolver solver;
            SolverOptions options;
            options.max_iterations = 10000;
            options.tolerance = 1e-6;
            solver.set_options(options);

            SolverResult solver_result = solver.solve(problem);

            // Build optimization result
            OptimizationResult result;
            result.weights = solver_result.solution;
            result.success = solver_result.success;
            result.message = solver_result.message;
            result.iterations = solver_result.iterations;
            result.objective_value = solver_result.objective_value;

            // Calculate portfolio statistics (no expected returns for min variance)
            result.expected_return = 0.0;
            double variance = result.weights.transpose() * covariance * result.weights;
            result.volatility = std::sqrt(std::max(0.0, variance));
            result.sharpe_ratio = 0.0;

            return result;
        }

        OptimizationResult MeanVarianceOptimizer::optimize_max_sharpe(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            const int n = expected_returns.size();

            // Max Sharpe is solved using a transformation
            // Original: max (mu^T * w - rf) / sqrt(w^T * Sigma * w)
            // Transform: Introduce y = w / sqrt(w^T * Sigma * w)
            // Then maximize: (mu - rf)^T * y subject to: y^T * Sigma * y = 1
            //
            // We use an approximation: optimize risk-aversion with adaptive lambda

            // Start with moderate risk aversion
            double best_sharpe = -std::numeric_limits<double>::max();
            OptimizationResult best_result;

            // Try different risk aversion levels
            std::vector<double> lambdas = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                           1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0,
                                           8.0, 9.0, 10.0, 15.0, 20.0};

            for (double lambda : lambdas)
            {
                // Temporarily set risk aversion
                double original_lambda = risk_aversion_;
                const_cast<MeanVarianceOptimizer *>(this)->risk_aversion_ = lambda;

                // Optimize with risk aversion
                OptimizationResult result = optimize_risk_aversion(
                    expected_returns, covariance, constraints, current_weights);

                // Restore original lambda
                const_cast<MeanVarianceOptimizer *>(this)->risk_aversion_ = original_lambda;

                // Check if this is better
                if (result.success && result.sharpe_ratio > best_sharpe)
                {
                    best_sharpe = result.sharpe_ratio;
                    best_result = result;
                }
            }

            if (best_sharpe > -std::numeric_limits<double>::max())
            {
                best_result.message = "Optimized via risk aversion grid search";
                return best_result;
            }

            // Fallback to single optimization with lambda = 1
            return optimize_risk_aversion(expected_returns, covariance, constraints, current_weights);
        }

        OptimizationResult MeanVarianceOptimizer::optimize_target_return(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            const int n = expected_returns.size();

            // Set up quadratic problem
            // Minimize: (1/2) * w^T * Sigma * w
            // Subject to: sum(w) = 1, mu^T * w = target_return, w_min <= w <= w_max
            QuadraticProblem problem;
            problem.P = covariance;
            problem.q = Eigen::VectorXd::Zero(n);

            // Tracking error (soft penalty) integration
            if (constraints.benchmark_weights.size() > 0)
            {
                if (constraints.benchmark_weights.size() != n)
                {
                    throw std::invalid_argument(
                        "benchmark_weights size mismatch: " + std::to_string(constraints.benchmark_weights.size()) +
                        " vs " + std::to_string(n));
                }

                double penalty = constraints.tracking_error_penalty;
                if (penalty <= 0.0)
                {
                    penalty = constraints.compute_tracking_error_penalty(covariance);
                }

                problem.P += 2.0 * penalty * covariance;
                problem.q -= 2.0 * penalty * covariance * constraints.benchmark_weights;
            }

            // Equality constraints
            int n_eq = constraints.sum_to_one ? 2 : 1;
            problem.A_eq = Eigen::MatrixXd(n_eq, n);
            problem.b_eq = Eigen::VectorXd(n_eq);

            int row = 0;
            if (constraints.sum_to_one)
            {
                problem.A_eq.row(row) = Eigen::RowVectorXd::Ones(n);
                problem.b_eq(row) = 1.0;
                row++;
            }

            // Target return constraint
            problem.A_eq.row(row) = expected_returns.transpose();
            problem.b_eq(row) = target_return_;

            // Box constraints intersected with per-asset turnover
            compute_turnover_bounds(n, constraints, current_weights,
                                    problem.lower_bounds, problem.upper_bounds);

            // Solve
            OSQPSolver solver;
            SolverOptions options;
            options.max_iterations = 10000;
            options.tolerance = 1e-3;
            options.step_size = 1.0;
            solver.set_options(options);

            SolverResult solver_result = solver.solve(problem);

            // Build result
            OptimizationResult result = calculate_statistics(
                solver_result.solution, expected_returns, covariance, risk_free_rate_);

            result.success = solver_result.success;
            result.message = solver_result.message;
            result.iterations = solver_result.iterations;

            return result;
        }

        OptimizationResult MeanVarianceOptimizer::optimize_risk_aversion(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            const int n = expected_returns.size();

            // Set up quadratic problem
            // Minimize: (1/2) * w^T * Sigma * w - (1/lambda) * mu^T * w
            // Subject to: sum(w) = 1, w_min <= w <= w_max
            QuadraticProblem problem;
            problem.P = covariance;
            problem.q = -expected_returns / risk_aversion_;

            // Tracking error (soft penalty) integration
            if (constraints.benchmark_weights.size() > 0)
            {
                if (constraints.benchmark_weights.size() != n)
                {
                    throw std::invalid_argument(
                        "benchmark_weights size mismatch: " + std::to_string(constraints.benchmark_weights.size()) +
                        " vs " + std::to_string(n));
                }

                double penalty = constraints.tracking_error_penalty;
                if (penalty <= 0.0)
                {
                    penalty = constraints.compute_tracking_error_penalty(covariance);
                }

                problem.P += 2.0 * penalty * covariance;
                problem.q -= 2.0 * penalty * covariance * constraints.benchmark_weights;
            }

            // Equality constraint: sum(w) = 1
            if (constraints.sum_to_one)
            {
                problem.A_eq = Eigen::MatrixXd::Ones(1, n);
                problem.b_eq = Eigen::VectorXd::Ones(1);
            }

            // Box constraints intersected with per-asset turnover
            compute_turnover_bounds(n, constraints, current_weights,
                                    problem.lower_bounds, problem.upper_bounds);

            // Solve
            OSQPSolver solver;
            SolverOptions options;
            options.max_iterations = 10000;
            options.tolerance = 1e-3;
            options.step_size = 1.0;
            solver.set_options(options);

            SolverResult solver_result = solver.solve(problem);

            // Build result
            OptimizationResult result = calculate_statistics(
                solver_result.solution, expected_returns, covariance, risk_free_rate_);

            result.success = solver_result.success;
            result.message = solver_result.message;
            result.iterations = solver_result.iterations;

            return result;
        }

        OptimizationResult MeanVarianceOptimizer::optimize_max_sharpe_direct(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights) const
        {
            // Validate inputs and constraints
            validate_inputs(expected_returns, covariance);
            constraints.validate();

            const int n = expected_returns.size();
            if (covariance.rows() != n || covariance.cols() != n)
            {
                throw std::invalid_argument(
                    "Covariance dimensions do not match expected_returns size");
            }
            if (current_weights.size() != 0 && current_weights.size() != n)
            {
                throw std::invalid_argument(
                    "current_weights size does not match number of assets");
            }

            // Determine risk-free rate (use default 0.02 when not set)
            double rf = (risk_free_rate_ == 0.0) ? 0.02 : risk_free_rate_;

            // If user requested grid search only, use existing implementation
            if (max_sharpe_method_ == MaxSharpeMethod::GRID_SEARCH)
            {
                OptimizationResult res = optimize_max_sharpe(
                    expected_returns, covariance, constraints, current_weights);
                res.message = "GRID_SEARCH: used fallback grid search";
                return res;
            }

            // Coarse grid search for initialization
            std::vector<double> lambdas = {0.1, 0.5, 1.0, 2.0, 5.0};
            double best_sharpe = -std::numeric_limits<double>::infinity();
            OptimizationResult best_grid_result;
            Eigen::VectorXd best_grid_weights = Eigen::VectorXd::Zero(n);

            for (double lambda : lambdas)
            {
                double original_lambda = risk_aversion_;
                const_cast<MeanVarianceOptimizer *>(this)->risk_aversion_ = lambda;
                OptimizationResult r = optimize_risk_aversion(
                    expected_returns, covariance, constraints, current_weights);
                const_cast<MeanVarianceOptimizer *>(this)->risk_aversion_ = original_lambda;

                if (r.success && r.sharpe_ratio > best_sharpe)
                {
                    best_sharpe = r.sharpe_ratio;
                    best_grid_result = r;
                    best_grid_weights = r.weights;
                }
            }

            // If no successful grid result, fallback
            if (best_sharpe == -std::numeric_limits<double>::infinity())
            {
                OptimizationResult res = optimize_max_sharpe(
                    expected_returns, covariance, constraints, current_weights);
                res.message = "Fallback: grid search (no valid coarse initializations)";
                return res;
            }

            // If method is GRID_SEARCH we already returned; for DIRECT or HYBRID
            OptimizationResult final_result;

            if (max_sharpe_method_ == MaxSharpeMethod::DIRECT ||
                max_sharpe_method_ == MaxSharpeMethod::HYBRID)
            {
                // Attempt Schaible refinement
                try
                {
                    // Initialize y = w / sqrt(w^T Sigma w)
                    double quad = (best_grid_weights.transpose() * covariance * best_grid_weights)(0, 0);
                    double denom = std::sqrt(std::max(0.0, quad));
                    if (denom <= 0.0)
                    {
                        throw std::runtime_error("Invalid scaling during Schaible init: zero variance");
                    }
                    Eigen::VectorXd y = best_grid_weights / denom;

                    const int max_iters = 20;
                    const double tol = 1e-6;
                    bool converged = false;

                    for (int it = 0; it < max_iters; ++it)
                    {
                        Eigen::VectorXd y_new = schaibles_iteration(y, expected_returns, covariance, constraints, rf);

                        if (check_schaibles_convergence(y, y_new, tol))
                        {
                            y = y_new;
                            converged = true;
                            final_result.iterations = it + 1;
                            break;
                        }

                        y = y_new;
                    }

                    // Normalize back to weights
                    Eigen::VectorXd weights;
                    if (constraints.sum_to_one)
                    {
                        double s = y.sum();
                        if (std::abs(s) < 1e-12)
                        {
                            throw std::runtime_error("Normalization failed: sum(y) is zero");
                        }
                        weights = y / s;
                    }
                    else
                    {
                        weights = y; // preserve sum equal to initial y sum
                    }

                    // Compute stats
                    final_result = calculate_statistics(weights, expected_returns, covariance, rf);
                    final_result.weights = weights;
                    final_result.success = true;
                    final_result.message = converged ? "Schaible refinement converged" : "Schaible refinement completed without convergence";
                }
                catch (const std::exception &e)
                {
                    // Fallback to grid search result
                    best_grid_result.message = std::string("Schaible refinement failed: ") + e.what() + ". Used coarse grid result.";
                    return best_grid_result;
                }
            }
            else
            {
                // Should not reach: safe fallback
                return best_grid_result;
            }

            return final_result;
        }

        // ============================================================================
        // Per-Asset Turnover Constraint Helper
        // ============================================================================

        void MeanVarianceOptimizer::compute_turnover_bounds(
            int n,
            const OptimizationConstraints &constraints,
            const Eigen::VectorXd &current_weights,
            Eigen::VectorXd &lower_bounds,
            Eigen::VectorXd &upper_bounds) const
        {
            // Start from standard box constraints
            lower_bounds = Eigen::VectorXd::Constant(n, constraints.min_weight);
            upper_bounds = Eigen::VectorXd::Constant(n, constraints.max_weight);

            // Tighten with per-asset turnover when applicable:
            //   current_weights must be provided (non-empty)
            //   max_turnover must be restrictive (< 1.0 is always non-binding
            //   for weights in [0, 1])
            if (current_weights.size() == n && constraints.max_turnover < 1.0)
            {
                for (int i = 0; i < n; ++i)
                {
                    double turnover_lb = current_weights(i) - constraints.max_turnover;
                    double turnover_ub = current_weights(i) + constraints.max_turnover;

                    // Intersection: tighter of box and turnover wins per asset
                    lower_bounds(i) = std::max(lower_bounds(i), turnover_lb);
                    upper_bounds(i) = std::min(upper_bounds(i), turnover_ub);
                }
            }
        }

        /**
         * @brief Single iteration of Schaible's fractional programming method
         *
         * Solves the QP subproblem:
         *   minimize_y  y^T Sigma y
         *   subject to  (mu - rf)^T y = 1
         *               sum(y) = sum(y_current)
         *               min_weight <= y_i <= max_weight
         *
         * Returns the raw y vector produced by the QP solver. Caller is
         * responsible for normalization and convergence checks.
         */
        Eigen::VectorXd MeanVarianceOptimizer::schaibles_iteration(
            const Eigen::VectorXd &y_current,
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            double risk_free_rate) const
        {
            const int n = covariance.rows();

            // Basic dimension checks
            if (covariance.cols() != n)
            {
                throw std::invalid_argument("Covariance matrix must be square");
            }
            if (expected_returns.size() != n)
            {
                throw std::invalid_argument(
                    "expected_returns size (" + std::to_string(expected_returns.size()) +
                    ") does not match covariance dimension (" + std::to_string(n) + ")");
            }
            if (y_current.size() != n)
            {
                throw std::invalid_argument(
                    "y_current size (" + std::to_string(y_current.size()) +
                    ") does not match number of assets (" + std::to_string(n) + ")");
            }

            // Set up quadratic problem: minimize y^T Sigma y
            QuadraticProblem problem;
            problem.P = 2.0 * covariance; // follow Schaible subproblem formulation
            problem.q = Eigen::VectorXd::Zero(n);

            // Integrate tracking error penalty (if any). Note: P already set to 2*Sigma
            if (constraints.benchmark_weights.size() > 0)
            {
                if (constraints.benchmark_weights.size() != n)
                {
                    throw std::invalid_argument(
                        "benchmark_weights size mismatch: " + std::to_string(constraints.benchmark_weights.size()) +
                        " vs " + std::to_string(n));
                }

                double penalty = constraints.tracking_error_penalty;
                if (penalty <= 0.0)
                {
                    penalty = constraints.compute_tracking_error_penalty(covariance);
                }

                // Here objective is y^T Sigma y represented by (1/2) y^T P y with P=2*Sigma
                // Adding lambda_te * y^T Sigma y requires adding 2*lambda_te*Sigma to P
                problem.P += 2.0 * penalty * covariance;
                problem.q -= 2.0 * penalty * covariance * constraints.benchmark_weights;
            }

            // Equality constraints: (mu - rf)^T y = 1; sum(y) = sum(y_current)
            problem.A_eq = Eigen::MatrixXd(2, n);
            problem.b_eq = Eigen::VectorXd(2);

            Eigen::VectorXd mu_minus_rf = expected_returns - Eigen::VectorXd::Constant(n, risk_free_rate);
            problem.A_eq.row(0) = mu_minus_rf.transpose();
            problem.b_eq(0) = 1.0;

            problem.A_eq.row(1) = Eigen::RowVectorXd::Ones(n);
            problem.b_eq(1) = y_current.sum();

            // Box constraints: use simple box [min_weight, max_weight] in y-space
            problem.lower_bounds = Eigen::VectorXd::Constant(n, constraints.min_weight);
            problem.upper_bounds = Eigen::VectorXd::Constant(n, constraints.max_weight);

            // Solve with OSQP
            OSQPSolver solver;
            SolverOptions options;
            options.max_iterations = 10000;
            options.tolerance = 1e-6;
            options.step_size = 1.0;
            solver.set_options(options);

            SolverResult solver_result = solver.solve(problem);

            if (!solver_result.success)
            {
                throw std::runtime_error(
                    "Schaible iteration QP failed: " + solver_result.message);
            }

            return solver_result.solution;
        }

        bool MeanVarianceOptimizer::check_schaibles_convergence(
            const Eigen::VectorXd &y_old,
            const Eigen::VectorXd &y_new,
            double tolerance) const
        {
            // Use L2 norm for relative change
            const double old_norm = y_old.norm();
            const double diff_norm = (y_new - y_old).norm();

            if (old_norm < 1e-10)
            {
                // Fall back to absolute tolerance when previous vector is (near) zero
                return diff_norm < tolerance;
            }

            return (diff_norm / old_norm) < tolerance;
        }

    } // namespace optimizer
} // namespace portfolio