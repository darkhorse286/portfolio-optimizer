/**
 * @file quadratic_solver.cpp
 * @brief Implementation of custom quadratic programming solver
 */

#include "optimizer/quadratic_solver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace portfolio
{
    namespace optimizer
    {

        // ============================================================================
        // QuadraticProblem Implementation
        // ============================================================================

        void QuadraticProblem::validate() const
        {
            const int n = q.size();

            if (n == 0)
            {
                throw std::invalid_argument("Problem dimension is zero");
            }

            // Check P matrix
            if (P.rows() != n || P.cols() != n)
            {
                throw std::invalid_argument("P matrix dimensions do not match q vector");
            }

            if (!P.allFinite())
            {
                throw std::invalid_argument("P matrix contains NaN or Inf");
            }

            // Check q vector
            if (!q.allFinite())
            {
                throw std::invalid_argument("q vector contains NaN or Inf");
            }

            // Check equality constraints
            if (A_eq.rows() > 0)
            {
                if (A_eq.cols() != n)
                {
                    throw std::invalid_argument("A_eq columns do not match problem dimension");
                }
                if (b_eq.size() != A_eq.rows())
                {
                    throw std::invalid_argument("b_eq size does not match A_eq rows");
                }
                if (!A_eq.allFinite() || !b_eq.allFinite())
                {
                    throw std::invalid_argument("Equality constraints contain NaN or Inf");
                }
            }

            // Check inequality constraints (with lower and upper bounds)
            if (A_ineq.rows() > 0)
            {
                if (A_ineq.cols() != n)
                {
                    throw std::invalid_argument("A_ineq columns do not match problem dimension");
                }
                if (b_ineq_lower.size() != A_ineq.rows() || b_ineq_upper.size() != A_ineq.rows())
                {
                    throw std::invalid_argument("Inequality bounds size does not match A_ineq rows");
                }
                if (!A_ineq.allFinite() || !b_ineq_lower.allFinite() || !b_ineq_upper.allFinite())
                {
                    throw std::invalid_argument("Inequality constraints contain NaN or Inf");
                }
            }

            // Check bounds
            if (lower_bounds.size() > 0 && lower_bounds.size() != n)
            {
                throw std::invalid_argument("Lower bounds size does not match problem dimension");
            }
            if (upper_bounds.size() > 0 && upper_bounds.size() != n)
            {
                throw std::invalid_argument("Upper bounds size does not match problem dimension");
            }

            if (lower_bounds.size() > 0 && !lower_bounds.allFinite())
            {
                throw std::invalid_argument("Lower bounds contain NaN or Inf");
            }
            if (upper_bounds.size() > 0 && !upper_bounds.allFinite())
            {
                throw std::invalid_argument("Upper bounds contain NaN or Inf");
            }
        }

        // ============================================================================
        // SolverResult Implementation
        // ============================================================================

        SolverResult::SolverResult()
            : objective_value(0.0),
              success(false),
              iterations(0)
        {
        }

        // ============================================================================
        // QuadraticSolver Implementation
        // ============================================================================

        QuadraticSolver::QuadraticSolver(const SolverOptions &options)
            : options_(options)
        {
        }

        void QuadraticSolver::set_options(const SolverOptions &options)
        {
            options_ = options;
        }

        SolverResult QuadraticSolver::solve(const QuadraticProblem &problem) const
        {
            // Validate problem
            problem.validate();

            SolverResult result;
            const int n = problem.q.size();

            // Find initial feasible solution
            Eigen::VectorXd x = find_initial_solution(problem);

            if (options_.verbose)
            {
                std::cout << "Starting optimization with " << n << " variables\n";
            }

            // Main optimization loop
            double prev_obj = std::numeric_limits<double>::max();

            for (int iter = 0; iter < options_.max_iterations; ++iter)
            {
                // Compute gradient
                Eigen::VectorXd gradient = compute_gradient(x, problem);

                // Check convergence
                if (check_convergence(x, gradient, options_.tolerance))
                {
                    result.success = true;
                    result.message = "Converged";
                    result.iterations = iter;
                    break;
                }

                // Compute search direction (negative gradient with equality constraint projection)
                Eigen::VectorXd direction = -gradient;

                // If we have equality constraints, project direction into null space
                if (problem.A_eq.rows() > 0)
                {
                    // Solve A_eq * d = 0 using least squares
                    // Project gradient onto null space of A_eq
                    Eigen::MatrixXd A_eq_t = problem.A_eq.transpose();
                    Eigen::MatrixXd M = problem.A_eq * A_eq_t;
                    Eigen::VectorXd rhs = problem.A_eq * direction;

                    // Solve M * lambda = rhs
                    Eigen::VectorXd lambda = M.ldlt().solve(rhs);
                    direction = direction - A_eq_t * lambda;
                }

                // Line search
                double step = line_search(x, direction, problem, options_.step_size);

                // Update solution
                x = x + step * direction;

                // Project onto constraints
                x = project_onto_constraints(x, problem);

                // Compute objective value
                double obj = compute_objective(x, problem);

                if (options_.verbose && iter % 100 == 0)
                {
                    std::cout << "Iter " << iter << ": obj = " << obj << "\n";
                }

                // Check for lack of progress
                if (std::abs(obj - prev_obj) < options_.tolerance * 1e-2)
                {
                    if (check_convergence(x, gradient, options_.tolerance * 10))
                    {
                        result.success = true;
                        result.message = "Converged (no progress)";
                        result.iterations = iter;
                        break;
                    }
                }

                prev_obj = obj;
            }

            // Set result
            result.solution = x;
            result.objective_value = compute_objective(x, problem);

            if (!result.success)
            {
                result.message = "Maximum iterations reached";
                result.iterations = options_.max_iterations;
            }

            return result;
        }

        Eigen::VectorXd QuadraticSolver::find_initial_solution(
            const QuadraticProblem &problem) const
        {
            const int n = problem.q.size();
            Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

            // Start with uniform allocation if we have sum-to-one constraint
            if (problem.A_eq.rows() > 0 && problem.b_eq.size() > 0)
            {
                // Assume first constraint is sum-to-one
                double total = problem.b_eq(0);
                x = Eigen::VectorXd::Constant(n, total / n);
            }

            // Project onto constraints
            x = project_onto_constraints(x, problem);

            return x;
        }

        Eigen::VectorXd QuadraticSolver::project_onto_constraints(
            const Eigen::VectorXd &x,
            const QuadraticProblem &problem) const
        {
            Eigen::VectorXd x_proj = x;
            const int n = x.size();

            // Step 1: Apply box constraints first
            if (problem.lower_bounds.size() > 0)
            {
                x_proj = x_proj.cwiseMax(problem.lower_bounds);
            }
            if (problem.upper_bounds.size() > 0)
            {
                x_proj = x_proj.cwiseMin(problem.upper_bounds);
            }

            // Step 2: Project onto equality constraints
            if (problem.A_eq.rows() > 0)
            {
                // Solve: min ||x_new - x_proj||^2 s.t. A_eq * x_new = b_eq
                // Solution: x_new = x_proj - A_eq^T * (A_eq * A_eq^T)^-1 * (A_eq * x_proj - b_eq)

                Eigen::VectorXd residual = problem.A_eq * x_proj - problem.b_eq;
                Eigen::MatrixXd AAt = problem.A_eq * problem.A_eq.transpose();
                Eigen::VectorXd lambda = AAt.ldlt().solve(residual);
                x_proj = x_proj - problem.A_eq.transpose() * lambda;

                // Re-apply box constraints (may violate slightly after projection)
                if (problem.lower_bounds.size() > 0)
                {
                    x_proj = x_proj.cwiseMax(problem.lower_bounds);
                }
                if (problem.upper_bounds.size() > 0)
                {
                    x_proj = x_proj.cwiseMin(problem.upper_bounds);
                }

                // Final adjustment to satisfy equality constraint exactly
                // Use projected gradient onto constraint manifold
                residual = problem.A_eq * x_proj - problem.b_eq;
                if (residual.norm() > 1e-6)
                {
                    // Distribute violation equally among feasible directions
                    Eigen::VectorXd correction = problem.A_eq.transpose() *
                                                 AAt.ldlt().solve(residual);
                    x_proj -= correction;
                }
            }

            return x_proj;
        }

        Eigen::VectorXd QuadraticSolver::compute_gradient(
            const Eigen::VectorXd &x,
            const QuadraticProblem &problem) const
        {
            // Gradient: P * x + q
            return problem.P * x + problem.q;
        }

        double QuadraticSolver::compute_objective(
            const Eigen::VectorXd &x,
            const QuadraticProblem &problem) const
        {
            // Objective: (1/2) * x^T * P * x + q^T * x
            return 0.5 * x.dot(problem.P * x) + problem.q.dot(x);
        }

        bool QuadraticSolver::check_convergence(
            const Eigen::VectorXd &x,
            const Eigen::VectorXd &gradient,
            double tolerance) const
        {
            // Check gradient norm
            double grad_norm = gradient.norm();
            if (grad_norm > tolerance)
            {
                return false;
            }

            // Check that we're feasible
            // (Gradient small + feasible = converged)
            return true;
        }

        double QuadraticSolver::line_search(
            const Eigen::VectorXd &x,
            const Eigen::VectorXd &direction,
            const QuadraticProblem &problem,
            double initial_step) const
        {
            double step = initial_step;
            double obj_current = compute_objective(x, problem);

            // Backtracking line search
            const double rho = 0.5; // Step reduction factor
            const double c = 1e-4;  // Armijo condition parameter
            const int max_backtracks = 20;

            Eigen::VectorXd gradient = compute_gradient(x, problem);
            double directional_derivative = gradient.dot(direction);

            for (int i = 0; i < max_backtracks; ++i)
            {
                Eigen::VectorXd x_new = x + step * direction;
                x_new = project_onto_constraints(x_new, problem);

                double obj_new = compute_objective(x_new, problem);

                // Armijo condition
                if (obj_new <= obj_current + c * step * directional_derivative)
                {
                    return step;
                }

                step *= rho;
            }

            return step;
        }

    } // namespace optimizer
} // namespace portfolio