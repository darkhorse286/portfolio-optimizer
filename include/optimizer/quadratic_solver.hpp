/**
 * @file quadratic_solver.hpp
 * @brief Custom quadratic programming solver
 *
 * Implements a simple but effective quadratic programming solver for
 * portfolio optimization problems. Uses active set method with projections.
 *
 * Solves problems of the form:
 *
 * Minimize:     (1/2) * x^T * P * x + q^T * x
 * Subject to:   A_eq * x = b_eq      (equality constraints)
 *              b_lower <= A_ineq * x <= b_upper  (inequality constraints with bounds)
 *              l <= x <= u           (box constraints)
 *
 * Algorithm: Active set method with projected gradient descent
 *
 * Performance: Fast for small to medium problems (N < 200)
 * Accuracy: Within 1e-6 of analytical solutions
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <optional>

namespace portfolio
{
    namespace optimizer
    {

        /**
         * @struct QuadraticProblem
         * @brief Quadratic programming problem specification
         */
        struct QuadraticProblem
        {
            Eigen::MatrixXd P; ///< Quadratic term (N x N)
            Eigen::VectorXd q; ///< Linear term (N x 1)

            Eigen::MatrixXd A_eq; ///< Equality constraint matrix
            Eigen::VectorXd b_eq; ///< Equality constraint values

            Eigen::MatrixXd A_ineq; ///< Inequality constraint matrix
            Eigen::VectorXd b_ineq_lower; ///< Inequality lower bounds (for A_ineq * x)
            Eigen::VectorXd b_ineq_upper; ///< Inequality upper bounds (for A_ineq * x)

            Eigen::VectorXd lower_bounds; ///< Lower bounds
            Eigen::VectorXd upper_bounds; ///< Upper bounds

            /**
             * @brief Default constructor
             */
            QuadraticProblem() = default;

            /**
             * @brief Validate problem specification
             * @throws std::invalid_argument if problem is ill-formed
             */
            void validate() const;
        };

        /**
         * @struct SolverOptions
         * @brief Options for quadratic solver
         */
        struct SolverOptions
        {
            int max_iterations = 1000; ///< Maximum iterations
            double tolerance = 1e-6;   ///< Convergence tolerance
            double step_size = 0.1;    ///< Initial step size
            bool verbose = false;      ///< Print progress

            /**
             * @brief Default constructor
             */
            SolverOptions() = default;
        };

        /**
         * @struct SolverResult
         * @brief Result from quadratic solver
         */
        struct SolverResult
        {
            Eigen::VectorXd solution; ///< Optimal solution
            double objective_value;   ///< Final objective value
            bool success;             ///< Convergence achieved
            int iterations;           ///< Number of iterations
            std::string message;      ///< Status message

            /**
             * @brief Default constructor
             */
            SolverResult();
        };

        /**
         * @class QuadraticSolver
         * @brief Custom quadratic programming solver
         *
         * Implements active set method for solving quadratic programs.
         * Designed specifically for portfolio optimization problems.
         *
         * Algorithm Steps:
         * 1. Initialize with feasible solution
         * 2. Identify active constraints
         * 3. Compute search direction in null space
         * 4. Line search with constraint projections
         * 5. Update active set
         * 6. Repeat until convergence
         *
         * Usage Example:
         * @code
         * QuadraticProblem problem;
         * problem.P = covariance;
         * problem.q = -expected_returns;
         * problem.A_eq = Eigen::MatrixXd::Ones(1, n);
         * problem.b_eq = Eigen::VectorXd::Ones(1);
         *
         * QuadraticSolver solver;
         * auto result = solver.solve(problem);
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class QuadraticSolver
        {
        public:
            /**
             * @brief Constructor
             * @param options Solver options
             */
            explicit QuadraticSolver(const SolverOptions &options = SolverOptions());

            /**
             * @brief Destructor
             */
            ~QuadraticSolver() = default;

            /**
             * @brief Solve quadratic program
             * @param problem Problem specification
             * @return Solver result
             * @throws std::invalid_argument if problem is invalid
             * @throws std::runtime_error if solver fails
             */
            SolverResult solve(const QuadraticProblem &problem) const;

            /**
             * @brief Set solver options
             */
            void set_options(const SolverOptions &options);

            /**
             * @brief Get solver options
             */
            const SolverOptions &get_options() const { return options_; }

        private:
            SolverOptions options_; ///< Solver configuration

            /**
             * @brief Find initial feasible solution
             */
            Eigen::VectorXd find_initial_solution(
                const QuadraticProblem &problem) const;

            /**
             * @brief Project solution onto constraints
             */
            Eigen::VectorXd project_onto_constraints(
                const Eigen::VectorXd &x,
                const QuadraticProblem &problem) const;

            /**
             * @brief Compute gradient at current point
             */
            Eigen::VectorXd compute_gradient(
                const Eigen::VectorXd &x,
                const QuadraticProblem &problem) const;

            /**
             * @brief Compute objective value
             */
            double compute_objective(
                const Eigen::VectorXd &x,
                const QuadraticProblem &problem) const;

            /**
             * @brief Check convergence
             */
            bool check_convergence(
                const Eigen::VectorXd &x,
                const Eigen::VectorXd &gradient,
                double tolerance) const;

            /**
             * @brief Perform line search
             */
            double line_search(
                const Eigen::VectorXd &x,
                const Eigen::VectorXd &direction,
                const QuadraticProblem &problem,
                double initial_step) const;
        };

    } // namespace optimizer
} // namespace portfolio