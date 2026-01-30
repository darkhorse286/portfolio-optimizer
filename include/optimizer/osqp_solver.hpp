/**
 * @file osqp_solver.hpp
 * @brief OSQP-based quadratic programming solver
 *
 * Wraps the OSQP library for solving portfolio optimization problems.
 * OSQP (Operator Splitting Quadratic Program) is a modern, efficient
 * solver designed specifically for embedded and real-time applications.
 *
 * Problem formulation:
 *   minimize     (1/2) x^T Q x + q^T x
 *   subject to   A_eq x = b_eq    (equality constraints)
 *                lb <= x <= ub     (box constraints)
 *
 * Performance: Typically converges in 10-100 iterations for portfolio problems.
 * Much faster and more robust than custom projected gradient descent.
 */

#pragma once

#include "quadratic_solver.hpp" // Reuse problem/result structures
#include <Eigen/Dense>
#include <osqp/osqp.h>
#include <vector>
#include <memory>

namespace portfolio
{
    namespace optimizer
    {

        /**
         * @class OSQPSolver
         * @brief Quadratic programming solver using OSQP library
         *
         * Solves QP problems of the form:
         *   min (1/2) x^T Q x + q^T x
         *   s.t. A_eq x = b_eq, lb <= x <= ub
         *
         * Features:
         * - Handles box constraints + equality constraints efficiently
         * - Typically converges in 10-100 iterations
         * - Warm-start capable (for efficient rebalancing)
         * - Industry-standard solver (used by Bloomberg, etc.)
         *
         * Usage Example:
         * @code
         * OSQPSolver solver;
         * SolverOptions options;
         * options.max_iterations = 10000;
         * options.tolerance = 1e-6;
         * solver.set_options(options);
         *
         * QuadraticProblem problem = ...;
         * SolverResult result = solver.solve(problem);
         *
         * if (result.success) {
         *     std::cout << "Optimal weights: " << result.x.transpose() << "\n";
         *     std::cout << "Converged in " << result.iterations << " iterations\n";
         * }
         * @endcode
         */
        class OSQPSolver
        {
        public:
            /**
             * @brief Default constructor
             */
            OSQPSolver();

            /**
             * @brief Destructor
             */
            ~OSQPSolver() = default;

            /**
             * @brief Set solver options
             * @param options Solver configuration
             */
            void set_options(const SolverOptions &options);

            /**
             * @brief Solve quadratic programming problem
             * @param problem QP problem specification
             * @return Solution with status, iterations, and objective value
             *
             * Time complexity: O(n³) worst case, typically much better in practice
             * Space complexity: O(n²) for storing problem matrices
             */
            SolverResult solve(const QuadraticProblem &problem);

        private:
            SolverOptions options_; ///< Solver configuration

            /**
             * @brief Convert Eigen dense matrix to OSQP sparse CSC format
             * @param dense Dense matrix (Eigen)
             * @param P_data Output: non-zero values
             * @param P_indices Output: row indices
             * @param P_indptr Output: column pointers
             * @param upper_triangular_only Only store upper triangle (for symmetric matrices)
             *
             * OSQP requires matrices in Compressed Sparse Column (CSC) format.
             * For symmetric matrices (like covariance), only upper triangle needed.
             */
            void convert_to_csc(
                const Eigen::MatrixXd &dense,
                std::vector<OSQPFloat> &P_data,
                std::vector<OSQPInt> &P_indices,
                std::vector<OSQPInt> &P_indptr,
                bool upper_triangular_only = false) const;

            /**
             * @brief Build constraint matrix for OSQP
             * @param problem QP problem
             * @param A_data Output: constraint matrix values
             * @param A_indices Output: row indices
             * @param A_indptr Output: column pointers
             * @param l Output: lower bounds for constraints
             * @param u Output: upper bounds for constraints
             * @return Number of constraints (m)
             *
             * Constructs constraint matrix as:
             *   A = [A_eq; I; I]  (equality + lower bounds + upper bounds)
             *   l = [b_eq; lb; -inf]
             *   u = [b_eq; inf; ub]
             */
            OSQPInt build_constraint_matrix(
                const QuadraticProblem &problem,
                std::vector<OSQPFloat> &A_data,
                std::vector<OSQPInt> &A_indices,
                std::vector<OSQPInt> &A_indptr,
                std::vector<OSQPFloat> &l,
                std::vector<OSQPFloat> &u) const;
        };

    } // namespace optimizer
} // namespace portfolio