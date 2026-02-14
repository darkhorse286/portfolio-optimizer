/**
 * @file osqp_solver.cpp
 * @brief Implementation of OSQP solver wrapper
 */

#include "optimizer/osqp_solver.hpp"
#include <iostream>
#include <cmath>

namespace portfolio
{
    namespace optimizer
    {

        OSQPSolver::OSQPSolver()
        {
            // Set default options
            options_.max_iterations = 10000;
            options_.tolerance = 1e-6;
            options_.step_size = 1.0; // Not used by OSQP, but kept for API consistency
        }

        void OSQPSolver::set_options(const SolverOptions &options)
        {
            options_ = options;
        }

        void OSQPSolver::convert_to_csc(
            const Eigen::MatrixXd &dense,
            std::vector<OSQPFloat> &data,
            std::vector<OSQPInt> &indices,
            std::vector<OSQPInt> &indptr,
            bool upper_triangular_only) const
        {
            const int rows = dense.rows();
            const int cols = dense.cols();

            data.clear();
            indices.clear();
            indptr.clear();
            indptr.reserve(cols + 1);

            indptr.push_back(0);

            // Iterate over columns (CSC format)
            for (int j = 0; j < cols; ++j)
            {
                int row_limit = upper_triangular_only ? (j + 1) : rows;

                for (int i = 0; i < row_limit; ++i)
                {
                    double val = dense(i, j);
                    if (std::abs(val) > 1e-14) // Skip near-zero values
                    {
                        data.push_back(val);
                        indices.push_back(i);
                    }
                }
                indptr.push_back(data.size());
            }
        }

        OSQPInt OSQPSolver::build_constraint_matrix(
            const QuadraticProblem &problem,
            std::vector<OSQPFloat> &A_data,
            std::vector<OSQPInt> &A_indices,
            std::vector<OSQPInt> &A_indptr,
            std::vector<OSQPFloat> &l,
            std::vector<OSQPFloat> &u) const
        {
            const int n = problem.q.size();
            const int n_eq = (problem.A_eq.size() > 0) ? problem.A_eq.rows() : 0;
            const int n_ineq = (problem.A_ineq.size() > 0) ? problem.A_ineq.rows() : 0;

            // Total constraints: n_eq equality + n_ineq inequality + n lower bounds + n upper bounds
            const int m = n_eq + n_ineq + n + n;

            A_data.clear();
            A_indices.clear();
            A_indptr.clear();
            l.resize(m);
            u.resize(m);

            A_indptr.push_back(0);

            // Build constraint matrix column by column (CSC format)
            for (int j = 0; j < n; ++j)
            {
                // 1. Equality constraints: A_eq * x = b_eq
                for (int i = 0; i < n_eq; ++i)
                {
                    double val = problem.A_eq(i, j);
                    if (std::abs(val) > 1e-14)
                    {
                        A_data.push_back(val);
                        A_indices.push_back(i);
                    }
                }

                // 2. Inequality constraints: A_ineq * x <= b_ineq_upper and >= b_ineq_lower
                for (int i = 0; i < n_ineq; ++i)
                {
                    double val = problem.A_ineq(i, j);
                    if (std::abs(val) > 1e-14)
                    {
                        A_data.push_back(val);
                        A_indices.push_back(n_eq + i);
                    }
                }

                // 3. Lower bound constraint row for x_j: x_j >= lb_j
                A_data.push_back(1.0);
                A_indices.push_back(n_eq + n_ineq + j);

                // 4. Upper bound constraint row for x_j: x_j <= ub_j
                A_data.push_back(1.0);
                A_indices.push_back(n_eq + n_ineq + n + j);

                A_indptr.push_back(A_data.size());
            }

            // Set constraint bounds
            // Equality constraints: l = u = b_eq
            for (int i = 0; i < n_eq; ++i)
            {
                l[i] = problem.b_eq(i);
                u[i] = problem.b_eq(i);
            }

            // Inequality constraints: l = b_ineq_lower, u = b_ineq_upper
            for (int i = 0; i < n_ineq; ++i)
            {
                l[n_eq + i] = (problem.b_ineq_lower.size() > 0) ? problem.b_ineq_lower(i) : -OSQP_INFTY;
                u[n_eq + i] = (problem.b_ineq_upper.size() > 0) ? problem.b_ineq_upper(i) : OSQP_INFTY;
            }

            // Lower bounds: x >= lb
            for (int i = 0; i < n; ++i)
            {
                l[n_eq + n_ineq + i] = problem.lower_bounds(i);
                u[n_eq + n_ineq + i] = OSQP_INFTY;
            }

            // Upper bounds: x <= ub
            for (int i = 0; i < n; ++i)
            {
                l[n_eq + n_ineq + n + i] = -OSQP_INFTY;
                u[n_eq + n_ineq + n + i] = problem.upper_bounds(i);
            }

            return m;
        }

        SolverResult OSQPSolver::solve(const QuadraticProblem &problem)
        {
            SolverResult result;
            const int n = problem.q.size();

            // Validate problem
            if (n == 0)
            {
                result.success = false;
                result.message = "Empty problem";
                return result;
            }

            if (problem.P.rows() != n || problem.P.cols() != n)
            {
                result.success = false;
                result.message = "P matrix size mismatch";
                return result;
            }

            // Convert P matrix to CSC format
            std::vector<OSQPFloat> P_data;
            std::vector<OSQPInt> P_indices;
            std::vector<OSQPInt> P_indptr;
            convert_to_csc(problem.P, P_data, P_indices, P_indptr, true);

            // Linear term q
            std::vector<OSQPFloat> q(n);
            for (int i = 0; i < n; ++i)
            {
                q[i] = problem.q(i);
            }

            // Build constraint matrix A and bounds l, u
            std::vector<OSQPFloat> A_data;
            std::vector<OSQPInt> A_indices;
            std::vector<OSQPInt> A_indptr;
            std::vector<OSQPFloat> l;
            std::vector<OSQPFloat> u;

            OSQPInt m = build_constraint_matrix(problem, A_data, A_indices, A_indptr, l, u);

            // Setup OSQP CSC matrices with CORRECT field order
            OSQPCscMatrix P_csc;
            P_csc.m = static_cast<OSQPInt>(n);
            P_csc.n = static_cast<OSQPInt>(n);
            P_csc.p = P_indptr.data();
            P_csc.i = P_indices.data();
            P_csc.x = P_data.data();
            P_csc.nzmax = static_cast<OSQPInt>(P_data.size());
            P_csc.nz = -1;  // -1 means CSC format (not triplet)

            OSQPCscMatrix A_csc;
            A_csc.m = m;
            A_csc.n = static_cast<OSQPInt>(n);
            A_csc.p = A_indptr.data();
            A_csc.i = A_indices.data();
            A_csc.x = A_data.data();
            A_csc.nzmax = static_cast<OSQPInt>(A_data.size());
            A_csc.nz = -1;

            // Configure OSQP settings
            OSQPSettings* settings = (OSQPSettings*)malloc(sizeof(OSQPSettings));
            osqp_set_default_settings(settings);
            settings->verbose = 0;
            settings->eps_abs = options_.tolerance;
            settings->eps_rel = options_.tolerance;
            settings->max_iter = options_.max_iterations;
            settings->polishing = 1;  // Note: 'polishing' not 'polish'

            // Setup OSQP workspace (use :: to disambiguate from our class)
            ::OSQPSolver* osqp_work = nullptr;
            OSQPInt exit_flag = osqp_setup(&osqp_work, &P_csc, q.data(), &A_csc,
                                           l.data(), u.data(), m, n, settings);

            free(settings);  // Free settings after setup

            if (exit_flag != 0)
            {
                result.success = false;
                result.message = "OSQP setup failed";
                return result;
            }

            // Solve the problem
            osqp_solve(osqp_work);

            // Extract solution
            result.solution = Eigen::VectorXd::Zero(n);
            for (int i = 0; i < n; ++i)
            {
                result.solution(i) = osqp_work->solution->x[i];
            }

            result.iterations = osqp_work->info->iter;
            result.objective_value = osqp_work->info->obj_val;

            // Check convergence status
            if (osqp_work->info->status_val == OSQP_SOLVED ||
                osqp_work->info->status_val == OSQP_SOLVED_INACCURATE)
            {
                result.success = true;
                result.message = osqp_work->info->status;
            }
            else
            {
                result.success = false;
                result.message = std::string(osqp_work->info->status);
            }

            // Cleanup
            osqp_cleanup(osqp_work);

            return result;
        }

    } // namespace optimizer
} // namespace portfolio