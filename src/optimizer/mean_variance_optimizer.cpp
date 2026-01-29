/**
 * @file mean_variance_optimizer.cpp
 * @brief Implementation of mean-variance portfolio optimizer
 */

#include "optimizer/mean_variance_optimizer.hpp"
#include <cmath>
#include <iostream>

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
            std::string obj_str;
            switch (objective_)
            {
            case ObjectiveType::MIN_VARIANCE:
                obj_str = "min_variance";
                break;
            case ObjectiveType::MAX_SHARPE:
                obj_str = "max_sharpe";
                break;
            case ObjectiveType::TARGET_RETURN:
                obj_str = "target_return";
                break;
            case ObjectiveType::RISK_AVERSION:
                obj_str = "risk_aversion";
                break;
            }

            return nlohmann::json{
                {"objective", obj_str},
                {"risk_free_rate", risk_free_rate_},
                {"risk_aversion", risk_aversion_},
                {"target_return", target_return_}};
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

            // Equality constraint: sum(w) = 1
            if (constraints.sum_to_one)
            {
                problem.A_eq = Eigen::MatrixXd::Ones(1, n);
                problem.b_eq = Eigen::VectorXd::Ones(1);
            }

            // Box constraints
            problem.lower_bounds = Eigen::VectorXd::Constant(n, constraints.min_weight);
            problem.upper_bounds = Eigen::VectorXd::Constant(n, constraints.max_weight);

            // Solve
            QuadraticSolver solver;
            SolverOptions options;
            options.max_iterations = 10000; 
            options.tolerance = 1e-3;      
            options.step_size = 1.0;      
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
            std::vector<double> lambdas = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};

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

            // Box constraints
            problem.lower_bounds = Eigen::VectorXd::Constant(n, constraints.min_weight);
            problem.upper_bounds = Eigen::VectorXd::Constant(n, constraints.max_weight);

            // Solve
            QuadraticSolver solver;
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

            // Equality constraint: sum(w) = 1
            if (constraints.sum_to_one)
            {
                problem.A_eq = Eigen::MatrixXd::Ones(1, n);
                problem.b_eq = Eigen::VectorXd::Ones(1);
            }

            // Box constraints
            problem.lower_bounds = Eigen::VectorXd::Constant(n, constraints.min_weight);
            problem.upper_bounds = Eigen::VectorXd::Constant(n, constraints.max_weight);

            // Solve
            QuadraticSolver solver;
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

    } // namespace optimizer
} // namespace portfolio