#include <catch2/catch_test_macros.hpp>
#include "optimizer/osqp_solver.hpp"
#include <cmath>

using namespace portfolio;

TEST_CASE("OSQP inequality constraint support", "[OSQP][Critical]") {
    // 2 variables: x, y
    // Objective: min x^2 + y^2
    // Constraint: x + y <= 0.5 (inequality)
    // Expected: x = y = 0.25 (on boundary)

    portfolio::optimizer::QuadraticProblem problem;
    // Small quadratic regularizer to choose center solution among maximizers
    problem.P = 2e-6 * Eigen::MatrixXd::Identity(2, 2);
    problem.q = -Eigen::VectorXd::Ones(2);

    // Inequality: x + y <= 0.5
    problem.A_ineq = Eigen::MatrixXd::Ones(1, 2);  // [1, 1]
    problem.b_ineq_lower = Eigen::VectorXd::Constant(1, -1e10);  // No lower bound
    problem.b_ineq_upper = Eigen::VectorXd::Constant(1, 0.5);

    problem.lower_bounds = Eigen::VectorXd::Zero(2);
    problem.upper_bounds = Eigen::VectorXd::Constant(2, 1.0);

    portfolio::optimizer::OSQPSolver solver;
    auto result = solver.solve(problem);

    REQUIRE(result.success);
    REQUIRE(std::abs(result.solution.sum() - 0.5) < 1e-6);
    REQUIRE(std::abs(result.solution(0) - 0.25) < 1e-6);
}
