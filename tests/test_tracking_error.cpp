#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "optimizer/mean_variance_optimizer.hpp"

using namespace portfolio::optimizer;
using Catch::Matchers::WithinAbs;

TEST_CASE("Tracking error: compute penalty heuristic", "[TrackingError][Penalty]")
{
    OptimizationConstraints c;
    c.max_tracking_error = 0.10; // 10%
    c.benchmark_weights = Eigen::VectorXd(2);
    c.benchmark_weights << 0.5, 0.5;

    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(2, 2) * 0.04; // variances 4%

    double penalty = c.compute_tracking_error_penalty(cov);

    // Expected: 1.0 / (max_te^2 * trace(Σ) / n)
    double expected = 1.0 / (c.max_tracking_error * c.max_tracking_error * (cov.trace() / 2.0));

    REQUIRE_THAT(penalty, WithinAbs(expected, 1e-12));
}

TEST_CASE("Tracking error: validate throws when max_tracking_error missing", "[TrackingError][Validation]")
{
    OptimizationConstraints c;
    c.benchmark_weights = Eigen::VectorXd(2);
    c.benchmark_weights << 1.0, 0.0;
    c.max_tracking_error = 0.0; // invalid per validation

    REQUIRE_THROWS_AS(c.validate(), std::invalid_argument);
}

TEST_CASE("Tracking error: warning allowed for benchmark sum != 1", "[TrackingError][Validation]")
{
    OptimizationConstraints c;
    c.benchmark_weights = Eigen::VectorXd(3);
    c.benchmark_weights << 0.3, 0.3, 0.3; // sums to 0.9
    c.max_tracking_error = 0.05; // valid

    // Should not throw (warning only)
    REQUIRE_NOTHROW(c.validate());
}

TEST_CASE("Min-variance optimization pulls toward benchmark with penalty", "[TrackingError][Integration]")
{
    // Two-asset identity covariance -> unpenalized min-variance is equal weights
    Eigen::VectorXd expected_returns = Eigen::VectorXd::Zero(2);
    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(2, 2);

    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;

    // Baseline (no penalty)
    auto base = opt.optimize(expected_returns, cov, constraints);
    REQUIRE(base.success);
    REQUIRE_THAT(base.weights.sum(), WithinAbs(1.0, 1e-6));
    // Expect approximately equal weights
    REQUIRE_THAT(base.weights(0), WithinAbs(0.5, 1e-6));
    REQUIRE_THAT(base.weights(1), WithinAbs(0.5, 1e-6));

    // Now add a strong tracking error penalty toward benchmark [1,0]
    constraints.benchmark_weights = Eigen::VectorXd(2);
    constraints.benchmark_weights << 1.0, 0.0;
    constraints.tracking_error_penalty = 1e3; // very large
    constraints.max_tracking_error = 0.1; // set so validation passes

    auto penalized = opt.optimize(expected_returns, cov, constraints);
    REQUIRE(penalized.success);
    REQUIRE_THAT(penalized.weights.sum(), WithinAbs(1.0, 1e-6));

    // With strong penalty, weight should be pulled close to benchmark (first asset)
    REQUIRE(penalized.weights(0) > 0.9);
    REQUIRE(penalized.weights(1) < 0.1);
}
