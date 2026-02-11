#include <catch2/catch_test_macros.hpp>
#include "optimizer/mean_variance_optimizer.hpp"

using namespace portfolio::optimizer;

TEST_CASE("Schaible refinement converges on simple problem", "[Schaible]")
{
    Eigen::VectorXd mu(3);
    mu << 0.05, 0.10, 0.08; // expected returns

    Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(3, 3) * 0.02; // small positive-definite

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.sum_to_one = true;

    MeanVarianceOptimizer optimizer(ObjectiveType::MAX_SHARPE, 0.02);

    OptimizationResult result;

    REQUIRE_NOTHROW(result = optimizer.optimize_max_sharpe_direct(mu, cov, constraints));

    REQUIRE(result.success);
    REQUIRE(result.weights.size() == 3);
    REQUIRE(result.sharpe_ratio > 0.0);
    REQUIRE(result.message.find("Schaible") != std::string::npos);
}
