#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "backtest/backtest_engine.hpp"
#include <Eigen/Dense>

using namespace portfolio::backtest;
using namespace portfolio::optimizer;
using Catch::Approx;

TEST_CASE("BacktestEngine::project_weights", "[project_weights]") {
    BacktestParams params;
    BacktestEngine engine(params);

    SECTION("Weights summing to more than 1.0 get renormalized") {
        Eigen::VectorXd raw(3);
        raw << 0.5, 0.5, 0.3;  // sums to 1.3

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        // Should be normalized to sum to 1.0
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
        // Ratios should be preserved
        REQUIRE(projected[0] == Approx(0.5 / 1.3).epsilon(1e-9));
        REQUIRE(projected[1] == Approx(0.5 / 1.3).epsilon(1e-9));
        REQUIRE(projected[2] == Approx(0.3 / 1.3).epsilon(1e-9));
    }

    SECTION("Single asset at 1.0 clips to max_weight and distributes remainder") {
        Eigen::VectorXd raw(3);
        raw << 1.0, 0.0, 0.0;

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 0.5;  // cap at 50%

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
        REQUIRE(projected[0] <= 0.5 + 1e-9);  // clipped to max_weight
    }

    SECTION("All-zero input falls back to equal weight") {
        Eigen::VectorXd raw(3);
        raw << 0.0, 0.0, 0.0;

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        // Should default to equal weight
        REQUIRE(projected[0] == Approx(1.0 / 3.0).epsilon(1e-9));
        REQUIRE(projected[1] == Approx(1.0 / 3.0).epsilon(1e-9));
        REQUIRE(projected[2] == Approx(1.0 / 3.0).epsilon(1e-9));
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
    }

    SECTION("Already-feasible weights pass through unchanged") {
        Eigen::VectorXd raw(3);
        raw << 0.2, 0.3, 0.5;  // already sums to 1.0

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        REQUIRE(projected[0] == Approx(0.2).epsilon(1e-9));
        REQUIRE(projected[1] == Approx(0.3).epsilon(1e-9));
        REQUIRE(projected[2] == Approx(0.5).epsilon(1e-9));
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
    }

    SECTION("Post-renorm clip handles edge case where renorm pushes weight above cap") {
        // Create a scenario where:
        // 1. Some weights are clipped to [0, max_weight]
        // 2. Some get zeroed because they're below min_weight
        // 3. After renorm, remaining weights exceed max_weight
        // 4. Second clip fixes this

        Eigen::VectorXd raw(4);
        raw << 0.6, 0.6, 0.0, 0.0;  // sums to 1.2

        OptimizationConstraints constraints;
        constraints.min_weight = 0.1;  // min 10%
        constraints.max_weight = 0.5;  // max 50%

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        // All weights should respect constraints
        for (int i = 0; i < projected.size(); ++i) {
            REQUIRE(projected[i] >= -1e-9);  // non-negative (or near-zero due to numerics)
            REQUIRE(projected[i] <= 0.5 + 1e-9);  // respect max_weight
            if (projected[i] > 1e-9)  // if non-zero, should be >= min_weight
                REQUIRE(projected[i] >= 0.1 - 1e-9);
        }
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
    }

    SECTION("Negative weights get clipped to zero") {
        Eigen::VectorXd raw(3);
        raw << -0.1, 0.6, 0.5;

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        REQUIRE(projected[0] >= 0.0);  // negative clipped to zero
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
    }

    SECTION("Weights below min_weight are zeroed out") {
        Eigen::VectorXd raw(4);
        raw << 0.25, 0.25, 0.25, 0.25;

        OptimizationConstraints constraints;
        constraints.min_weight = 0.3;  // 30% minimum
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        // All weights below 30% should be zeroed out
        for (int i = 0; i < projected.size(); ++i) {
            REQUIRE((projected[i] == 0.0 || projected[i] >= 0.3 - 1e-9));
        }
        REQUIRE(projected.sum() == Approx(1.0).epsilon(1e-9));
    }

    SECTION("Empty weight vector returns empty") {
        Eigen::VectorXd raw(0);

        OptimizationConstraints constraints;
        constraints.min_weight = 0.0;
        constraints.max_weight = 1.0;

        Eigen::VectorXd projected = engine.project_weights(raw, constraints);

        REQUIRE(projected.size() == 0);
    }
}
