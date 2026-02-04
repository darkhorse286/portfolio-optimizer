/**
 * @file test_turnover_constraints.cpp
 * @brief Unit tests for per-asset turnover constraints in MeanVarianceOptimizer
 *
 * Validates that the per-asset turnover constraint correctly tightens box
 * bounds based on current_weights and max_turnover, and that all three
 * objective types (MIN_VARIANCE, TARGET_RETURN, RISK_AVERSION) enforce
 * the constraint consistently.
 *
 * Constraint semantics:
 *   |w_i - current_weights_i| <= max_turnover   for each asset i
 *
 * Implementation: bounds intersection in compute_turnover_bounds().
 * Box and turnover are intersected independently per asset; the effective
 * bound is the tighter of the two.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>

#include "optimizer/mean_variance_optimizer.hpp"

using namespace portfolio;
using namespace portfolio::optimizer;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Test Fixture
// ============================================================================

/**
 * @class TurnoverConstraintFixture
 * @brief Shared test data for turnover constraint tests
 *
 * Provides the same 2-asset and 10-asset market data used by the
 * mean-variance optimizer tests, plus pre-computed starting portfolios
 * used as current_weights across the test suite.
 */
class TurnoverConstraintFixture
{
protected:
    // --- Market data (identical to OptimizerTestFixture) ---
    Eigen::VectorXd returns_2asset_;
    Eigen::MatrixXd cov_2asset_;

    Eigen::VectorXd returns_10asset_;
    Eigen::MatrixXd cov_10asset_;

    // --- Starting portfolios for turnover tests ---
    Eigen::VectorXd equal_weights_2_;       ///< [0.5, 0.5]
    Eigen::VectorXd equal_weights_10_;      ///< [0.1, 0.1, ..., 0.1]
    Eigen::VectorXd concentrated_10_;       ///< [0.46, 0.06, 0.06, ..., 0.06]

    TurnoverConstraintFixture()
    {
        setup_2asset_data();
        setup_10asset_data();
        setup_starting_portfolios();
    }

    void setup_2asset_data()
    {
        returns_2asset_ = Eigen::VectorXd(2);
        returns_2asset_ << 0.10, 0.08;

        cov_2asset_ = Eigen::MatrixXd(2, 2);
        cov_2asset_ << 0.04, 0.01,
                       0.01, 0.02;
    }

    void setup_10asset_data()
    {
        returns_10asset_ = Eigen::VectorXd(10);
        returns_10asset_ << 0.10, 0.08, 0.12, 0.09, 0.11,
                            0.07, 0.13, 0.08, 0.10, 0.09;

        cov_10asset_ = Eigen::MatrixXd::Identity(10, 10) * 0.04;
        for (int i = 0; i < 10; ++i)
        {
            for (int j = i + 1; j < 10; ++j)
            {
                cov_10asset_(i, j) = 0.3 * 0.04;
                cov_10asset_(j, i) = 0.3 * 0.04;
            }
        }
    }

    void setup_starting_portfolios()
    {
        equal_weights_2_ = Eigen::VectorXd::Constant(2, 0.5);

        equal_weights_10_ = Eigen::VectorXd::Constant(10, 0.1);

        // Concentrated: asset 0 holds 46%, remaining 9 assets hold 6% each
        // Sum = 0.46 + 9 * 0.06 = 0.46 + 0.54 = 1.0
        concentrated_10_ = Eigen::VectorXd::Constant(10, 0.06);
        concentrated_10_(0) = 0.46;
    }
};

// ============================================================================
// Helper: compute actual per-asset turnover from two weight vectors
// ============================================================================

/**
 * @brief Compute maximum per-asset weight change between two portfolios
 * @param w_new New weights
 * @param w_old Previous weights
 * @return max_i |w_new_i - w_old_i|
 */
static double max_per_asset_turnover(const Eigen::VectorXd &w_new,
                                     const Eigen::VectorXd &w_old)
{
    return (w_new - w_old).cwiseAbs().maxCoeff();
}

// ============================================================================
// Test 1: Turnover disabled by default (no current_weights)
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Turnover disabled when current_weights not provided",
                 "[TurnoverConstraints][Default]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.05;  // Set a tight turnover limit

    // Omit current_weights entirely (uses default empty vector)
    auto result = opt.optimize(returns_2asset_, cov_2asset_, constraints);

    REQUIRE(result.success);
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

    // With no current_weights, turnover constraint must be inactive.
    // Weights are free across the full [0, 1] box, so the min-variance
    // solution should NOT be artificially pinned near [0.5, 0.5].
    // The 2-asset min-variance analytical solution heavily favors asset 1
    // (lower variance), so we expect a significant deviation from equal.
    REQUIRE(std::abs(result.weights(0) - 0.5) > 0.05);
}

// ============================================================================
// Test 2: Turnover limits weight movement - 2 asset
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Per-asset turnover limits weight movement",
                 "[TurnoverConstraints][Basic]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.10;  // Each weight can move at most 10%

    // Start from equal weights [0.5, 0.5]
    auto result = opt.optimize(returns_2asset_, cov_2asset_, constraints, equal_weights_2_);

    REQUIRE(result.success);
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

    // Each weight must stay within [0.4, 0.6] (equal +/- 0.10)
    for (int i = 0; i < 2; ++i)
    {
        INFO("Asset " << i << " weight: " << result.weights(i));
        REQUIRE(result.weights(i) >= equal_weights_2_(i) - 0.10 - 1e-6);
        REQUIRE(result.weights(i) <= equal_weights_2_(i) + 0.10 + 1e-6);
    }

    // Verify using the helper
    double actual_turnover = max_per_asset_turnover(result.weights, equal_weights_2_);
    INFO("Max per-asset turnover: " << actual_turnover);
    REQUIRE(actual_turnover <= 0.10 + 1e-6);
}

// ============================================================================
// Test 3: Zero turnover locks portfolio to current weights
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Zero turnover locks portfolio to current weights",
                 "[TurnoverConstraints][ZeroTurnover]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.0;  // No movement allowed

    // equal_weights_10_ sums to 1.0 and satisfies all box constraints,
    // so it is the unique feasible point.
    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

    REQUIRE(result.success);

    // Result must be (within solver tolerance) identical to starting weights
    for (int i = 0; i < 10; ++i)
    {
        INFO("Asset " << i << ": expected " << equal_weights_10_(i)
                            << ", got " << result.weights(i));
        REQUIRE_THAT(result.weights(i), WithinAbs(equal_weights_10_(i), 1e-4));
    }
}

// ============================================================================
// Test 4: Turnover from a concentrated starting portfolio
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Turnover from concentrated starting portfolio",
                 "[TurnoverConstraints][Concentrated]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.05;  // Tight: each asset moves at most 5%

    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, concentrated_10_);

    REQUIRE(result.success);
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

    // Every weight must stay within 5% of its starting value
    for (int i = 0; i < 10; ++i)
    {
        INFO("Asset " << i << ": start=" << concentrated_10_(i)
                            << " result=" << result.weights(i));
        REQUIRE(result.weights(i) >= concentrated_10_(i) - 0.05 - 1e-6);
        REQUIRE(result.weights(i) <= concentrated_10_(i) + 0.05 + 1e-6);
    }

    // Asset 0 started at 0.46 - it must remain >= 0.41
    REQUIRE(result.weights(0) >= 0.41 - 1e-6);
}

// ============================================================================
// Test 5: Box constraints tighter than turnover (box wins)
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Box constraints tighter than turnover are preserved",
                 "[TurnoverConstraints][ConstraintInteraction]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.08;   // Tight box: each asset in [8%, 12%]
    constraints.max_weight = 0.12;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.50; // Very loose turnover - should not bind

    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

    REQUIRE(result.success);
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

    // Box is [0.08, 0.12]; turnover from 0.1 is [0.1-0.5, 0.1+0.5] = [-0.4, 0.6].
    // Intersection is [0.08, 0.12]. All weights must be in box.
    for (int i = 0; i < 10; ++i)
    {
        INFO("Asset " << i << " weight: " << result.weights(i));
        REQUIRE(result.weights(i) >= 0.08 - 1e-6);
        REQUIRE(result.weights(i) <= 0.12 + 1e-6);
    }
}

// ============================================================================
// Test 6: Turnover tighter than box constraints (turnover wins)
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Turnover tighter than box constraints is enforced",
                 "[TurnoverConstraints][ConstraintInteraction]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;    // Very loose box
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.02; // Very tight turnover: each weight moves <= 2%

    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

    REQUIRE(result.success);
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

    // Box is [0, 1]; turnover from 0.1 is [0.08, 0.12].
    // Turnover wins. All weights must be in [0.08, 0.12].
    for (int i = 0; i < 10; ++i)
    {
        INFO("Asset " << i << " weight: " << result.weights(i));
        REQUIRE(result.weights(i) >= 0.08 - 1e-6);
        REQUIRE(result.weights(i) <= 0.12 + 1e-6);
    }

    double actual_turnover = max_per_asset_turnover(result.weights, equal_weights_10_);
    INFO("Max per-asset turnover: " << actual_turnover);
    REQUIRE(actual_turnover <= 0.02 + 1e-6);
}

// ============================================================================
// Test 7: current_weights size mismatch throws
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Mismatched current_weights size throws invalid_argument",
                 "[TurnoverConstraints][Validation]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.sum_to_one = true;
    constraints.max_turnover = 0.10;

    // 2-asset problem but 5-element current_weights
    Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant(5, 0.2);

    REQUIRE_THROWS_AS(
        opt.optimize(returns_2asset_, cov_2asset_, constraints, wrong_size),
        std::invalid_argument);
}

// ============================================================================
// Test 8: Negative max_turnover fails validation
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Negative max_turnover throws on constraint validation",
                 "[TurnoverConstraints][Validation]")
{
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.max_turnover = -0.1;  // Invalid

    REQUIRE_THROWS_AS(constraints.validate(), std::invalid_argument);
}

// ============================================================================
// Test 9: Turnover enforced consistently across all objective types
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Turnover constraint enforced across MIN_VARIANCE, TARGET_RETURN, RISK_AVERSION",
                 "[TurnoverConstraints][AllObjectives]")
{
    const double TURNOVER_LIMIT = 0.08;

    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    constraints.max_turnover = TURNOVER_LIMIT;

    // --- MIN_VARIANCE ---
    SECTION("MIN_VARIANCE respects turnover")
    {
        MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

        REQUIRE(result.success);
        REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

        double turnover = max_per_asset_turnover(result.weights, equal_weights_10_);
        INFO("MIN_VARIANCE max turnover: " << turnover);
        REQUIRE(turnover <= TURNOVER_LIMIT + 1e-6);
    }

    // --- TARGET_RETURN ---
    SECTION("TARGET_RETURN respects turnover")
    {
        MeanVarianceOptimizer opt(ObjectiveType::TARGET_RETURN);

        // Target return must be achievable within turnover bounds.
        // Equal weights on 10-asset data: return = mean(returns) = 0.097.
        // With turnover=0.08 from equal, weights stay near 0.1 each,
        // so achievable returns are close to 0.097.
        opt.set_target_return(0.097);

        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

        REQUIRE(result.success);
        REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

        double turnover = max_per_asset_turnover(result.weights, equal_weights_10_);
        INFO("TARGET_RETURN max turnover: " << turnover);
        REQUIRE(turnover <= TURNOVER_LIMIT + 1e-6);
    }

    // --- RISK_AVERSION ---
    SECTION("RISK_AVERSION respects turnover")
    {
        MeanVarianceOptimizer opt(ObjectiveType::RISK_AVERSION);
        opt.set_risk_aversion(2.0);

        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);

        REQUIRE(result.success);
        REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));

        double turnover = max_per_asset_turnover(result.weights, equal_weights_10_);
        INFO("RISK_AVERSION max turnover: " << turnover);
        REQUIRE(turnover <= TURNOVER_LIMIT + 1e-6);
    }
}

// ============================================================================
// Test 10: Tighter turnover produces result closer to starting portfolio
// ============================================================================

TEST_CASE_METHOD(TurnoverConstraintFixture,
                 "Tighter turnover produces result closer to starting portfolio",
                 "[TurnoverConstraints][Monotonicity]")
{
    // Optimization should produce results that are monotonically closer
    // to the starting portfolio as max_turnover decreases. Run three
    // levels and verify the ordering.

    OptimizationConstraints base;
    base.min_weight = 0.0;
    base.max_weight = 1.0;
    base.long_only = true;
    base.sum_to_one = true;

    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);

    double turnovers[] = {0.20, 0.10, 0.03};
    double prev_distance = 1e10;  // Start impossibly large

    for (double t : turnovers)
    {
        OptimizationConstraints constraints = base;
        constraints.max_turnover = t;

        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints, equal_weights_10_);
        REQUIRE(result.success);

        // L2 distance from starting portfolio
        double distance = (result.weights - equal_weights_10_).norm();
        INFO("Turnover limit " << t << " -> L2 distance from start: " << distance);

        // Distance must be non-increasing as turnover tightens
        REQUIRE(distance <= prev_distance + 1e-6);
        prev_distance = distance;
    }
}