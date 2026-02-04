/**
 * @file test_risk_models.cpp
 * @brief Comprehensive unit tests for risk model implementations
 *
 * Tests all risk model implementations against reference computations
 * and validates edge cases, error handling, and numerical accuracy.
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "risk/sample_covariance.hpp"
#include "risk/ewma_covariance.hpp"
#include "risk/ledoit_wolf_shrinkage.hpp"
#include "risk/risk_model_factory.hpp"

using namespace portfolio::risk;
using Catch::Matchers::WithinAbs;

/**
 * @class RiskModelTestFixture
 * @brief Test fixture providing common test data and utilities
 */
class RiskModelTestFixture
{
protected:
    // Small test matrix for basic validation
    Eigen::MatrixXd returns_small_;

    // Larger test matrix for realistic scenarios
    Eigen::MatrixXd returns_large_;

    // Expected values computed with NumPy for validation
    Eigen::MatrixXd expected_sample_cov_;

    RiskModelTestFixture()
    {
        setup_small_data();
        setup_large_data();
        setup_expected_values();
    }

    void setup_small_data()
    {
        // 5 observations, 3 assets
        returns_small_ = Eigen::MatrixXd(5, 3);
        returns_small_ << 0.01, 0.02, -0.01,
            0.02, -0.01, 0.01,
            -0.01, 0.01, 0.02,
            0.00, 0.01, -0.01,
            0.01, 0.00, 0.00;
    }

    void setup_large_data()
    {
        // Generate synthetic returns: 100 observations, 5 assets
        // Using known seed for reproducibility
        returns_large_ = Eigen::MatrixXd(100, 5);
        returns_large_.setRandom();
        returns_large_ *= 0.02; // Scale to realistic return magnitude
    }

    void setup_expected_values()
    {
        // Pre-computed sample covariance for small data (with bias correction)
        // These values match NumPy's np.cov(returns_small_, rowvar=False)
        expected_sample_cov_ = Eigen::MatrixXd(3, 3);
        expected_sample_cov_ << 0.000125, -0.000050, -0.000050,
            -0.000050, 0.000150, 0.000025,
            -0.000050, 0.000025, 0.000150;
    }

    // Utility: Compute condition number
    double compute_condition_number(const Eigen::MatrixXd &matrix) const
    {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
        auto singular_values = svd.singularValues();
        double max_sv = singular_values(0);
        double min_sv = singular_values(singular_values.size() - 1);
        return max_sv / std::max(min_sv, 1e-10); // Avoid division by zero
    }

    // Utility: Check if matrix is symmetric
    bool is_symmetric(const Eigen::MatrixXd &matrix, double tol = 1e-10) const
    {
        return (matrix - matrix.transpose()).cwiseAbs().maxCoeff() < tol;
    }

    // Utility: Check if matrix is positive semi-definite
    bool is_positive_semidefinite(const Eigen::MatrixXd &matrix, double tol = 1e-10) const
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
        return solver.eigenvalues().minCoeff() >= -tol;
    }
};

// ========================
// SampleCovariance Tests
// ========================

TEST_CASE_METHOD(RiskModelTestFixture, "SampleCovariance construction", "[RiskModel][SampleCovariance]")
{

    SECTION("Default construction uses bias correction")
    {
        SampleCovariance cov;
        REQUIRE(cov.uses_bias_correction() == true);
        REQUIRE(cov.get_name() == "SampleCovariance");
    }

    SECTION("Explicit bias correction")
    {
        SampleCovariance cov_biased(false);
        REQUIRE(cov_biased.uses_bias_correction() == false);

        SampleCovariance cov_unbiased(true);
        REQUIRE(cov_unbiased.uses_bias_correction() == true);
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "SampleCovariance basic estimation", "[RiskModel][SampleCovariance]")
{
    SampleCovariance cov(true);

    SECTION("Estimate covariance from small data")
    {
        auto result = cov.estimate_covariance(returns_small_);

        // Check dimensions
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 3);

        // Check symmetry
        REQUIRE(is_symmetric(result));

        // Check positive semi-definite
        REQUIRE(is_positive_semidefinite(result));

        // Check diagonal elements are positive
        for (int i = 0; i < result.rows(); ++i)
        {
            REQUIRE(result(i, i) > 0.0);
        }
    }

    SECTION("Match expected values")
    {
        auto result = cov.estimate_covariance(returns_small_);

        // Compare element-wise with tolerance
        for (int i = 0; i < result.rows(); ++i)
        {
            for (int j = 0; j < result.cols(); ++j)
            {
                REQUIRE_THAT(result(i, j),
                             WithinAbs(expected_sample_cov_(i, j), 1e-4));
            }
        }
    }

    SECTION("Bias correction effect")
    {
        SampleCovariance cov_biased(false);
        SampleCovariance cov_unbiased(true);

        auto biased = cov_biased.estimate_covariance(returns_small_);
        auto unbiased = cov_unbiased.estimate_covariance(returns_small_);

        // Unbiased should be larger by factor n/(n-1)
        int n = returns_small_.rows();
        double expected_ratio = static_cast<double>(n) / (n - 1);

        // Check all elements
        for (int i = 0; i < biased.rows(); ++i)
        {
            for (int j = 0; j < biased.cols(); ++j)
            {
                double ratio = unbiased(i, j) / biased(i, j);
                REQUIRE_THAT(ratio, WithinAbs(expected_ratio, 1e-6));
            }
        }
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "SampleCovariance error handling", "[RiskModel][SampleCovariance]")
{
    SampleCovariance cov;

    SECTION("Empty matrix throws")
    {
        Eigen::MatrixXd empty(0, 0);
        REQUIRE_THROWS_AS(cov.estimate_covariance(empty), std::invalid_argument);
    }

    SECTION("Single observation throws")
    {
        Eigen::MatrixXd single(1, 3);
        single << 0.01, 0.02, 0.03;
        REQUIRE_THROWS_AS(cov.estimate_covariance(single), std::invalid_argument);
    }

    SECTION("Matrix with NaN throws")
    {
        Eigen::MatrixXd with_nan = returns_small_;
        with_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
        REQUIRE_THROWS_AS(cov.estimate_covariance(with_nan), std::invalid_argument);
    }

    SECTION("Matrix with Inf throws")
    {
        Eigen::MatrixXd with_inf = returns_small_;
        with_inf(0, 0) = std::numeric_limits<double>::infinity();
        REQUIRE_THROWS_AS(cov.estimate_covariance(with_inf), std::invalid_argument);
    }
}

// =====================
// EWMACovariance Tests
// =====================

TEST_CASE_METHOD(RiskModelTestFixture, "EWMACovariance construction", "[RiskModel][EWMACovariance]")
{

    SECTION("Valid lambda values")
    {
        EWMACovariance ewma1(0.94);
        REQUIRE_THAT(ewma1.get_lambda(), WithinAbs(0.94, 1e-10));

        EWMACovariance ewma2(0.5);
        REQUIRE_THAT(ewma2.get_lambda(), WithinAbs(0.5, 1e-10));

        EWMACovariance ewma3(0.99);
        REQUIRE_THAT(ewma3.get_lambda(), WithinAbs(0.99, 1e-10));
    }

    SECTION("Invalid lambda throws")
    {
        REQUIRE_THROWS_AS(EWMACovariance(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(1.0), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(-0.5), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(1.5), std::invalid_argument);
    }

    SECTION("Effective window calculation")
    {
        EWMACovariance ewma_94(0.94);
        REQUIRE_THAT(ewma_94.get_effective_window(), WithinAbs(16.67, 0.01));

        EWMACovariance ewma_97(0.97);
        REQUIRE_THAT(ewma_97.get_effective_window(), WithinAbs(33.33, 0.01));

        EWMACovariance ewma_99(0.99);
        REQUIRE_THAT(ewma_99.get_effective_window(), WithinAbs(100.0, 0.01));
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "EWMACovariance basic estimation", "[RiskModel][EWMACovariance]")
{
    EWMACovariance ewma(0.94);

    SECTION("Estimate covariance")
    {
        auto result = ewma.estimate_covariance(returns_large_);

        // Check dimensions
        REQUIRE(result.rows() == 5);
        REQUIRE(result.cols() == 5);

        // Check symmetry
        REQUIRE(is_symmetric(result));

        // Check positive semi-definite
        REQUIRE(is_positive_semidefinite(result));
    }

    SECTION("Weight function")
    {
        EWMACovariance ewma(0.94);

        // Weights should decrease exponentially
        double w0 = ewma.get_weight(0);
        double w1 = ewma.get_weight(1);
        double w10 = ewma.get_weight(10);

        REQUIRE(w0 > w1);
        REQUIRE(w1 > w10);

        // Check formula: w_i = (1-λ) * λ^i
        REQUIRE_THAT(w0, WithinAbs(0.06, 1e-10));
        REQUIRE_THAT(w1, WithinAbs(0.06 * 0.94, 1e-10));
    }

    SECTION("EWMA adapts to recent data")
    {
        // Create returns with increasing volatility
        Eigen::MatrixXd increasing_vol(100, 2);
        for (int i = 0; i < 100; ++i)
        {
            double scale = 1.0 + i / 50.0; // Volatility triples
            increasing_vol(i, 0) = scale * 0.01;
            increasing_vol(i, 1) = scale * 0.01;
        }

        EWMACovariance ewma(0.94);
        auto ewma_cov = ewma.estimate_covariance(increasing_vol);

        SampleCovariance sample;
        auto sample_cov = sample.estimate_covariance(increasing_vol);

        // EWMA should show higher variance (recent data more volatile)
        REQUIRE(ewma_cov(0, 0) > sample_cov(0, 0));
    }
}

// ============================================================================
// LedoitWolfShrinkage Tests
// ============================================================================

TEST_CASE_METHOD(RiskModelTestFixture, "LedoitWolfShrinkage construction", "[RiskModel][LedoitWolfShrinkage]")
{

    SECTION("Default construction")
    {
        LedoitWolfShrinkage lw;
        REQUIRE(lw.get_target() == ShrinkageTarget::CONSTANT_CORRELATION);
        REQUIRE(lw.get_name() == "LedoitWolfShrinkage");
    }

    SECTION("Different targets")
    {
        LedoitWolfShrinkage lw_id(ShrinkageTarget::IDENTITY);
        REQUIRE(lw_id.get_target() == ShrinkageTarget::IDENTITY);

        LedoitWolfShrinkage lw_cv(ShrinkageTarget::CONSTANT_VARIANCE);
        REQUIRE(lw_cv.get_target() == ShrinkageTarget::CONSTANT_VARIANCE);

        LedoitWolfShrinkage lw_cc(ShrinkageTarget::CONSTANT_CORRELATION);
        REQUIRE(lw_cc.get_target() == ShrinkageTarget::CONSTANT_CORRELATION);
    }

    SECTION("Fixed shrinkage")
    {
        LedoitWolfShrinkage lw_fixed(ShrinkageTarget::IDENTITY, 0.5);
        auto result = lw_fixed.estimate_covariance(returns_small_);
        REQUIRE_THAT(lw_fixed.get_shrinkage_intensity(), WithinAbs(0.5, 1e-10));
    }

    SECTION("Invalid fixed shrinkage throws")
    {
        REQUIRE_THROWS_AS(
            LedoitWolfShrinkage(ShrinkageTarget::IDENTITY, 1.5),
            std::invalid_argument);
        REQUIRE_THROWS_AS(
            LedoitWolfShrinkage(ShrinkageTarget::IDENTITY, -0.5),
            std::invalid_argument);
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "LedoitWolfShrinkage basic estimation", "[RiskModel][LedoitWolfShrinkage]")
{

    SECTION("Shrinkage intensity in valid range")
    {
        LedoitWolfShrinkage lw;
        auto result = lw.estimate_covariance(returns_large_);

        double shrinkage = lw.get_shrinkage_intensity();
        REQUIRE(shrinkage >= 0.0);
        REQUIRE(shrinkage <= 1.0);
    }

    SECTION("Shrunk matrix properties")
    {
        LedoitWolfShrinkage lw;
        auto result = lw.estimate_covariance(returns_large_);

        // Check dimensions
        REQUIRE(result.rows() == 5);
        REQUIRE(result.cols() == 5);

        // Check symmetry
        REQUIRE(is_symmetric(result));

        // Check positive semi-definite
        REQUIRE(is_positive_semidefinite(result));
    }

    SECTION("Shrinkage improves conditioning for small samples")
    {
        // Use very small sample: 10 obs, 5 assets (T=2N)
        Eigen::MatrixXd small_sample = returns_large_.topRows(10);

        SampleCovariance sample;
        auto sample_cov = sample.estimate_covariance(small_sample);

        LedoitWolfShrinkage lw;
        auto shrunk_cov = lw.estimate_covariance(small_sample);

        // Shrunk matrix should have better condition number
        double cond_sample = compute_condition_number(sample_cov);
        double cond_shrunk = compute_condition_number(shrunk_cov);

        REQUIRE(cond_shrunk < cond_sample);

        // Shrinkage should be significant for small sample
        REQUIRE(lw.get_shrinkage_intensity() > 0.1);
    }

    SECTION("Different targets produce different results")
    {
        LedoitWolfShrinkage lw_id(ShrinkageTarget::IDENTITY);
        LedoitWolfShrinkage lw_cv(ShrinkageTarget::CONSTANT_VARIANCE);
        LedoitWolfShrinkage lw_cc(ShrinkageTarget::CONSTANT_CORRELATION);

        auto cov_id = lw_id.estimate_covariance(returns_small_);
        auto cov_cv = lw_cv.estimate_covariance(returns_small_);
        auto cov_cc = lw_cc.estimate_covariance(returns_small_);

        // Results should differ
        REQUIRE((cov_id - cov_cv).norm() > 1e-6);
        REQUIRE((cov_cv - cov_cc).norm() > 1e-6);
    }
}

// ============================================================================
// Correlation Estimation Tests
// ============================================================================

TEST_CASE_METHOD(RiskModelTestFixture, "Correlation matrix estimation", "[RiskModel]")
{
    SampleCovariance cov;
    auto corr = cov.estimate_correlation(returns_small_);

    SECTION("Diagonal elements are 1.0")
    {
        for (int i = 0; i < corr.rows(); ++i)
        {
            REQUIRE_THAT(corr(i, i), WithinAbs(1.0, 1e-10));
        }
    }

    SECTION("Symmetry")
    {
        REQUIRE(is_symmetric(corr));
    }

    SECTION("Off-diagonal in [-1, 1]")
    {
        for (int i = 0; i < corr.rows(); ++i)
        {
            for (int j = 0; j < corr.cols(); ++j)
            {
                if (i != j)
                {
                    REQUIRE(corr(i, j) >= -1.0);
                    REQUIRE(corr(i, j) <= 1.0);
                }
            }
        }
    }
}

// ============================================================================
// RiskModelFactory Tests
// ============================================================================

TEST_CASE("RiskModelFactory creation", "[RiskModel][Factory]")
{

    SECTION("Create sample covariance")
    {
        auto model = RiskModelFactory::create_sample_covariance(true);
        REQUIRE(model != nullptr);
        REQUIRE(model->get_name() == "SampleCovariance");
    }

    SECTION("Create EWMA")
    {
        auto model = RiskModelFactory::create_ewma_covariance(0.94);
        REQUIRE(model != nullptr);
        REQUIRE(model->get_name() == "EWMACovariance");
    }

    SECTION("Create Ledoit-Wolf")
    {
        auto model = RiskModelFactory::create_ledoit_wolf();
        REQUIRE(model != nullptr);
        REQUIRE(model->get_name() == "LedoitWolfShrinkage");
    }

    SECTION("Create from configuration")
    {
        RiskModelConfig config;
        config.type = "ewma";
        config.ewma_lambda = 0.97;

        auto model = RiskModelFactory::create(config);
        REQUIRE(model != nullptr);
        REQUIRE(model->get_name() == "EWMACovariance");
    }

    SECTION("Create from JSON")
    {
        nlohmann::json params = {
            {"ewma_lambda", 0.94}};

        auto model = RiskModelFactory::create("ewma", params);
        REQUIRE(model != nullptr);
    }

    SECTION("Unknown type throws")
    {
        RiskModelConfig config;
        config.type = "unknown_type";
        REQUIRE_THROWS_AS(RiskModelFactory::create(config), std::invalid_argument);
    }
}

TEST_CASE("RiskModelConfig JSON serialization", "[RiskModel][Factory]")
{

    SECTION("Round-trip JSON conversion")
    {
        RiskModelConfig original;
        original.type = "ewma";
        original.ewma_lambda = 0.94;
        original.bias_correction = true;

        auto json = original.to_json();
        auto restored = RiskModelConfig::from_json(json);

        REQUIRE(restored.type == original.type);
        REQUIRE_THAT(restored.ewma_lambda, WithinAbs(original.ewma_lambda, 1e-10));
        REQUIRE(restored.bias_correction == original.bias_correction);
    }

    SECTION("JSON with defaults")
    {
        nlohmann::json minimal = {{"type", "sample"}};
        auto config = RiskModelConfig::from_json(minimal);

        REQUIRE(config.type == "sample");
        REQUIRE(config.bias_correction == true);                  // Default
        REQUIRE_THAT(config.ewma_lambda, WithinAbs(0.94, 1e-10)); // Default
    }
}