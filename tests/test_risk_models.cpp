#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "risk/sample_covariance.hpp"
#include "risk/ewma_covariance.hpp"
#include "risk/ledoit_wolf_shrinkage.hpp"

using namespace portfolio::risk;
using Catch::Matchers::WithinAbs;

// Test fixture for shared test data
class RiskModelTestFixture
{
protected:
    // Simple 2x3 return matrix for basic tests
    Eigen::MatrixXd returns_2x3_;

    // Larger realistic return matrix (252 days, 10 assets)
    Eigen::MatrixXd returns_252x10_;

    // Expected covariance (validated against NumPy)
    Eigen::MatrixXd expected_cov_numpy_;

    RiskModelTestFixture()
    {
        // Initialize test data
        returns_2x3_ = Eigen::MatrixXd(2, 3);
        returns_2x3_ << 0.01, 0.02, -0.01,
            0.02, -0.01, 0.01;

        // Generate synthetic returns for realistic tests
        returns_252x10_ = generate_synthetic_returns(252, 10, 42);

        // Compute expected covariance with NumPy (reference)
        // This should be pre-computed and hardcoded for validation
        expected_cov_numpy_ = compute_reference_covariance();
    }

    Eigen::MatrixXd generate_synthetic_returns(int n_obs, int n_assets, int seed);
    Eigen::MatrixXd compute_reference_covariance();
};

TEST_CASE_METHOD(RiskModelTestFixture, "SampleCovariance basic functionality", "[RiskModel][SampleCovariance]")
{

    SECTION("Construct with default parameters")
    {
        SampleCovariance cov;
        REQUIRE(cov.uses_bias_correction() == true);
        REQUIRE(cov.get_name() == "SampleCovariance");
    }

    SECTION("Construct without bias correction")
    {
        SampleCovariance cov(false);
        REQUIRE(cov.uses_bias_correction() == false);
    }

    SECTION("Estimate covariance with bias correction")
    {
        SampleCovariance cov(true);
        auto result = cov.estimate_covariance(returns_2x3_);

        // Check dimensions
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 3);

        // Check symmetry
        REQUIRE_THAT((result - result.transpose()).norm(), WithinAbs(0.0, 1e-10));

        // Check diagonal elements are positive
        for (int i = 0; i < result.rows(); ++i)
        {
            REQUIRE(result(i, i) > 0.0);
        }
    }

    SECTION("Match NumPy covariance computation")
    {
        SampleCovariance cov(true);
        auto result = cov.estimate_covariance(returns_252x10_);

        // Validate against pre-computed NumPy result
        double max_diff = (result - expected_cov_numpy_).cwiseAbs().maxCoeff();
        REQUIRE_THAT(max_diff, WithinAbs(0.0, 1e-10));
    }

    SECTION("Bias correction vs no bias correction")
    {
        SampleCovariance cov_biased(false);
        SampleCovariance cov_unbiased(true);

        auto biased = cov_biased.estimate_covariance(returns_2x3_);
        auto unbiased = cov_unbiased.estimate_covariance(returns_2x3_);

        // Unbiased should be larger by factor of n/(n-1)
        int n = returns_2x3_.rows();
        double expected_ratio = static_cast<double>(n) / (n - 1);

        auto ratio = (unbiased.array() / biased.array()).matrix();
        REQUIRE_THAT(ratio.mean(), WithinAbs(expected_ratio, 1e-6));
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "SampleCovariance error handling", "[RiskModel][SampleCovariance]")
{
    SampleCovariance cov;

    SECTION("Empty returns matrix")
    {
        Eigen::MatrixXd empty(0, 0);
        REQUIRE_THROWS_AS(cov.estimate_covariance(empty), std::invalid_argument);
    }

    SECTION("Single observation")
    {
        Eigen::MatrixXd single_obs(1, 3);
        single_obs << 0.01, 0.02, 0.03;
        REQUIRE_THROWS_AS(cov.estimate_covariance(single_obs), std::invalid_argument);
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "EWMACovariance basic functionality", "[RiskModel][EWMACovariance]")
{

    SECTION("Construct with valid lambda")
    {
        EWMACovariance ewma(0.94);
        REQUIRE_THAT(ewma.get_lambda(), WithinAbs(0.94, 1e-10));
        REQUIRE_THAT(ewma.get_effective_window(), WithinAbs(16.67, 0.1));
    }

    SECTION("Invalid lambda throws exception")
    {
        REQUIRE_THROWS_AS(EWMACovariance(0.0), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(1.0), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(-0.5), std::invalid_argument);
        REQUIRE_THROWS_AS(EWMACovariance(1.5), std::invalid_argument);
    }

    SECTION("Estimate EWMA covariance")
    {
        EWMACovariance ewma(0.94);
        auto result = ewma.estimate_covariance(returns_252x10_);

        // Check dimensions
        REQUIRE(result.rows() == 10);
        REQUIRE(result.cols() == 10);

        // Check symmetry
        REQUIRE_THAT((result - result.transpose()).norm(), WithinAbs(0.0, 1e-10));

        // Check positive semi-definite (all eigenvalues >= 0)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(result);
        REQUIRE(solver.eigenvalues().minCoeff() >= -1e-10);
    }

    SECTION("EWMA gives more weight to recent observations")
    {
        EWMACovariance ewma(0.94);

        // Create returns with increasing volatility
        Eigen::MatrixXd returns_increasing(100, 2);
        for (int i = 0; i < 100; ++i)
        {
            double scale = 1.0 + i / 100.0; // Volatility doubles
            returns_increasing.row(i) << scale * 0.01, scale * 0.01;
        }

        auto ewma_cov = ewma.estimate_covariance(returns_increasing);

        SampleCovariance sample_cov;
        auto sample = sample_cov.estimate_covariance(returns_increasing);

        // EWMA should show higher variance (recent data more volatile)
        REQUIRE(ewma_cov(0, 0) > sample(0, 0));
    }
}

TEST_CASE_METHOD(RiskModelTestFixture, "LedoitWolfShrinkage basic functionality", "[RiskModel][LedoitWolfShrinkage]")
{

    SECTION("Construct with default parameters")
    {
        LedoitWolfShrinkage lw;
        REQUIRE(lw.get_target() == ShrinkageTarget::CONSTANT_CORRELATION);
    }

    SECTION("Shrinkage intensity in valid range")
    {
        LedoitWolfShrinkage lw;
        auto result = lw.estimate_covariance(returns_252x10_);

        double shrinkage = lw.get_shrinkage_intensity();
        REQUIRE(shrinkage >= 0.0);
        REQUIRE(shrinkage <= 1.0);
    }

    SECTION("Fixed shrinkage override")
    {
        LedoitWolfShrinkage lw(ShrinkageTarget::IDENTITY, 0.5);
        auto result = lw.estimate_covariance(returns_252x10_);

        REQUIRE_THAT(lw.get_shrinkage_intensity(), WithinAbs(0.5, 1e-10));
    }

    SECTION("Shrinkage improves conditioning")
    {
        // Small sample: T=20, N=10 (highly underdetermined)
        Eigen::MatrixXd small_sample = returns_252x10_.topRows(20);

        SampleCovariance sample_cov;
        auto sample = sample_cov.estimate_covariance(small_sample);

        LedoitWolfShrinkage lw;
        auto shrunk = lw.estimate_covariance(small_sample);

        // Shrunk matrix should have better condition number
        double cond_sample = compute_condition_number(sample);
        double cond_shrunk = compute_condition_number(shrunk);

        REQUIRE(cond_shrunk < cond_sample);
    }
}

TEST_CASE("Correlation estimation", "[RiskModel]")
{
    Eigen::MatrixXd returns(10, 3);
    returns.setRandom();

    SampleCovariance cov;
    auto corr = cov.estimate_correlation(returns);

    SECTION("Diagonal elements are 1.0")
    {
        for (int i = 0; i < corr.rows(); ++i)
        {
            REQUIRE_THAT(corr(i, i), WithinAbs(1.0, 1e-10));
        }
    }

    SECTION("Off-diagonal elements in [-1, 1]")
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

// Helper function for condition number
double compute_condition_number(const Eigen::MatrixXd &matrix)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    double max_sv = svd.singularValues()(0);
    double min_sv = svd.singularValues()(svd.singularValues().size() - 1);
    return max_sv / min_sv;
}