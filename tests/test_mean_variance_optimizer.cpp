/**
 * @file test_mean_variance_optimizer.cpp
 * @brief Unit tests and convergence diagnostics for MeanVarianceOptimizer
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>

#include "optimizer/mean_variance_optimizer.hpp"
#include "risk/sample_covariance.hpp"
#include "data/data_loader.hpp"
#include "data/market_data.hpp"


using namespace portfolio;
using namespace portfolio::optimizer;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Test Fixture
// ============================================================================

class OptimizerTestFixture
{
protected:
    // Simple test data
    Eigen::VectorXd returns_2asset_;
    Eigen::MatrixXd cov_2asset_;
    
    // Medium complexity
    Eigen::VectorXd returns_10asset_;
    Eigen::MatrixXd cov_10asset_;
    
    OptimizerTestFixture()
    {
        setup_2asset_data();
        setup_10asset_data();
    }
    
    void setup_2asset_data()
    {
        // Two asset case with known analytical solution
        returns_2asset_ = Eigen::VectorXd(2);
        returns_2asset_ << 0.10, 0.08;  // 10% and 8% expected returns
        
        cov_2asset_ = Eigen::MatrixXd(2, 2);
        cov_2asset_ << 0.04, 0.01,   // Asset 1: 20% vol
                       0.01, 0.02;   // Asset 2: 14.1% vol, correlation ~0.35
    }
    
    void setup_10asset_data()
    {
        // 10 assets with varying characteristics
        returns_10asset_ = Eigen::VectorXd(10);
        returns_10asset_ << 0.10, 0.08, 0.12, 0.09, 0.11,
                           0.07, 0.13, 0.08, 0.10, 0.09;
        
        // Create reasonable covariance matrix
        // Start with diagonal (individual variances)
        cov_10asset_ = Eigen::MatrixXd::Identity(10, 10);
        cov_10asset_ *= 0.04;  // 20% volatility for all
        
        // Add correlation structure
        for (int i = 0; i < 10; ++i) {
            for (int j = i + 1; j < 10; ++j) {
                double correlation = 0.3;  // Moderate correlation
                cov_10asset_(i, j) = correlation * 0.04;
                cov_10asset_(j, i) = correlation * 0.04;
            }
        }
    }
};

// ============================================================================
// PHASE 1: Basic Unit Tests
// ============================================================================

TEST_CASE_METHOD(OptimizerTestFixture, "MeanVarianceOptimizer construction", 
                 "[MeanVarianceOptimizer][Basic]")
{
    SECTION("Default construction") {
        MeanVarianceOptimizer opt;
        REQUIRE(opt.get_name() == "MeanVarianceOptimizer");
        REQUIRE(opt.get_objective() == ObjectiveType::MAX_SHARPE);
        REQUIRE_THAT(opt.get_risk_aversion(), WithinAbs(1.0, 1e-10));
    }
    
    SECTION("Construction with parameters") {
        MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE, 0.02);
        REQUIRE(opt.get_objective() == ObjectiveType::MIN_VARIANCE);
    }
    
    SECTION("Invalid risk-free rate throws") {
        REQUIRE_THROWS_AS(
            MeanVarianceOptimizer(ObjectiveType::MAX_SHARPE, -0.01),
            std::invalid_argument
        );
    }
}

TEST_CASE_METHOD(OptimizerTestFixture, "OptimizationConstraints validation",
                 "[MeanVarianceOptimizer][Constraints]")
{
    SECTION("Valid constraints") {
        OptimizationConstraints c;
        c.min_weight = 0.0;
        c.max_weight = 0.5;
        c.long_only = true;
        c.sum_to_one = true;
        
        REQUIRE_NOTHROW(c.validate());
    }
    
    SECTION("min_weight > max_weight throws") {
        OptimizationConstraints c;
        c.min_weight = 0.6;
        c.max_weight = 0.4;
        
        REQUIRE_THROWS_AS(c.validate(), std::invalid_argument);
    }
    
    SECTION("Negative min_weight with long_only throws") {
        OptimizationConstraints c;
        c.min_weight = -0.1;
        c.long_only = true;
        
        REQUIRE_THROWS_AS(c.validate(), std::invalid_argument);
    }
}

TEST_CASE_METHOD(OptimizerTestFixture, "Minimum variance - simple case",
                 "[MeanVarianceOptimizer][MinVariance]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    SECTION("Two-asset portfolio") {
        auto result = opt.optimize(returns_2asset_, cov_2asset_, constraints);
        
        // Basic checks
        REQUIRE(result.success);
        INFO("Optimization message: " << result.message);
        REQUIRE(result.iterations < 5000);  // Should converge well before max
        
        // Weights sum to 1
        REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));
        
        // Weights in bounds
        for (int i = 0; i < result.weights.size(); ++i) {
            REQUIRE(result.weights(i) >= -1e-6);
            REQUIRE(result.weights(i) <= 1.0 + 1e-6);
        }
        
        // For min variance, should favor lower variance asset (asset 2)
        REQUIRE(result.weights(1) > result.weights(0));
        
        // Volatility should be reasonable
        REQUIRE(result.volatility > 0.0);
        REQUIRE(result.volatility < 0.25);  // Less than either individual asset
    }
    
    SECTION("10-asset portfolio") {
        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
        
        REQUIRE(result.success);
        INFO("Optimization message: " << result.message);
        REQUIRE(result.iterations < 5000);
        
        // Constraints satisfied
        REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));
        
        for (int i = 0; i < result.weights.size(); ++i) {
            REQUIRE(result.weights(i) >= -1e-6);
            REQUIRE(result.weights(i) <= 1.0 + 1e-6);
        }
    }
}

TEST_CASE_METHOD(OptimizerTestFixture, "Box constraints enforcement",
                 "[MeanVarianceOptimizer][Constraints]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.05;  // At least 5% per asset
    constraints.max_weight = 0.30;  // At most 30% per asset
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
    
    REQUIRE(result.success);
    INFO("Optimization message: " << result.message);
    
    // Check every weight is in bounds
    for (int i = 0; i < result.weights.size(); ++i) {
        INFO("Weight " << i << ": " << result.weights(i));
        REQUIRE(result.weights(i) >= constraints.min_weight - 1e-4);
        REQUIRE(result.weights(i) <= constraints.max_weight + 1e-4);
    }
    
    // Sum to one
    REQUIRE_THAT(result.weights.sum(), WithinAbs(1.0, 1e-3));
}

// ============================================================================
// PHASE 2: Convergence Diagnostics
// ============================================================================

TEST_CASE_METHOD(OptimizerTestFixture, "Convergence diagnostics - constraint tightness",
                 "[MeanVarianceOptimizer][Convergence][Diagnostics]")
{
    struct TestCase {
        std::string name;
        double min_weight;
        double max_weight;
    };
    
    std::vector<TestCase> cases = {
        {"Loose constraints",      0.0,   1.0},
        {"Moderate constraints",   0.02,  0.30},
        {"Tight constraints",      0.05,  0.25},
        {"Very tight constraints", 0.08,  0.15}
    };
    
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    std::cout << "\n=== Convergence Diagnostics: Constraint Tightness ===\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::setw(25) << "Test Case" 
              << std::setw(12) << "Success" 
              << std::setw(12) << "Iterations"
              << std::setw(15) << "Volatility"
              << std::setw(20) << "Message\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& test : cases) {
        OptimizationConstraints constraints;
        constraints.min_weight = test.min_weight;
        constraints.max_weight = test.max_weight;
        constraints.long_only = true;
        constraints.sum_to_one = true;
        
        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
        
        std::cout << std::setw(25) << test.name
                  << std::setw(12) << (result.success ? "YES" : "NO")
                  << std::setw(12) << result.iterations
                  << std::setw(14) << std::fixed << std::setprecision(4)
                  << result.volatility * 100 << "%"
                  << "  " << result.message.substr(0, 30) << "\n";
        
        // All should converge
        REQUIRE(result.success);
        
        // Should not hit max iterations
        if (result.iterations >= 4900) {
            WARN("Near max iterations for: " << test.name);
        }
    }
    std::cout << std::string(80, '-') << "\n\n";
}

TEST_CASE_METHOD(OptimizerTestFixture, "Convergence diagnostics - problem size",
                 "[MeanVarianceOptimizer][Convergence][Diagnostics]")
{
    std::vector<int> sizes = {2, 5, 10, 20};
    
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    std::cout << "\n=== Convergence Diagnostics: Problem Size ===\n";
    std::cout << std::string(70, '-') << "\n";
    std::cout << std::setw(15) << "Num Assets"
              << std::setw(12) << "Success"
              << std::setw(15) << "Iterations"
              << std::setw(15) << "Time (ms)"
              << std::setw(15) << "Message\n";
    std::cout << std::string(70, '-') << "\n";
    
    for (int n : sizes) {
        // Generate problem of size n
        Eigen::VectorXd returns = (Eigen::VectorXd::Random(n) * 0.05).array() + 0.08;
        
        // Create covariance matrix
        Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(n, n) * 0.04;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                cov(i, j) = 0.01;
                cov(j, i) = 0.01;
            }
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        auto result = opt.optimize(returns, cov, constraints);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::cout << std::setw(15) << n
                  << std::setw(12) << (result.success ? "YES" : "NO")
                  << std::setw(15) << result.iterations
                  << std::setw(15) << duration
                  << "  " << result.message.substr(0, 20) << "\n";
        
        REQUIRE(result.success);
    }
    std::cout << std::string(70, '-') << "\n\n";
}

TEST_CASE_METHOD(OptimizerTestFixture, "Max Sharpe convergence test",
                 "[MeanVarianceOptimizer][MaxSharpe][Convergence]")
{
    MeanVarianceOptimizer opt(ObjectiveType::MAX_SHARPE, 0.02);
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 0.30;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    std::cout << "\n=== Max Sharpe Optimization ===\n";
    
    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
    
    std::cout << "Success: " << (result.success ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Sharpe Ratio: " << result.sharpe_ratio << "\n";
    std::cout << "Message: " << result.message << "\n";
    std::cout << "================================\n\n";
    
    REQUIRE(result.success);
    
    // Max Sharpe uses grid search, so iterations are cumulative
    // But each sub-optimization should converge
    REQUIRE(result.sharpe_ratio > 0.0);
    REQUIRE(result.sharpe_ratio < 10.0);  // Reasonable bound
}

TEST_CASE_METHOD(OptimizerTestFixture, "Risk aversion optimization",
                 "[MeanVarianceOptimizer][RiskAversion]")
{
    std::vector<double> lambdas = {0.5, 1.0, 2.0, 5.0};
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 0.30;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    std::cout << "\n=== Risk Aversion Tests ===\n";
    std::cout << std::string(70, '-') << "\n";
    std::cout << std::setw(12) << "Lambda"
              << std::setw(12) << "Success"
              << std::setw(12) << "Iterations"
              << std::setw(15) << "Return"
              << std::setw(15) << "Volatility\n";
    std::cout << std::string(70, '-') << "\n";
    
    for (double lambda : lambdas) {
        MeanVarianceOptimizer opt(ObjectiveType::RISK_AVERSION, 0.02);
        opt.set_risk_aversion(lambda);
        
        auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
        
        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << lambda
                  << std::setw(12) << (result.success ? "YES" : "NO")
                  << std::setw(12) << result.iterations
                  << std::setw(14) << std::setprecision(4) 
                  << result.expected_return * 100 << "%"
                  << std::setw(14) << result.volatility * 100 << "%\n";
        
        REQUIRE(result.success);
    }
    std::cout << std::string(70, '-') << "\n\n";
}

// ============================================================================
// PHASE 3: Real Data Test
// ============================================================================

TEST_CASE("Real market data convergence", "[MeanVarianceOptimizer][Integration]")
{
    // Load configuration and data
    auto config = DataLoader::load_config("../data/config/portfolio_config.json");
    
    // Load a subset of data for faster testing
    std::vector<std::string> test_tickers = {"AAPL", "MSFT", "JPM", "JNJ", "XOM"};
    auto data = DataLoader::load_csv(config.data.data_file, test_tickers);
    
    // Use last 252 days (1 year)
    if (data.num_dates() > 252) {
        std::string end_date = data.get_dates().back();
        std::string start_date = data.get_dates()[data.num_dates() - 252];
        data = data.filter_by_date(start_date, end_date);
    }
    
    auto returns_matrix = data.calculate_returns(ReturnType::SIMPLE);
    auto mean_returns = data.mean_returns(ReturnType::SIMPLE);
    
    // Use sample covariance
    risk::SampleCovariance cov_estimator(true);
    auto covariance = cov_estimator.estimate_covariance(returns_matrix);
    
    std::cout << "\n=== Real Data Optimization ===\n";
    std::cout << "Assets: " << data.num_assets() << "\n";
    std::cout << "Periods: " << data.num_dates() << "\n\n";
    
    // Test min variance
    {
        MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
        
        OptimizationConstraints constraints;
        constraints.min_weight = 0.05;
        constraints.max_weight = 0.30;
        constraints.long_only = true;
        constraints.sum_to_one = true;
        
        auto result = opt.optimize(mean_returns, covariance, constraints);
        
        std::cout << "Min Variance Portfolio:\n";
        std::cout << "  Success: " << (result.success ? "YES" : "NO") << "\n";
        std::cout << "  Iterations: " << result.iterations << "\n";
        std::cout << "  Volatility: " << result.volatility * std::sqrt(252) * 100 << "% (ann.)\n";
        std::cout << "  Message: " << result.message << "\n\n";
        
        REQUIRE(result.success);
        REQUIRE(result.iterations < 5000);
    }
    
    // Test max Sharpe
    {
        MeanVarianceOptimizer opt(ObjectiveType::MAX_SHARPE, 0.02);
        
        OptimizationConstraints constraints;
        constraints.min_weight = 0.05;
        constraints.max_weight = 0.30;
        constraints.long_only = true;
        constraints.sum_to_one = true;
        
        auto result = opt.optimize(mean_returns, covariance, constraints);
        
        std::cout << "Max Sharpe Portfolio:\n";
        std::cout << "  Success: " << (result.success ? "YES" : "NO") << "\n";
        std::cout << "  Iterations: " << result.iterations << "\n";
        std::cout << "  Sharpe: " << result.sharpe_ratio << "\n";
        std::cout << "  Return: " << result.expected_return * 252 * 100 << "% (ann.)\n";
        std::cout << "  Volatility: " << result.volatility * std::sqrt(252) * 100 << "% (ann.)\n";
        std::cout << "  Message: " << result.message << "\n\n";
        
        REQUIRE(result.success);
    }
    
    std::cout << "================================\n\n";
}

// ============================================================================
// PHASE 4: Stress Tests - Identify Failure Modes
// ============================================================================

TEST_CASE_METHOD(OptimizerTestFixture, "Ill-conditioned covariance matrix",
                 "[MeanVarianceOptimizer][StressTest]")
{
    // Create nearly singular covariance matrix
    Eigen::MatrixXd ill_cov = Eigen::MatrixXd::Identity(5, 5) * 0.04;
    
    // Make highly correlated
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 5; ++j) {
            ill_cov(i, j) = 0.039;  // correlation ~0.975
            ill_cov(j, i) = 0.039;
        }
    }
    
    Eigen::VectorXd returns = Eigen::VectorXd(5);
    returns << 0.10, 0.10, 0.10, 0.10, 0.10;
    
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.0;
    constraints.max_weight = 1.0;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    auto result = opt.optimize(returns, ill_cov, constraints);
    
    std::cout << "\n=== Ill-Conditioned Matrix Test ===\n";
    std::cout << "Success: " << (result.success ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Message: " << result.message << "\n";
    std::cout << "====================================\n\n";
    
    // This might fail or take many iterations
    if (!result.success) {
        WARN("Failed with ill-conditioned matrix - expected behavior");
    }
}

TEST_CASE_METHOD(OptimizerTestFixture, "Conflicting tight constraints",
                 "[MeanVarianceOptimizer][StressTest]")
{
    // 10 assets, each must be between 9-11%
    // This is barely feasible (sum to 100%)
    
    OptimizationConstraints constraints;
    constraints.min_weight = 0.09;
    constraints.max_weight = 0.11;
    constraints.long_only = true;
    constraints.sum_to_one = true;
    
    MeanVarianceOptimizer opt(ObjectiveType::MIN_VARIANCE);
    
    auto result = opt.optimize(returns_10asset_, cov_10asset_, constraints);
    
    std::cout << "\n=== Conflicting Constraints Test ===\n";
    std::cout << "Success: " << (result.success ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Weights sum: " << result.weights.sum() << "\n";
    std::cout << "Message: " << result.message << "\n";
    std::cout << "====================================\n\n";
    
    // Might fail or struggle
    if (!result.success) {
        WARN("Failed with tight constraints - may need better initialization");
    }
}