#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "quantum/benchmark_runner.hpp"
#include "quantum/simulated_annealing_solver.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
#include "data/data_loader.hpp"
#include "backtest/backtest_engine.hpp"

#include <memory>
#include <fstream>

std::vector<std::string> create_test_tickers(int num_assets) {
    std::vector<std::string> tickers;
    for (int i = 0; i < num_assets; ++i) {
        tickers.push_back("ASSET_" + std::to_string(i));
    }
    return tickers;
}

TEST_CASE("BenchmarkRunner registers solvers", "[BenchmarkRunner]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::quantum::BenchmarkRunner runner(data, params);

    auto classical = std::make_shared<portfolio::optimizer::MeanVarianceOptimizer>(
        portfolio::optimizer::ObjectiveType::MAX_SHARPE, 0.02);
    runner.add_classical_solver("markowitz", classical);

    portfolio::quantum::SimulatedAnnealingConfig sa_config;
    sa_config.max_iterations = 100; // Reduce for testing
    auto quantum = std::make_shared<portfolio::quantum::SimulatedAnnealingSolver>(sa_config);
    runner.add_quantum_solver("sa", "quantum_inspired", quantum);

    auto result = runner.run(1);
    REQUIRE(result.results.size() == 2);
}

TEST_CASE("BenchmarkRunner run with SA only produces valid result", "[BenchmarkRunner]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::quantum::BenchmarkRunner runner(data, params);

    portfolio::quantum::SimulatedAnnealingConfig sa_config;
    sa_config.max_iterations = 100;
    auto sa = std::make_shared<portfolio::quantum::SimulatedAnnealingSolver>(sa_config);
    runner.add_quantum_solver("sa", "quantum_inspired", sa);

    auto result = runner.run(1);
    REQUIRE(result.results.size() == 1);
    REQUIRE(result.results[0].successful_runs == 1);
    REQUIRE(std::isfinite(result.results[0].mean_sharpe));
}

TEST_CASE("BenchmarkRunner run with classical produces valid result", "[BenchmarkRunner]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::quantum::BenchmarkRunner runner(data, params);

    auto mv = std::make_shared<portfolio::optimizer::MeanVarianceOptimizer>(
        portfolio::optimizer::ObjectiveType::MAX_SHARPE, 0.02);
    runner.add_classical_solver("markowitz", mv);

    auto result = runner.run(1);
    REQUIRE(result.results.size() == 1);
    REQUIRE(result.results[0].solution_quality_vs_classical == Catch::Approx(1.0));
}

TEST_CASE("BenchmarkRunner comparison_results.json is valid JSON", "[BenchmarkRunner]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::quantum::BenchmarkRunner runner(data, params);

    auto mv = std::make_shared<portfolio::optimizer::MeanVarianceOptimizer>(
        portfolio::optimizer::ObjectiveType::MAX_SHARPE, 0.02);
    runner.add_classical_solver("markowitz", mv);

    auto result = runner.run(1);
    result.export_comparison_json("/tmp/test_comparison.json");

    std::ifstream in("/tmp/test_comparison.json");
    REQUIRE(in.good());
    nlohmann::json j = nlohmann::json::parse(in);
    REQUIRE(j.contains("schema_version"));
    REQUIRE(j["schema_version"] == "1.0");
    REQUIRE(j.contains("runs"));
    REQUIRE(j["runs"].is_array());
}

TEST_CASE("BenchmarkRunner print_summary does not throw", "[BenchmarkRunner]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::quantum::BenchmarkRunner runner(data, params);

    auto mv = std::make_shared<portfolio::optimizer::MeanVarianceOptimizer>(
        portfolio::optimizer::ObjectiveType::MAX_SHARPE, 0.02);
    runner.add_classical_solver("markowitz", mv);

    auto result = runner.run(1);
    REQUIRE_NOTHROW(result.print_summary());
}

TEST_CASE("BacktestEngine run with QuantumOptimizer overload", "[BacktestEngine][quantum]")
{
    auto tickers = create_test_tickers(10);
    auto data = portfolio::DataLoader::generate_synthetic_data(tickers, 252);
    portfolio::backtest::BacktestParams params;
    params.initial_capital = 1000000.0;
    params.min_history = 60;

    portfolio::backtest::BacktestEngine engine(params);
    portfolio::quantum::SimulatedAnnealingConfig sa_config;
    sa_config.max_iterations = 100;
    portfolio::quantum::SimulatedAnnealingSolver sa(sa_config);

    auto result = engine.run(data, sa);
    REQUIRE(result.success == true);
    REQUIRE(!result.nav_series.empty());
}