#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <iostream>
#include "backtest/backtest_engine.hpp"
#include "data/data_loader.hpp"

using namespace portfolio;
using namespace portfolio::backtest;

TEST_CASE("BacktestEngine construction", "[BacktestEngine]") {
    BacktestParams p;
    REQUIRE(p.initial_capital == 1000000.0);

    // from PortfolioConfig
    portfolio::PortfolioConfig cfg;
    cfg.backtest.initial_capital = 500000.0;
    cfg.backtest.lookback_window = 20;
    cfg.backtest.min_history = 5;
    cfg.backtest.rebalance_frequency = "monthly";
    BacktestParams p2 = BacktestParams::from_config(cfg);
    REQUIRE(p2.initial_capital == 500000.0);
    REQUIRE(p2.lookback_window == 20);
    REQUIRE(p2.min_history == 5);
}

TEST_CASE("BacktestEngine basic run", "[BacktestEngine]") {
    std::vector<std::string> tickers = {"A","B","C"};
    size_t days = 60;
    auto data = DataLoader::generate_synthetic_data(tickers, days, "2020-01-01");

    BacktestParams params;
    params.initial_capital = 100000.0;
    params.min_history = 2; // need at least 2 return observations
    params.lookback_window = 10;

    BacktestEngine engine(params);
    auto result = engine.run(data);

    REQUIRE(result.success == true);
    std::cerr << "TEST DEBUG: snapshots.size=" << result.snapshots.size() << " nav_series.size=" << result.nav_series.size() << "\n";
    REQUIRE(result.snapshots.size() == days - 1);
    REQUIRE(result.nav_series.size() == days - 1);
    REQUIRE(result.nav_series.front() == params.initial_capital);
    REQUIRE(result.dates.front() < result.dates.back());
}

TEST_CASE("BacktestEngine edge cases", "[BacktestEngine]") {
    std::vector<std::string> tickers = {"A"};
    size_t days = 10;
    auto data = DataLoader::generate_synthetic_data(tickers, days, "2020-01-01");

    BacktestParams params;
    params.min_history = 20; // larger than available
    BacktestEngine engine(params);
    REQUIRE_THROWS_AS(engine.run(data), std::invalid_argument);

    // Single asset runs with small history
    params.min_history = 2;
    BacktestEngine engine2(params);
    auto res = engine2.run(data);
    REQUIRE(res.success);
}

TEST_CASE("BacktestEngine CSV export", "[BacktestEngine]") {
    std::vector<std::string> tickers = {"A","B"};
    size_t days = 30;
    auto data = DataLoader::generate_synthetic_data(tickers, days, "2020-01-01");

    BacktestParams params;
    params.min_history = 2;
    BacktestEngine engine(params);
    auto res = engine.run(data);

    std::string tmp = "/tmp/backtest_nav.csv";
    res.export_nav_to_csv(tmp);

    std::ifstream in(tmp);
    REQUIRE(in.good());
    std::string line;
    int count = 0;
    while (std::getline(in, line)) ++count;
    // header + rows == nav_series.size()
    REQUIRE(count == static_cast<int>(res.nav_series.size()) + 1);
}
