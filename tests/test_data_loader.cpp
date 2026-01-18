/**
 * @file test_data_loader.cpp
 * @brief Unit tests for DataLoader and MarketData classes
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "data/data_loader.hpp"
#include "data/market_data.hpp"
#include <Eigen/Dense>

using namespace portfolio;
using Catch::Matchers::WithinAbs;

TEST_CASE("MarketData construction", "[MarketData]") {
    SECTION("Default constructor") {
        MarketData data(10, 5);
        REQUIRE(data.num_dates() == 10);
        REQUIRE(data.num_assets() == 5);
    }
    
    SECTION("Constructor with data") {
        Eigen::MatrixXd prices(3, 2);
        prices << 100.0, 150.0,
                  101.0, 151.5,
                  102.0, 152.0;
        
        std::vector<std::string> dates = {"2020-01-01", "2020-01-02", "2020-01-03"};
        std::vector<std::string> tickers = {"AAPL", "MSFT"};
        
        MarketData data(prices, dates, tickers);
        
        REQUIRE(data.num_dates() == 3);
        REQUIRE(data.num_assets() == 2);
        REQUIRE(data.get_dates() == dates);
        REQUIRE(data.get_tickers() == tickers);
    }
}

TEST_CASE("Return calculations", "[MarketData]") {
    Eigen::MatrixXd prices(4, 2);
    prices << 100.0, 200.0,
              110.0, 210.0,
              105.0, 220.0,
              115.0, 215.0;
    
    std::vector<std::string> dates = {"2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"};
    std::vector<std::string> tickers = {"AAPL", "MSFT"};
    
    MarketData data(prices, dates, tickers);
    
    SECTION("Simple returns") {
        auto returns = data.calculate_returns(ReturnType::SIMPLE);
        
        REQUIRE(returns.rows() == 3);
        REQUIRE(returns.cols() == 2);
        
        // First return for AAPL: (110-100)/100 = 0.10
        REQUIRE_THAT(returns(0, 0), WithinAbs(0.10, 0.0001));
        
        // First return for MSFT: (210-200)/200 = 0.05
        REQUIRE_THAT(returns(0, 1), WithinAbs(0.05, 0.0001));
    }
    
    SECTION("Log returns") {
        auto returns = data.calculate_returns(ReturnType::LOG);
        
        REQUIRE(returns.rows() == 3);
        REQUIRE(returns.cols() == 2);
        
        // First log return for AAPL: log(110/100) ≈ 0.0953
        REQUIRE_THAT(returns(0, 0), WithinAbs(0.0953, 0.0001));
    }
}

TEST_CASE("Data filtering", "[MarketData]") {
    Eigen::MatrixXd prices(5, 2);
    prices << 100.0, 200.0,
              110.0, 210.0,
              105.0, 220.0,
              115.0, 215.0,
              120.0, 225.0;
    
    std::vector<std::string> dates = {"2020-01-01", "2020-01-02", "2020-01-03", 
                                      "2020-01-04", "2020-01-05"};
    std::vector<std::string> tickers = {"AAPL", "MSFT"};
    
    MarketData data(prices, dates, tickers);
    
    SECTION("Filter by date range") {
        auto filtered = data.filter_by_date("2020-01-02", "2020-01-04");
        
        REQUIRE(filtered.num_dates() == 3);
        REQUIRE(filtered.get_dates()[0] == "2020-01-02");
        REQUIRE(filtered.get_dates()[2] == "2020-01-04");
    }
    
    SECTION("Select specific asset") {
        std::vector<std::string> selected = {"MSFT"};
        auto filtered = data.select_assets(selected);
        
        REQUIRE(filtered.num_assets() == 1);
        REQUIRE(filtered.get_tickers()[0] == "MSFT");
    }
}

TEST_CASE("Statistical methods", "[MarketData]") {
    Eigen::MatrixXd prices(6, 2);
    prices << 100.0, 200.0,
              110.0, 210.0,
              105.0, 220.0,
              115.0, 215.0,
              120.0, 225.0,
              118.0, 230.0;
    
    std::vector<std::string> dates = {"2020-01-01", "2020-01-02", "2020-01-03", 
                                      "2020-01-04", "2020-01-05", "2020-01-06"};
    std::vector<std::string> tickers = {"AAPL", "MSFT"};
    
    MarketData data(prices, dates, tickers);
    
    SECTION("Mean returns") {
        auto mean_ret = data.mean_returns(ReturnType::SIMPLE);
        
        REQUIRE(mean_ret.size() == 2);
        REQUIRE(mean_ret(0) > 0.0); // AAPL has positive average return
        REQUIRE(mean_ret(1) > 0.0); // MSFT has positive average return
    }
    
    SECTION("Volatilities") {
        auto vols = data.volatilities(ReturnType::SIMPLE);
        
        REQUIRE(vols.size() == 2);
        REQUIRE(vols(0) > 0.0);
        REQUIRE(vols(1) > 0.0);
    }
}

TEST_CASE("DataLoader synthetic data", "[DataLoader]") {
    std::vector<std::string> tickers = {"AAPL", "MSFT", "GOOGL"};
    
    auto data = DataLoader::generate_synthetic_data(tickers, 100, "2020-01-01");
    
    REQUIRE(data.num_dates() == 100);
    REQUIRE(data.num_assets() == 3);
    REQUIRE(data.is_valid());
}

TEST_CASE("JSON configuration loading", "[DataLoader]") {
    SECTION("Load configuration") {
        // This will throw if the config file doesn't exist
        REQUIRE_THROWS(DataLoader::load_config("nonexistent_config.json"));
    }
}
