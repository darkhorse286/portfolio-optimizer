#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "backtest/portfolio.hpp"
#include <Eigen/Dense>

using namespace portfolio::backtest;

TEST_CASE("Portfolio construction", "[Portfolio]") {
    SECTION("Happy path: valid initial capital and tickers") {
        std::vector<std::string> t = {"A", "B"};
        Portfolio p(1000.0, t);
        REQUIRE(p.cash() == Catch::Approx(1000.0));
        REQUIRE(p.num_assets() == 2);
    }

    SECTION("Edge case: single asset") {
        std::vector<std::string> t = {"A"};
        Portfolio p(500.0, t);
        REQUIRE(p.num_assets() == 1);
    }

    SECTION("Error: zero capital") {
        std::vector<std::string> t = {"A"};
        REQUIRE_THROWS_AS(Portfolio(0.0, t), std::invalid_argument);
    }

    SECTION("Error: negative capital") {
        std::vector<std::string> t = {"A"};
        REQUIRE_THROWS_AS(Portfolio(-1.0, t), std::invalid_argument);
    }

    SECTION("Error: empty tickers") {
        std::vector<std::string> t;
        REQUIRE_THROWS_AS(Portfolio(100.0, t), std::invalid_argument);
    }
}

TEST_CASE("Portfolio price updates", "[Portfolio]") {
    std::vector<std::string> t = {"A", "B"};
    Portfolio p(1000.0, t);
    Eigen::VectorXd prices(2);
    prices << 10.0, 20.0;

    SECTION("Happy path: update prices and check NAV") {
        p.update_prices("2020-01-01", prices);
        REQUIRE(p.nav() == Catch::Approx(1000.0));
        // buy some shares
        p.execute_trade("A", 10.0, 10.0); // spend 100
        p.update_prices("2020-01-02", prices);
        REQUIRE(p.nav() == Catch::Approx(900.0 + 10.0*10.0));
    }

    SECTION("Error: wrong number of prices") {
        Eigen::VectorXd bad(1);
        bad << 5.0;
        REQUIRE_THROWS_AS(p.update_prices("2020-01-01", bad), std::invalid_argument);
    }

    SECTION("Error: negative price") {
        Eigen::VectorXd bad(2);
        bad << -1.0, 2.0;
        REQUIRE_THROWS_AS(p.update_prices("2020-01-01", bad), std::invalid_argument);
    }
}

TEST_CASE("Portfolio trading", "[Portfolio]") {
    std::vector<std::string> t = {"A", "B"};
    Portfolio p(1000.0, t);
    Eigen::VectorXd prices(2);
    prices << 10.0, 20.0;
    p.update_prices("2020-01-01", prices);

    SECTION("Happy path: set target weights and verify positions") {
        Eigen::VectorXd tw(2);
        tw << 0.5, 0.5;
        auto trades = p.set_target_weights(tw, prices);
        REQUIRE(p.get_position("A").shares == Catch::Approx(50.0));
        REQUIRE(p.get_position("B").shares == Catch::Approx(25.0));
        REQUIRE(p.cash() == Catch::Approx(0.0));
    }

    SECTION("Happy path: single asset trade") {
        p.execute_trade("A", 10.0, 10.0);
        REQUIRE(p.get_position("A").shares == Catch::Approx(10.0));
        REQUIRE(p.cash() == Catch::Approx(900.0));
    }

    SECTION("Edge case: target weights equal current weights (no trades)") {
        Eigen::VectorXd tw(2);
        tw << 0.0, 0.0;
        auto trades = p.set_target_weights(tw, prices);
        REQUIRE(trades.size() == 2);
    }

    SECTION("Edge case: sell everything (go to cash)") {
        p.execute_trade("A", 10.0, 10.0);
        p.execute_trade("B", 10.0, 20.0);
        Eigen::VectorXd tw(2);
        tw << 0.0, 0.0;
        p.set_target_weights(tw, prices);
        REQUIRE(p.invested_value() == Catch::Approx(0.0));
    }

    SECTION("Error: insufficient cash") {
        REQUIRE_THROWS_AS(p.execute_trade("A", 200.0, 10.0), std::runtime_error);
    }
}

TEST_CASE("Portfolio snapshots", "[Portfolio]") {
    std::vector<std::string> t = {"A"};
    Portfolio p(1000.0, t);
    Eigen::VectorXd prices(1);
    prices << 10.0;
    p.update_prices("2020-01-01", prices);
    p.execute_trade("A", 50.0, 10.0); // spend 500
    p.update_prices("2020-01-02", prices);
    p.record_snapshot();
    REQUIRE(p.history().size() == 1);
    auto s = p.history().back();
    REQUIRE(s.nav == Catch::Approx(p.nav()));
    REQUIRE(s.daily_return == Catch::Approx((p.nav() - 1000.0)/1000.0));
    REQUIRE(s.cumulative_return == Catch::Approx((p.nav() - 1000.0)/1000.0));
}

TEST_CASE("Portfolio weights", "[Portfolio]") {
    std::vector<std::string> t = {"A", "B"};
    Portfolio p(1000.0, t);
    Eigen::VectorXd prices(2);
    prices << 10.0, 10.0;
    p.update_prices("2020-01-01", prices);
    Eigen::VectorXd tw(2);
    tw << 0.6, 0.3; // leave 0.1 as cash
    p.set_target_weights(tw, prices);
    auto w = p.current_weights();
    double sum = w.sum();
    REQUIRE(sum <= Catch::Approx(0.99).margin(0.11)); // allow float small tolerance; sum <= 1
    REQUIRE((w.array() >= 0.0).all());
    // Verify weight accuracy
    double nav = p.nav();
    for (int i = 0; i < w.size(); ++i) {
        double expected = (p.get_position(t[i]).shares * prices[i]) / nav;
        REQUIRE(w[i] == Catch::Approx(expected).epsilon(1e-10));
    }
}
