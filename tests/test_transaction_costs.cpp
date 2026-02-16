// Unit tests for TransactionCostModel

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "backtest/transaction_cost_model.hpp"

using namespace portfolio::backtest;
using Catch::Matchers::WithinAbs;

TEST_CASE("TransactionCostModel construction", "[TransactionCostModel]") {
    SECTION("Default config") {
        TransactionCostModel m;
        REQUIRE(m.get_name() == "TransactionCostModel");
        REQUIRE(m.config().commission_rate == 0.0);
        REQUIRE(m.config().market_impact_model == "none");
    }

    SECTION("Custom config") {
        TransactionCostConfig cfg;
        cfg.commission_rate = 0.001;
        cfg.slippage_bps = 5.0;
        cfg.market_impact_model = "linear";
        cfg.market_impact_coeff = 0.002;
        TransactionCostModel m(cfg);
        REQUIRE(m.config().commission_rate == 0.001);
    }

    SECTION("Error: negative commission rate") {
        TransactionCostConfig cfg = TransactionCostConfig::default_config();
        cfg.commission_rate = -0.1;
        REQUIRE_THROWS_AS(TransactionCostModel(cfg), std::invalid_argument);
    }

    SECTION("Error: negative slippage") {
        TransactionCostConfig cfg = TransactionCostConfig::default_config();
        cfg.slippage_bps = -1.0;
        REQUIRE_THROWS_AS(TransactionCostModel(cfg), std::invalid_argument);
    }
}

TEST_CASE("Commission calculation", "[TransactionCostModel]") {
    TransactionCostConfig cfg = TransactionCostConfig::default_config();
    cfg.commission_rate = 0.001; // 10 bps
    TransactionCostModel m(cfg);

    TradeOrder buy{"ABC", 100.0, 100.0}; // $10,000
    TradeCost c = m.calculate_cost(buy);
    REQUIRE_THAT(c.commission, WithinAbs(10000.0 * 0.001, 1e-12));

    TradeOrder zero{"ABC", 0.0, 100.0};
    REQUIRE_THAT(m.calculate_cost(zero).commission, WithinAbs(0.0, 1e-12));

    TradeOrder sell{"ABC", -100.0, 100.0};
    REQUIRE_THAT(m.calculate_cost(sell).commission, WithinAbs(c.commission, 1e-12));
}

TEST_CASE("Slippage calculation", "[TransactionCostModel]") {
    TransactionCostConfig cfg = TransactionCostConfig::default_config();
    cfg.slippage_bps = 5.0; // 5 bps
    TransactionCostModel m(cfg);

    TradeOrder o{"XYZ", 100.0, 100.0}; // $10,000
    TradeCost tc = m.calculate_cost(o);
    REQUIRE_THAT(tc.slippage, WithinAbs(10000.0 * (5.0 / 10000.0), 1e-12));

    cfg.slippage_bps = 0.0;
    TransactionCostModel m2(cfg);
    REQUIRE_THAT(m2.calculate_cost(o).slippage, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Market impact models", "[TransactionCostModel]") {
    TradeOrder o{"MKT", 100.0, 100.0}; // $10,000
    double notional = o.notional();

    SECTION("Linear impact: known value") {
        TransactionCostConfig cfg = TransactionCostConfig::default_config();
        cfg.market_impact_model = "linear";
        cfg.market_impact_coeff = 0.002;
        TransactionCostModel m(cfg);
        double expected = notional * cfg.market_impact_coeff;
        REQUIRE_THAT(m.calculate_cost(o).market_impact, WithinAbs(expected, 1e-12));
    }

    SECTION("Sqrt impact: known value") {
        TransactionCostConfig cfg = TransactionCostConfig::default_config();
        cfg.market_impact_model = "sqrt";
        cfg.market_impact_coeff = 1e-4;
        TransactionCostModel m(cfg);
        double expected = cfg.market_impact_coeff * std::sqrt(notional * o.price);
        REQUIRE_THAT(m.calculate_cost(o).market_impact, WithinAbs(expected, 1e-10));
    }

    SECTION("None: zero cost") {
        TransactionCostModel m;
        REQUIRE_THAT(m.calculate_cost(o).market_impact, WithinAbs(0.0, 1e-12));
    }

    SECTION("Sqrt increases sub-linearly with size") {
        TransactionCostConfig cfg = TransactionCostConfig::default_config();
        cfg.market_impact_model = "sqrt";
        cfg.market_impact_coeff = 1e-6;
        TransactionCostModel m(cfg);
        TradeOrder small{"A", 100.0, 100.0};
        TradeOrder big{"A", 400.0, 100.0};
        double c_small = m.calculate_cost(small).market_impact;
        double c_big = m.calculate_cost(big).market_impact;
        REQUIRE(c_big / c_small < 4.0); // sub-linear: quadrupling notional should give <4x cost
    }
}

TEST_CASE("Total cost calculation", "[TransactionCostModel]") {
    TransactionCostConfig cfg;
    cfg.commission_rate = 0.001;
    cfg.slippage_bps = 5.0;
    cfg.market_impact_model = "linear";
    cfg.market_impact_coeff = 0.002;
    TransactionCostModel m(cfg);

    std::vector<TradeOrder> orders{
        {"A", 100.0, 100.0},
        {"B", -50.0, 200.0}
    };
    double sum_components = 0.0;
    for (const auto& o : orders) sum_components += m.calculate_cost(o).total();
    REQUIRE_THAT(m.calculate_total_cost(orders), WithinAbs(sum_components, 1e-12));

    // Vector interface
    Eigen::VectorXd shares(2); shares << 100.0, -50.0;
    Eigen::VectorXd prices(2); prices << 100.0, 200.0;
    std::vector<std::string> tickers{"A", "B"};
    REQUIRE_THAT(m.calculate_total_cost(shares, prices, tickers), WithinAbs(sum_components, 1e-10));

    // Empty vector returns 0
    std::vector<TradeOrder> empty;
    REQUIRE_THAT(m.calculate_total_cost(empty), WithinAbs(0.0, 1e-12));
}
