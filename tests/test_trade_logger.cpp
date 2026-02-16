#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "backtest/trade_logger.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <filesystem>

using namespace portfolio::backtest;

TEST_CASE("TradeLogger logging", "[TradeLogger]") {
    TradeLogger logger;

    // single trade
    TradeRecord r;
    r.trade_id = -1;
    r.date = "2020-01-01";
    r.ticker = "AAA";
    r.shares = 5.0;
    r.price = 10.0;
    r.notional = r.shares * r.price;
    r.commission = 0.5;
    r.slippage = 0.2;
    r.market_impact = 0.1;
    r.total_cost = r.commission + r.slippage + r.market_impact;
    r.pre_trade_weight = 0.0;
    r.post_trade_weight = 0.05;

    logger.log_trade(r);
    REQUIRE(logger.num_trades() == 1);
    REQUIRE(logger.trades().front().trade_id == 0);

    // rebalance batch: one zero-share should be skipped
    std::vector<std::string> tickers = {"A","B","C"};
    Eigen::VectorXd shares(3);
    shares << 10.0, 0.0, -5.0;
    Eigen::VectorXd prices(3);
    prices << 1.0, 2.0, 3.0;
    std::vector<TradeCost> costs(3);
    costs[0].commission = 0.1; costs[0].slippage = 0.0; costs[0].market_impact = 0.0;
    costs[1].commission = 0.0; costs[1].slippage = 0.0; costs[1].market_impact = 0.0;
    costs[2].commission = 0.2; costs[2].slippage = 0.1; costs[2].market_impact = 0.0;
    Eigen::VectorXd pre_w(3); pre_w << 0.1, 0.2, 0.3;
    Eigen::VectorXd post_w(3); post_w << 0.2, 0.2, 0.1;

    logger.log_rebalance("2020-01-02", tickers, shares, prices, costs, pre_w, post_w);
    REQUIRE(logger.num_trades() == 3); // 1 + 2 (middle skipped)
    // trade ids should auto-increment
    auto all = logger.trades();
    REQUIRE(all[0].trade_id == 0);
    REQUIRE(all[1].trade_id == 1);
    REQUIRE(all[2].trade_id == 2);
}

TEST_CASE("TradeLogger queries", "[TradeLogger]") {
    TradeLogger logger;
    TradeRecord r1{0, "2020-01-01", "AAA", 5.0, 10.0, 50.0, 0.5,0.2,0.1,0.8,0.0,0.05};
    TradeRecord r2{1, "2020-01-02", "BBB", -2.0, 20.0, -40.0, 0.2,0.0,0.0,0.2,0.1,0.0};
    logger.log_trade(r1);
    logger.log_trade(r2);

    auto d1 = logger.trades_for_date("2020-01-01");
    REQUIRE(d1.size() == 1);
    REQUIRE(d1.front().ticker == "AAA");

    auto tbb = logger.trades_for_ticker("BBB");
    REQUIRE(tbb.size() == 1);
    REQUIRE(tbb.front().date == "2020-01-02");

    auto none = logger.trades_for_date("1999-01-01");
    REQUIRE(none.empty());
}

TEST_CASE("TradeLogger summary", "[TradeLogger]") {
    TradeLogger logger;
    TradeRecord r1{0, "2020-01-01", "AAA", 5.0, 10.0, 50.0, 0.5,0.2,0.1,0.8,0.0,0.05};
    TradeRecord r2{1, "2020-01-02", "BBB", -2.0, 20.0, -40.0, 0.2,0.0,0.0,0.2,0.1,0.0};
    logger.log_trade(r1);
    logger.log_trade(r2);

    // add one rebalance to exercise turnover calculation
    std::vector<std::string> tickers = {"A","B"};
    Eigen::VectorXd shares(2); shares << 1.0, -1.0;
    Eigen::VectorXd prices(2); prices << 1.0, 2.0;
    std::vector<TradeCost> costs(2);
    Eigen::VectorXd pre_w(2); pre_w << 0.1, 0.2;
    Eigen::VectorXd post_w(2); post_w << 0.2, 0.1;
    logger.log_rebalance("2020-01-03", tickers, shares, prices, costs, pre_w, post_w);

    auto s = logger.get_summary();
    REQUIRE(s.total_trades == logger.num_trades());
    REQUIRE(s.buy_trades == 2); // r1 and rebalance buy
    REQUIRE(s.sell_trades == 2); // r2 and rebalance sell
    REQUIRE(s.total_notional == Catch::Approx(std::abs(r1.notional) + std::abs(r2.notional) + std::abs(1.0*1.0) + std::abs(-1.0*2.0)));
    REQUIRE(s.total_costs == Catch::Approx(r1.total_cost + r2.total_cost));

    // turnover: single rebalance contributed (|0.2-0.1| + |0.1-0.2|)/2 = (0.1+0.1)/2 = 0.1
    REQUIRE(s.turnover == Catch::Approx(0.1));
}

TEST_CASE("TradeLogger export", "[TradeLogger]") {
    TradeLogger logger;
    TradeRecord r{0, "2020-01-01", "AAA", 5.0, 10.0, 50.0, 0.5,0.2,0.1,0.8,0.0,0.05};
    logger.log_trade(r);

    const std::string out = "build/tmp/trades_test/trades.csv";
    // ensure parent dir removed first
    std::filesystem::remove_all(std::filesystem::path(out).parent_path());

    REQUIRE_NOTHROW(logger.export_to_csv(out));
    REQUIRE(std::filesystem::exists(out));

    std::ifstream f(out);
    REQUIRE(f.is_open());
    std::string header;
    std::getline(f, header);
    REQUIRE(header == "trade_id,date,ticker,shares,price,notional,commission,slippage,market_impact,total_cost,pre_trade_weight,post_trade_weight");

    std::string row;
    std::getline(f, row);
    REQUIRE_FALSE(row.empty());
}
