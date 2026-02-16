#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "backtest/rebalance_scheduler.hpp"
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

using namespace portfolio::backtest;

TEST_CASE("RebalanceScheduler calendar triggers", "[RebalanceScheduler]") {
    RebalanceConfig cfg;
    cfg.min_days_between = 0;

    // DAILY
    cfg.frequency = RebalanceFrequency::DAILY;
    RebalanceScheduler s_daily(cfg);
    Eigen::VectorXd a(1); a << 0.5;
    Eigen::VectorXd b(1); b << 0.5;
    REQUIRE(s_daily.should_rebalance("2022-02-02", a, b));

    // WEEKLY: Monday triggers
    cfg.frequency = RebalanceFrequency::WEEKLY;
    RebalanceScheduler s_weekly(cfg);
    // 2021-11-01 is Monday
    REQUIRE(s_weekly.is_calendar_trigger("2021-11-01"));
    // 2021-11-03 is Wednesday -> not a Monday and without prior rebalance should be false
    REQUIRE(!s_weekly.is_calendar_trigger("2021-11-03"));

    // MONTHLY
    cfg.frequency = RebalanceFrequency::MONTHLY;
    RebalanceScheduler s_monthly(cfg);
    // simulate last rebalance in previous month
    s_monthly.record_rebalance("2021-10-29");
    REQUIRE(s_monthly.is_calendar_trigger("2021-11-01"));
    REQUIRE(!s_monthly.is_calendar_trigger("2021-11-15"));

    // QUARTERLY
    cfg.frequency = RebalanceFrequency::QUARTERLY;
    RebalanceScheduler s_quarter(cfg);
    s_quarter.record_rebalance("2021-12-31");
    REQUIRE(s_quarter.is_calendar_trigger("2022-01-03"));

    // ANNUALLY
    cfg.frequency = RebalanceFrequency::ANNUALLY;
    RebalanceScheduler s_annual(cfg);
    s_annual.record_rebalance("2021-12-31");
    REQUIRE(s_annual.is_calendar_trigger("2022-01-03"));
}

TEST_CASE("RebalanceScheduler drift triggers", "[RebalanceScheduler]") {
    RebalanceConfig cfg;
    cfg.frequency = RebalanceFrequency::DAILY;
    cfg.drift_threshold = 0.05;
    RebalanceScheduler s(cfg);

    Eigen::VectorXd t(3); t << 0.4, 0.3, 0.3;
    Eigen::VectorXd same = t;
    REQUIRE(!s.is_drift_trigger(same, t));

    Eigen::VectorXd within(3); within << 0.42, 0.28, 0.30;
    REQUIRE(!s.is_drift_trigger(within, t));

    Eigen::VectorXd exceed(3); exceed << 0.5, 0.25, 0.25;
    REQUIRE(s.is_drift_trigger(exceed, t));

    Eigen::VectorXd single(1); single << 0.2;
    Eigen::VectorXd single_t(1); single_t << 0.0;
    cfg.drift_threshold = 0.1;
    RebalanceScheduler s2(cfg);
    REQUIRE(s2.is_drift_trigger(single, single_t));

    cfg.drift_threshold = 0.0;
    RebalanceScheduler s3(cfg);
    REQUIRE(!s3.is_drift_trigger(exceed, t));
}

TEST_CASE("RebalanceScheduler combined triggers", "[RebalanceScheduler]") {
    RebalanceConfig cfg;
    cfg.frequency = RebalanceFrequency::MONTHLY;
    cfg.drift_threshold = 0.02;
    cfg.min_days_between = 5;
    RebalanceScheduler s(cfg);

    Eigen::VectorXd t(2); t << 0.5, 0.5;
    Eigen::VectorXd drift(2); drift << 0.55, 0.45;

    // Not a calendar date but drift exceeds threshold -> rebalance
    REQUIRE(s.should_rebalance("2021-11-15", drift, t));

    // Min days between honored
    s.record_rebalance("2021-11-10");
    // now asking within 5 days should return false even though drift present
    REQUIRE(!s.should_rebalance("2021-11-12", drift, t));
}

TEST_CASE("RebalanceScheduler configuration", "[RebalanceScheduler]") {
    REQUIRE(RebalanceConfig::parse_frequency("daily") == RebalanceFrequency::DAILY);
    REQUIRE(RebalanceConfig::parse_frequency("Weekly") == RebalanceFrequency::WEEKLY);
    REQUIRE(RebalanceConfig::parse_frequency("MONTHLY") == RebalanceFrequency::MONTHLY);
    REQUIRE(RebalanceConfig::parse_frequency("quarterly") == RebalanceFrequency::QUARTERLY);
    REQUIRE(RebalanceConfig::parse_frequency("annually") == RebalanceFrequency::ANNUALLY);

    nlohmann::json j = { {"frequency","monthly"}, {"drift_threshold", 0.05}, {"min_days_between", 3} };
    RebalanceConfig cfg = RebalanceConfig::from_json(j);
    REQUIRE(cfg.frequency == RebalanceFrequency::MONTHLY);
    REQUIRE(cfg.drift_threshold == Catch::Approx(0.05));
    REQUIRE(cfg.min_days_between == 3);

    REQUIRE_THROWS_AS(RebalanceConfig::parse_frequency("bogus"), std::invalid_argument);
}
