#pragma once

#include <string>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

namespace portfolio {
namespace backtest {

enum class RebalanceFrequency {
    DAILY,
    WEEKLY,
    MONTHLY,
    QUARTERLY,
    ANNUALLY
};

struct RebalanceConfig {
    RebalanceFrequency frequency = RebalanceFrequency::MONTHLY;
    double drift_threshold = 0.0;
    int min_days_between = 0;

    static RebalanceConfig from_json(const nlohmann::json& j);
    static RebalanceConfig from_string(const std::string& freq_str);
    static RebalanceFrequency parse_frequency(const std::string& freq_str);
};

class RebalanceScheduler {
public:
    explicit RebalanceScheduler(const RebalanceConfig& config);
    ~RebalanceScheduler() = default;

    bool should_rebalance(const std::string& date,
                          const Eigen::VectorXd& current_weights,
                          const Eigen::VectorXd& target_weights) const;

    bool is_calendar_trigger(const std::string& date) const;

    bool is_drift_trigger(const Eigen::VectorXd& current_weights,
                          const Eigen::VectorXd& target_weights) const;

    void record_rebalance(const std::string& date);
    void reset();

    const RebalanceConfig& config() const { return config_; }
    const std::string& last_rebalance_date() const { return last_rebalance_date_; }
    int rebalance_count() const { return rebalance_count_; }

private:
    RebalanceConfig config_;
    mutable std::string last_rebalance_date_;
    mutable int rebalance_count_;

    bool check_calendar(const std::string& date) const;
    bool check_min_days(const std::string& date) const;
    int days_between(const std::string& date1, const std::string& date2) const;
    int extract_year(const std::string& date) const;
    int extract_month(const std::string& date) const;
    int extract_day(const std::string& date) const;
    int day_of_week(const std::string& date) const;  // 0=Mon, 6=Sun
    long long days_since_epoch(const std::string& date) const;
};

} // namespace backtest
} // namespace portfolio
