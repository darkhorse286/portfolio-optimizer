#include "backtest/rebalance_scheduler.hpp"

#include <ctime>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cctype>

namespace portfolio {
namespace backtest {

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}

RebalanceConfig RebalanceConfig::from_json(const nlohmann::json& j) {
    RebalanceConfig cfg;
    if (j.contains("frequency")) {
        cfg.frequency = parse_frequency(j.at("frequency").get<std::string>());
    }
    if (j.contains("drift_threshold")) {
        cfg.drift_threshold = j.at("drift_threshold").get<double>();
    }
    if (j.contains("min_days_between")) {
        cfg.min_days_between = j.at("min_days_between").get<int>();
    }
    return cfg;
}

RebalanceConfig RebalanceConfig::from_string(const std::string& freq_str) {
    RebalanceConfig cfg;
    cfg.frequency = parse_frequency(freq_str);
    return cfg;
}

RebalanceFrequency RebalanceConfig::parse_frequency(const std::string& freq_str) {
    auto s = to_lower(freq_str);
    if (s == "daily" || s == "d") return RebalanceFrequency::DAILY;
    if (s == "weekly" || s == "w") return RebalanceFrequency::WEEKLY;
    if (s == "monthly" || s == "m") return RebalanceFrequency::MONTHLY;
    if (s == "quarterly" || s == "q") return RebalanceFrequency::QUARTERLY;
    if (s == "annually" || s == "annual" || s == "y" || s == "yearly") return RebalanceFrequency::ANNUALLY;
    throw std::invalid_argument("Invalid rebalance frequency: " + freq_str);
}

RebalanceScheduler::RebalanceScheduler(const RebalanceConfig& config)
    : config_(config), last_rebalance_date_(), rebalance_count_(0) {}

bool RebalanceScheduler::should_rebalance(const std::string& date,
                                         const Eigen::VectorXd& current_weights,
                                         const Eigen::VectorXd& target_weights) const {
    if (!check_min_days(date)) return false;
    if (check_calendar(date)) return true;
    return is_drift_trigger(current_weights, target_weights);
}

bool RebalanceScheduler::is_calendar_trigger(const std::string& date) const {
    return check_calendar(date);
}

bool RebalanceScheduler::is_drift_trigger(const Eigen::VectorXd& current_weights,
                                         const Eigen::VectorXd& target_weights) const {
    if (config_.drift_threshold <= 0.0) return false;
    if (current_weights.size() != target_weights.size())
        throw std::invalid_argument("Weight vectors must have same size");
    for (int i = 0; i < current_weights.size(); ++i) {
        if (std::abs(current_weights[i] - target_weights[i]) > config_.drift_threshold) return true;
    }
    return false;
}

void RebalanceScheduler::record_rebalance(const std::string& date) {
    last_rebalance_date_ = date;
    ++rebalance_count_;
}

void RebalanceScheduler::reset() {
    last_rebalance_date_.clear();
    rebalance_count_ = 0;
}

bool RebalanceScheduler::check_calendar(const std::string& date) const {
    int year = extract_year(date);
    int month = extract_month(date);
    int day = extract_day(date);

    switch (config_.frequency) {
        case RebalanceFrequency::DAILY:
            return true;
        case RebalanceFrequency::WEEKLY: {
            int dow = day_of_week(date);
            if (dow == 0) return true; // Monday
            if (last_rebalance_date_.empty()) return false;
            long long w1 = days_since_epoch(date) / 7;
            long long w0 = days_since_epoch(last_rebalance_date_) / 7;
            return w1 != w0;
        }
        case RebalanceFrequency::MONTHLY: {
            // Trigger only on the first day of a new month
            if (day != 1) return false;
            if (last_rebalance_date_.empty()) return true;
            int last_m = extract_month(last_rebalance_date_);
            return month != last_m;
        }
        case RebalanceFrequency::QUARTERLY: {
            auto is_quarter_month = [](int m){ return m == 1 || m == 4 || m == 7 || m == 10; };
            if (!is_quarter_month(month)) return false;
            if (last_rebalance_date_.empty()) return day == 1;
            int last_m = extract_month(last_rebalance_date_);
            int last_q = (last_m - 1) / 3;
            int q = (month - 1) / 3;
            return q != last_q;
        }
        case RebalanceFrequency::ANNUALLY: {
            if (last_rebalance_date_.empty()) return (month == 1 && day == 1);
            int last_y = extract_year(last_rebalance_date_);
            return year != last_y;
        }
    }
    return false;
}

bool RebalanceScheduler::check_min_days(const std::string& date) const {
    if (config_.min_days_between <= 0) return true;
    if (last_rebalance_date_.empty()) return true;
    return days_between(last_rebalance_date_, date) >= config_.min_days_between;
}

int RebalanceScheduler::days_between(const std::string& date1, const std::string& date2) const {
    std::tm tm1 = {};
    std::tm tm2 = {};
    tm1.tm_year = extract_year(date1) - 1900;
    tm1.tm_mon = extract_month(date1) - 1;
    tm1.tm_mday = extract_day(date1);
    tm1.tm_hour = 12;
    tm2.tm_year = extract_year(date2) - 1900;
    tm2.tm_mon = extract_month(date2) - 1;
    tm2.tm_mday = extract_day(date2);
    tm2.tm_hour = 12;
    std::time_t t1 = std::mktime(&tm1);
    std::time_t t2 = std::mktime(&tm2);
    if (t1 == (std::time_t)-1 || t2 == (std::time_t)-1) return 0;
    double diff = std::difftime(t2, t1);
    return static_cast<int>(std::llround(std::llround(diff) / 86400.0));
}

int RebalanceScheduler::extract_year(const std::string& date) const {
    return std::stoi(date.substr(0,4));
}

int RebalanceScheduler::extract_month(const std::string& date) const {
    return std::stoi(date.substr(5,2));
}

int RebalanceScheduler::extract_day(const std::string& date) const {
    return std::stoi(date.substr(8,2));
}

int RebalanceScheduler::day_of_week(const std::string& date) const {
    std::tm tm = {};
    tm.tm_year = extract_year(date) - 1900;
    tm.tm_mon = extract_month(date) - 1;
    tm.tm_mday = extract_day(date);
    tm.tm_hour = 12;
    std::mktime(&tm);
    // tm_wday: 0=Sun, 1=Mon ... 6=Sat
    int w = tm.tm_wday;
    int result = (w + 6) % 7; // convert to 0=Mon, 6=Sun
    return result;
}

long long RebalanceScheduler::days_since_epoch(const std::string& date) const {
    std::tm tm = {};
    tm.tm_year = extract_year(date) - 1900;
    tm.tm_mon = extract_month(date) - 1;
    tm.tm_mday = extract_day(date);
    tm.tm_hour = 12;
    std::time_t t = std::mktime(&tm);
    if (t == (std::time_t)-1) return 0;
    return static_cast<long long>(t / 86400);
}

} // namespace backtest
} // namespace portfolio
