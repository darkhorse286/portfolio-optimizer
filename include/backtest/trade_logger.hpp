#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include <Eigen/Dense>
#include "backtest/transaction_cost_model.hpp"

namespace portfolio {
namespace backtest {

struct TradeRecord {
    int trade_id;
    std::string date;
    std::string ticker;
    double shares;
    double price;
    double notional;
    double commission;
    double slippage;
    double market_impact;
    double total_cost;
    double pre_trade_weight;
    double post_trade_weight;
};

struct TradeSummary {
    int total_trades = 0;
    int buy_trades = 0;
    int sell_trades = 0;
    double total_notional = 0.0;
    double total_costs = 0.0;
    double avg_cost_per_trade = 0.0;
    int rebalance_count = 0;
    double turnover = 0.0;
};

class TradeLogger {
public:
    TradeLogger() = default;
    ~TradeLogger() = default;

    void log_trade(const TradeRecord& record);

    void log_rebalance(const std::string& date,
                       const std::vector<std::string>& tickers,
                       const Eigen::VectorXd& shares_traded,
                       const Eigen::VectorXd& prices,
                       const std::vector<TradeCost>& costs,
                       const Eigen::VectorXd& pre_weights,
                       const Eigen::VectorXd& post_weights);

    const std::vector<TradeRecord>& trades() const { return trades_; }
    std::vector<TradeRecord> trades_for_date(const std::string& date) const;
    std::vector<TradeRecord> trades_for_ticker(const std::string& ticker) const;
    TradeSummary get_summary() const;
    int num_trades() const { return static_cast<int>(trades_.size()); }

    void export_to_csv(const std::string& filepath) const;
    void print_summary() const;

    void clear() { trades_.clear(); next_trade_id_ = 0; rebalance_count_ = 0; }

private:
    std::vector<TradeRecord> trades_;
    int next_trade_id_ = 0;
    int rebalance_count_ = 0;
    double turnover_ = 0.0;
};

} // namespace backtest
} // namespace portfolio
