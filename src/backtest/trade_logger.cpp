/**
 * @file trade_logger.cpp
 * @brief Implementation of TradeLogger
 */

#include "backtest/trade_logger.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace portfolio {
namespace backtest {

void TradeLogger::log_trade(const TradeRecord &record)
{
    TradeRecord r = record;
    r.trade_id = next_trade_id_++;
    trades_.push_back(r);
}

void TradeLogger::log_rebalance(const std::string &date,
                                const std::vector<std::string> &tickers,
                                const Eigen::VectorXd &shares_traded,
                                const Eigen::VectorXd &prices,
                                const std::vector<TradeCost> &costs,
                                const Eigen::VectorXd &pre_weights,
                                const Eigen::VectorXd &post_weights)
{
    if ((int)tickers.size() != shares_traded.size() ||
        (int)tickers.size() != prices.size() ||
        (int)tickers.size() != (int)costs.size() ||
        (int)tickers.size() != pre_weights.size() ||
        (int)tickers.size() != post_weights.size())
    {
        throw std::invalid_argument("Mismatched rebalance vector sizes");
    }

    // compute turnover contribution: sum |delta| / 2
    double abs_sum = 0.0;
    for (int i = 0; i < pre_weights.size(); ++i)
    {
        abs_sum += std::abs(post_weights[i] - pre_weights[i]);
    }
    double turnover_contribution = abs_sum / 2.0;

    rebalance_count_ += 1;

    // record trades where shares changed (non-zero within tolerance)
    for (int i = 0; i < (int)tickers.size(); ++i)
    {
        double s = shares_traded[i];
        if (std::abs(s) <= 1e-10)
            continue;

        TradeRecord r;
        r.trade_id = next_trade_id_++;
        r.date = date;
        r.ticker = tickers[i];
        r.shares = s;
        r.price = prices[i];
        r.notional = s * prices[i];
        r.commission = costs[i].commission;
        r.slippage = costs[i].slippage;
        r.market_impact = costs[i].market_impact;
        r.total_cost = costs[i].total();
        r.pre_trade_weight = pre_weights[i];
        r.post_trade_weight = post_weights[i];

        trades_.push_back(r);
    }

    // accumulate turnover contribution in member variable
    turnover_ += turnover_contribution;
}

std::vector<TradeRecord> TradeLogger::trades_for_date(const std::string &date) const
{
    std::vector<TradeRecord> out;
    for (const auto &t : trades_)
    {
        if (t.date == date)
            out.push_back(t);
    }
    return out;
}

std::vector<TradeRecord> TradeLogger::trades_for_ticker(const std::string &ticker) const
{
    std::vector<TradeRecord> out;
    for (const auto &t : trades_)
    {
        if (t.ticker == ticker)
            out.push_back(t);
    }
    return out;
}

TradeSummary TradeLogger::get_summary() const
{
    TradeSummary s;
    s.total_trades = static_cast<int>(trades_.size());
    s.rebalance_count = rebalance_count_;

    for (const auto &t : trades_)
    {
        if (t.shares > 0)
            ++s.buy_trades;
        else if (t.shares < 0)
            ++s.sell_trades;

        s.total_notional += std::abs(t.notional);
        s.total_costs += t.total_cost;
    }

    if (s.total_trades > 0)
        s.avg_cost_per_trade = s.total_costs / s.total_trades;

    // Turnover must be computed from rebalance events. We approximate turnover by
    // summing |post - pre| / 2 across trades by grouping by rebalance date.
    // Simpler approach: iterate trades and accumulate |post-pre|, but that double-counts
    // per rebalance if multiple trades share same pre/post vectors. Since we don't
    // store rebalances separately, compute turnover by summing |post-pre| across all trades
    // and dividing by 2. This matches the spec if each rebalance included all assets.
    s.turnover = turnover_;

    return s;
}

void TradeLogger::export_to_csv(const std::string &filepath) const
{
    std::filesystem::path path(filepath);
    if (path.has_parent_path())
    {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream file(filepath);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file for writing: " + filepath);
    }

    file << "trade_id,date,ticker,shares,price,notional,commission,slippage,market_impact,total_cost,pre_trade_weight,post_trade_weight\n";
    file << std::fixed << std::setprecision(8);

    for (const auto &t : trades_)
    {
        file << t.trade_id << ","
             << t.date << ","
             << t.ticker << ","
             << t.shares << ","
             << t.price << ","
             << t.notional << ","
             << t.commission << ","
             << t.slippage << ","
             << t.market_impact << ","
             << t.total_cost << ","
             << t.pre_trade_weight << ","
             << t.post_trade_weight << "\n";
    }

    file.close();
}

void TradeLogger::print_summary() const
{
    auto s = get_summary();
    std::cout << "\n=== Trade Log Summary ===\n";
    std::cout << "Total trades: " << s.total_trades << "\n";
    std::cout << "Buys: " << s.buy_trades << "  Sells: " << s.sell_trades << "\n";
    std::cout << "Total notional: " << s.total_notional << "\n";
    std::cout << "Total costs: " << s.total_costs << "\n";
    std::cout << "Average cost/trade: " << s.avg_cost_per_trade << "\n";
    std::cout << "Rebalance events: " << s.rebalance_count << "\n";
    std::cout << "Turnover: " << s.turnover << "\n";
    std::cout << "==========================\n";
}

} // namespace backtest
} // namespace portfolio
