// SPDX-License-Identifier: MIT
#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "data/market_data.hpp"
#include "data/data_loader.hpp"
#include "backtest/portfolio.hpp"
#include "backtest/rebalance_scheduler.hpp"
#include "backtest/trade_logger.hpp"
#include "backtest/transaction_cost_model.hpp"
#include "risk/risk_model_factory.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
// forward-declarations above to avoid duplicate symbol definitions

namespace portfolio {
namespace backtest {

struct BacktestResult {
    bool success = false;
    std::string message;

    std::vector<PortfolioSnapshot> snapshots;
    std::vector<double> nav_series;
    std::vector<double> return_series;
    std::vector<std::string> dates;

    TradeSummary trade_summary;
    std::vector<TradeRecord> trades;

    double total_return() const;
    double annualized_return(int trading_days_per_year = 252) const;
    double annualized_volatility(int trading_days_per_year = 252) const;
    double sharpe_ratio(double risk_free_rate = 0.02,
                        int trading_days_per_year = 252) const;
    double max_drawdown() const;

    void print_summary() const;
    void export_nav_to_csv(const std::string& filepath) const;
};

struct BacktestParams {
    double initial_capital = 1000000.0;
    int lookback_window = 252;
    int min_history = 60;
    RebalanceConfig rebalance;
    TransactionCostConfig transaction_costs;
    optimizer::OptimizationConstraints constraints;
    std::string risk_model_type = "ewma";
    double risk_free_rate = 0.02;
    bool verbose = false;

    static BacktestParams from_config(const portfolio::PortfolioConfig& config);
};

class BacktestEngine {
public:
    explicit BacktestEngine(const BacktestParams& params);
    ~BacktestEngine() = default;

    BacktestResult run(const MarketData& market_data);
    BacktestResult run(const MarketData& market_data,
                       optimizer::MeanVarianceOptimizer& optimizer);

    const BacktestParams& params() const { return params_; }

private:
    BacktestParams params_;

    Eigen::MatrixXd extract_window(const MarketData& market_data,
                                    int date_index) const;

    Eigen::MatrixXd estimate_risk(const Eigen::MatrixXd& returns) const;

    Eigen::VectorXd estimate_returns(const Eigen::MatrixXd& returns) const;

    Eigen::VectorXd optimize(const Eigen::VectorXd& expected_returns,
                              const Eigen::MatrixXd& covariance,
                              const Eigen::VectorXd& current_weights,
                              optimizer::MeanVarianceOptimizer& optimizer) const;
};

} // namespace backtest
} // namespace portfolio
