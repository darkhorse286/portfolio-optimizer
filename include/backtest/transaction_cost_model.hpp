// transaction_cost_model.hpp
#pragma once

#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

namespace portfolio {
namespace backtest {

struct TradeOrder {
    std::string ticker;
    double shares{0.0};
    double price{0.0};
    double notional() const { return std::abs(shares * price); }
};

struct TradeCost {
    double commission{0.0};
    double slippage{0.0};
    double market_impact{0.0};
    double total() const { return commission + slippage + market_impact; }
};

struct TransactionCostConfig {
    double commission_rate{0.0};
    double slippage_bps{0.0};
    std::string market_impact_model{"none"};
    double market_impact_coeff{0.0};

    static TransactionCostConfig from_json(const nlohmann::json& j);
    static TransactionCostConfig default_config();
};

class TransactionCostModel {
public:
    explicit TransactionCostModel(const TransactionCostConfig& config);
    TransactionCostModel();
    ~TransactionCostModel() = default;

    TradeCost calculate_cost(const TradeOrder& order) const;
    double calculate_total_cost(const std::vector<TradeOrder>& orders) const;
    double calculate_total_cost(const Eigen::VectorXd& trade_shares,
                                 const Eigen::VectorXd& prices,
                                 const std::vector<std::string>& tickers) const;

    const TransactionCostConfig& config() const { return config_; }
    std::string get_name() const { return "TransactionCostModel"; }

    double commission_cost(double notional) const;
    double slippage_cost(double notional) const;
    double market_impact_cost(double notional, double price) const;

private:
    TransactionCostConfig config_;
    void validate_config() const;
};

} // namespace backtest
} // namespace portfolio
