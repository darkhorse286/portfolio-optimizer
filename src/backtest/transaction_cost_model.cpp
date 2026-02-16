#include "backtest/transaction_cost_model.hpp"

#include <cmath>
#include <sstream>

namespace portfolio {
namespace backtest {

TransactionCostConfig TransactionCostConfig::default_config() {
    TransactionCostConfig cfg;
    cfg.commission_rate = 0.0;
    cfg.slippage_bps = 0.0;
    cfg.market_impact_model = "none";
    cfg.market_impact_coeff = 0.0;
    return cfg;
}

TransactionCostConfig TransactionCostConfig::from_json(const nlohmann::json& j) {
    TransactionCostConfig cfg = default_config();
    if (j.is_object()) {
        cfg.commission_rate = j.value("commission_rate", cfg.commission_rate);
        cfg.slippage_bps = j.value("slippage_bps", cfg.slippage_bps);
        cfg.market_impact_model = j.value("market_impact_model", cfg.market_impact_model);
        cfg.market_impact_coeff = j.value("market_impact_coeff", cfg.market_impact_coeff);
    }
    return cfg;
}

TransactionCostModel::TransactionCostModel(const TransactionCostConfig& config)
    : config_(config) {
    validate_config();
}

TransactionCostModel::TransactionCostModel()
    : config_(TransactionCostConfig::default_config()) {
}

void TransactionCostModel::validate_config() const {
    if (config_.commission_rate < 0.0) {
        std::ostringstream ss; ss << config_.commission_rate;
        throw std::invalid_argument("Expected positive value for parameter 'commission_rate', got: " + ss.str());
    }
    if (config_.slippage_bps < 0.0) {
        std::ostringstream ss; ss << config_.slippage_bps;
        throw std::invalid_argument("Expected positive value for parameter 'slippage_bps', got: " + ss.str());
    }
    if (config_.market_impact_coeff < 0.0) {
        std::ostringstream ss; ss << config_.market_impact_coeff;
        throw std::invalid_argument("Expected positive value for parameter 'market_impact_coeff', got: " + ss.str());
    }
    if (!(config_.market_impact_model == "none" || config_.market_impact_model == "linear" || config_.market_impact_model == "sqrt")) {
        throw std::invalid_argument("Expected one of 'none','linear','sqrt' for parameter 'market_impact_model', got: " + config_.market_impact_model);
    }
}

TradeCost TransactionCostModel::calculate_cost(const TradeOrder& order) const {
    if (order.price <= 0.0) {
        std::ostringstream ss; ss << order.price;
        throw std::invalid_argument("Expected positive value for parameter 'price', got: " + ss.str());
    }
    double notional = order.notional();
    double commission = commission_cost(notional);
    double slippage = slippage_cost(notional);
    double impact = market_impact_cost(notional, order.price);
    // ensure non-negative
    if (commission < 0.0) commission = 0.0;
    if (slippage < 0.0) slippage = 0.0;
    if (impact < 0.0) impact = 0.0;
    return TradeCost{commission, slippage, impact};
}

double TransactionCostModel::calculate_total_cost(const std::vector<TradeOrder>& orders) const {
    double sum = 0.0;
    for (const auto& o : orders) {
        sum += calculate_cost(o).total();
    }
    return sum;
}

double TransactionCostModel::calculate_total_cost(const Eigen::VectorXd& trade_shares,
                                                  const Eigen::VectorXd& prices,
                                                  const std::vector<std::string>& tickers) const {
    if (trade_shares.size() != prices.size() || static_cast<size_t>(trade_shares.size()) != tickers.size()) {
        throw std::invalid_argument("Expected equal sizes for trade_shares, prices and tickers");
    }
    double sum = 0.0;
    for (int i = 0; i < trade_shares.size(); ++i) {
        double price = prices[i];
        if (price <= 0.0) {
            std::ostringstream ss; ss << price;
            throw std::invalid_argument("Expected positive value for parameter 'price', got: " + ss.str());
        }
        TradeOrder o;
        o.ticker = tickers[static_cast<size_t>(i)];
        o.shares = trade_shares[i];
        o.price = price;
        sum += calculate_cost(o).total();
    }
    return sum;
}

double TransactionCostModel::commission_cost(double notional) const {
    return std::abs(notional) * config_.commission_rate;
}

double TransactionCostModel::slippage_cost(double notional) const {
    return std::abs(notional) * (config_.slippage_bps / 10000.0);
}

double TransactionCostModel::market_impact_cost(double notional, double price) const {
    double abs_notional = std::abs(notional);
    if (config_.market_impact_model == "none") return 0.0;
    if (config_.market_impact_model == "linear") {
        return abs_notional * config_.market_impact_coeff;
    }
    // sqrt model: use a sub-linear form proportional to sqrt(notional * price)
    if (config_.market_impact_model == "sqrt") {
        if (price <= 0.0) return 0.0;
        double factor = std::sqrt(abs_notional * price);
        return config_.market_impact_coeff * factor;
    }
    return 0.0;
}

} // namespace backtest
} // namespace portfolio
