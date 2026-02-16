// ============================================================================
// Implementation of Portfolio
// ============================================================================

#include "backtest/portfolio.hpp"

#include <iostream>
#include <sstream>

namespace portfolio {
namespace backtest {

// ============================================================================
// Position
// ============================================================================

double Position::unrealized_pnl() const {
    return market_value - (shares * cost_basis);
}

double Position::weight(double total_nav) const {
    if (total_nav <= 0.0) return 0.0;
    return market_value / total_nav;
}

// ============================================================================
// Portfolio - helpers
// ============================================================================

void Portfolio::build_ticker_index() {
    ticker_index_.clear();
    for (size_t i = 0; i < tickers_.size(); ++i) {
        ticker_index_[tickers_[i]] = static_cast<int>(i);
    }
}

int Portfolio::find_ticker_index(const std::string& ticker) const {
    auto it = ticker_index_.find(ticker);
    if (it == ticker_index_.end()) return -1;
    return it->second;
}

void Portfolio::validate_prices(const Eigen::VectorXd& prices) const {
    if (static_cast<size_t>(prices.size()) != tickers_.size()) {
        std::ostringstream msg;
        msg << "prices size (" << prices.size() << ") != num assets (" << tickers_.size() << ")";
        throw std::invalid_argument(msg.str());
    }
    for (int i = 0; i < prices.size(); ++i) {
        if (!(prices[i] > 0.0)) {
            std::ostringstream msg;
            msg << "price at index " << i << " is not positive: " << prices[i];
            throw std::invalid_argument(msg.str());
        }
    }
}

// ============================================================================
// Portfolio - lifecycle
// ============================================================================

Portfolio::Portfolio(double initial_capital, const std::vector<std::string>& tickers)
    : initial_capital_(initial_capital), cash_(initial_capital), tickers_(tickers) {
    if (!(initial_capital > 0.0)) {
        throw std::invalid_argument("initial_capital must be > 0");
    }
    if (tickers.empty()) {
        throw std::invalid_argument("tickers must not be empty");
    }
    positions_.resize(tickers_.size());
    for (size_t i = 0; i < tickers_.size(); ++i) {
        positions_[i].ticker = tickers_[i];
    }
    prices_ = Eigen::VectorXd::Zero(static_cast<int>(tickers_.size()));
    build_ticker_index();
    prev_nav_ = initial_capital_;
}

// ============================================================================
// Portfolio - queries
// ============================================================================

double Portfolio::nav() const {
    double invested = 0.0;
    for (size_t i = 0; i < positions_.size(); ++i) {
        invested += positions_[i].shares * prices_[static_cast<int>(i)];
    }
    return cash_ + invested;
}

double Portfolio::cash() const { return cash_; }

double Portfolio::invested_value() const {
    double invested = 0.0;
    for (size_t i = 0; i < positions_.size(); ++i) {
        invested += positions_[i].shares * prices_[static_cast<int>(i)];
    }
    return invested;
}

size_t Portfolio::num_assets() const { return tickers_.size(); }

const std::vector<std::string>& Portfolio::tickers() const { return tickers_; }

Eigen::VectorXd Portfolio::current_weights() const {
    double n = nav();
    Eigen::VectorXd w(static_cast<int>(tickers_.size()));
    if (!(n > 0.0)) {
        w.setZero();
        return w;
    }
    for (size_t i = 0; i < positions_.size(); ++i) {
        double mv = positions_[i].shares * prices_[static_cast<int>(i)];
        w[static_cast<int>(i)] = mv / n;
    }
    return w;
}

Eigen::VectorXd Portfolio::current_shares() const {
    Eigen::VectorXd s(static_cast<int>(positions_.size()));
    for (size_t i = 0; i < positions_.size(); ++i) s[static_cast<int>(i)] = positions_[i].shares;
    return s;
}

Eigen::VectorXd Portfolio::current_prices() const { return prices_; }

const Position& Portfolio::get_position(const std::string& ticker) const {
    int idx = find_ticker_index(ticker);
    if (idx < 0) throw std::invalid_argument("ticker not in universe: " + ticker);
    return positions_[static_cast<size_t>(idx)];
}

bool Portfolio::has_position(const std::string& ticker) const {
    return find_ticker_index(ticker) >= 0;
}

// ============================================================================
// Portfolio - price updates
// ============================================================================

void Portfolio::update_prices(const std::string& date, const Eigen::VectorXd& prices) {
    validate_prices(prices);
    current_date_ = date;
    prices_ = prices;
    for (size_t i = 0; i < positions_.size(); ++i) {
        positions_[i].market_value = positions_[i].shares * prices_[static_cast<int>(i)];
    }
}

// ============================================================================
// Portfolio - trading
// ============================================================================

Eigen::VectorXd Portfolio::set_target_weights(const Eigen::VectorXd& target_weights,
                                              const Eigen::VectorXd& prices) {
    if (static_cast<size_t>(target_weights.size()) != tickers_.size()) {
        throw std::invalid_argument("target_weights size mismatch");
    }
    validate_prices(prices);

    // compute current market values and desired dollars
    double current_nav = cash_;
    Eigen::VectorXd current_dollars(static_cast<int>(tickers_.size()));
    for (size_t i = 0; i < tickers_.size(); ++i) {
        current_dollars[static_cast<int>(i)] = positions_[i].shares * prices[static_cast<int>(i)];
        current_nav += current_dollars[static_cast<int>(i)];
    }

    Eigen::VectorXd desired_dollars = target_weights * current_nav;

    Eigen::VectorXd share_trades = Eigen::VectorXd::Zero(static_cast<int>(tickers_.size()));

    // SELL first
    for (size_t i = 0; i < tickers_.size(); ++i) {
        double delta = desired_dollars[static_cast<int>(i)] - current_dollars[static_cast<int>(i)];
        if (delta < 0.0) {
            double shares = delta / prices[static_cast<int>(i)];
            // shares is negative -> execute sell
            execute_trade(tickers_[i], shares, prices[static_cast<int>(i)]);
            share_trades[static_cast<int>(i)] = shares;
            current_dollars[static_cast<int>(i)] = positions_[i].shares * prices[static_cast<int>(i)];
        }
    }

    // BUY second
    for (size_t i = 0; i < tickers_.size(); ++i) {
        double delta = desired_dollars[static_cast<int>(i)] - current_dollars[static_cast<int>(i)];
        if (delta > 0.0) {
            double shares = delta / prices[static_cast<int>(i)];
            execute_trade(tickers_[i], shares, prices[static_cast<int>(i)]);
            share_trades[static_cast<int>(i)] += shares;
        }
    }

    // update stored prices and market values
    update_prices(current_date_, prices);
    return share_trades;
}

void Portfolio::execute_trade(const std::string& ticker, double shares, double price) {
    int idx = find_ticker_index(ticker);
    if (idx < 0) throw std::invalid_argument("ticker not in universe: " + ticker);
    if (!(price > 0.0)) throw std::invalid_argument("price must be > 0");

    double cost = shares * price; // positive for buy, negative for sell

    if (shares > 0.0) {
        if (cost > cash_ + 1e-12) {
            throw std::runtime_error("insufficient cash for buy");
        }
    } else if (shares < 0.0) {
        double abs_shares = -shares;
        if (abs_shares > positions_[static_cast<size_t>(idx)].shares + 1e-12) {
            throw std::runtime_error("attempt to sell more shares than held");
        }
    }

    // update cost basis for buys
    Position& pos = positions_[static_cast<size_t>(idx)];
    if (shares > 0.0) {
        double new_total_shares = pos.shares + shares;
        if (new_total_shares > 0.0) {
            pos.cost_basis = ((pos.shares * pos.cost_basis) + (shares * price)) / new_total_shares;
        } else {
            pos.cost_basis = 0.0;
        }
    } else if (shares < 0.0) {
        if (pos.shares + shares <= 0.0) {
            pos.cost_basis = 0.0;
        }
    }

    pos.shares += shares;
    pos.market_value = pos.shares * price;

    cash_ -= cost;
}

void Portfolio::deduct_costs(double amount) {
    if (!(amount >= 0.0)) throw std::invalid_argument("amount must be >= 0");
    if (amount > cash_ + 1e-12) throw std::runtime_error("insufficient cash to deduct costs");
    cash_ -= amount;
}

// ============================================================================
// Portfolio - snapshots & history
// ============================================================================

PortfolioSnapshot Portfolio::take_snapshot() const {
    PortfolioSnapshot s;
    s.date = current_date_;
    s.nav = nav();
    s.cash = cash_;
    s.invested = invested_value();
    s.tickers = tickers_;
    s.weights = current_weights();
    if (prev_nav_ > 0.0) {
        s.daily_return = (s.nav - prev_nav_) / prev_nav_;
    } else {
        s.daily_return = 0.0;
    }
    if (initial_capital_ > 0.0) {
        s.cumulative_return = (s.nav - initial_capital_) / initial_capital_;
    } else {
        s.cumulative_return = 0.0;
    }
    return s;
}

const std::vector<PortfolioSnapshot>& Portfolio::history() const { return history_; }

void Portfolio::record_snapshot() {
    PortfolioSnapshot s = take_snapshot();
    history_.push_back(s);
    prev_nav_ = s.nav;
}

// ============================================================================
// Portfolio - utilities
// ============================================================================

void Portfolio::print_summary() const {
    std::cout << "Date: " << current_date_ << " NAV: " << nav() << " Cash: " << cash_ << "\n";
}

void Portfolio::reset(double initial_capital) {
    if (!(initial_capital > 0.0)) throw std::invalid_argument("initial_capital must be > 0");
    initial_capital_ = initial_capital;
    cash_ = initial_capital;
    prices_.setZero();
    for (auto& p : positions_) {
        p.shares = 0.0;
        p.cost_basis = 0.0;
        p.market_value = 0.0;
    }
    history_.clear();
    prev_nav_ = initial_capital_;
}

} // namespace backtest
} // namespace portfolio
