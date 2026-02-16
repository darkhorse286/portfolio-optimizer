// SPDX-License-Identifier: MIT
#ifndef PORTFOLIO_BACKTEST_PORTFOLIO_HPP
#define PORTFOLIO_BACKTEST_PORTFOLIO_HPP

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <Eigen/Dense>

namespace portfolio {
namespace backtest {

/**
 * @struct Position
 * @brief Represents a single asset position
 */
struct Position {
    std::string ticker;       ///< Asset identifier
    double shares = 0.0;      ///< Number of shares held
    double cost_basis = 0.0;  ///< Average cost per share
    double market_value = 0.0;///< Current market value

    double unrealized_pnl() const;
    double weight(double total_nav) const;
};

/**
 * @struct PortfolioSnapshot
 * @brief Point-in-time view of portfolio state
 */
struct PortfolioSnapshot {
    std::string date;                        ///< Snapshot date (YYYY-MM-DD)
    double nav = 0.0;                        ///< Net asset value
    double cash = 0.0;                       ///< Cash balance
    double invested = 0.0;                   ///< Total invested value
    Eigen::VectorXd weights;                 ///< Current asset weights
    std::vector<std::string> tickers;        ///< Asset names (same order as weights)
    double daily_return = 0.0;               ///< Return since last snapshot
    double cumulative_return = 0.0;          ///< Return since inception
};

/**
 * @class Portfolio
 * @brief Manages portfolio state: positions, cash, NAV, and trade execution
 */
class Portfolio {
public:
    explicit Portfolio(double initial_capital,
                       const std::vector<std::string>& tickers);

    ~Portfolio() = default;

    // -- State queries
    double nav() const;
    double cash() const;
    double invested_value() const;
    size_t num_assets() const;
    const std::vector<std::string>& tickers() const;

    Eigen::VectorXd current_weights() const;
    Eigen::VectorXd current_shares() const;
    Eigen::VectorXd current_prices() const;
    const Position& get_position(const std::string& ticker) const;
    bool has_position(const std::string& ticker) const;

    // -- Price updates
    void update_prices(const std::string& date, const Eigen::VectorXd& prices);

    // -- Trading
    Eigen::VectorXd set_target_weights(const Eigen::VectorXd& target_weights,
                                      const Eigen::VectorXd& prices);

    void execute_trade(const std::string& ticker, double shares, double price);
    void deduct_costs(double amount);

    // -- Snapshots
    PortfolioSnapshot take_snapshot() const;

    // -- History
    const std::vector<PortfolioSnapshot>& history() const;
    void record_snapshot();

    // -- Utilities
    void print_summary() const;
    void reset(double initial_capital);

private:
    double initial_capital_;
    double cash_;
    std::string current_date_;
    std::vector<std::string> tickers_;
    std::map<std::string, int> ticker_index_;
    std::vector<Position> positions_;
    Eigen::VectorXd prices_;
    std::vector<PortfolioSnapshot> history_;
    double prev_nav_ = 0.0;

    void build_ticker_index();
    int find_ticker_index(const std::string& ticker) const;
    void validate_prices(const Eigen::VectorXd& prices) const;
};

} // namespace backtest
} // namespace portfolio

#endif // PORTFOLIO_BACKTEST_PORTFOLIO_HPP
