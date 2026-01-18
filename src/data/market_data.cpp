/**
 * @file market_data.cpp
 * @brief Implementation of MarketData class
 */

#include "data/market_data.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>

namespace portfolio
{

    // ============================================================================
    // Constructors
    // ============================================================================

    MarketData::MarketData(size_t num_dates, size_t num_assets) : prices_(num_dates, num_assets)
    {
        prices_.setZero();
    }

    MarketData::MarketData(const Eigen::MatrixXd &prices,
                           const std::vector<std::string> &dates,
                           const std::vector<std::string> &tickers)
        : prices_(prices), dates_(dates), tickers_(tickers)
    {

        // Validate dimensions
        if (prices_.rows() != static_cast<int>(dates_.size()))
        {
            throw std::invalid_argument("Price matrix rows must match dates vector size");
        }
        if (prices_.cols() != static_cast<int>(tickers_.size()))
        {
            throw std::invalid_argument("Price matrix columns must match tickers vector size");
        }

        build_index_maps();
    }

    // ============================================================================
    // Data Access Methods
    // ============================================================================

    Eigen::VectorXd MarketData::get_prices(const std::string &ticker) const
    {
        int idx = find_ticker_index(ticker);
        if (idx < 0)
        {
            throw std::invalid_argument("Ticker not found: " + ticker);
        }
        return prices_.col(idx);
    }

    Eigen::VectorXd MarketData::get_prices_at_date(const std::string &date) const
    {
        int idx = find_date_index(date);
        if (idx < 0)
        {
            throw std::invalid_argument("Date not found: " + date);
        }
        return prices_.row(idx);
    }

    double MarketData::get_price(const std::string &ticker, const std::string &date) const
    {
        int ticker_idx = find_ticker_index(ticker);
        int date_idx = find_date_index(date);

        if (ticker_idx < 0)
        {
            throw std::invalid_argument("Ticker not found: " + ticker);
        }
        if (date_idx < 0)
        {
            throw std::invalid_argument("Date not found: " + date);
        }

        return prices_(date_idx, ticker_idx);
    }

    // ==========================
    // Data Modification Methods
    // ==========================

    void MarketData::set_prices(const Eigen::MatrixXd &prices)
    {
        if (prices.rows() != prices_.rows() || prices.cols() != prices_.cols())
        {
            throw std::invalid_argument("New price matrix must have same dimensions");
        }
        prices_ = prices;
    }

    void MarketData::set_price(const std::string &ticker, const std::string &date, double price)
    {
        int ticker_idx = find_ticker_index(ticker);
        int date_idx = find_date_index(date);

        if (ticker_idx < 0)
        {
            throw std::invalid_argument("Ticker not found: " + ticker);
        }
        if (date_idx < 0)
        {
            throw std::invalid_argument("Date not found: " + date);
        }

        prices_(date_idx, ticker_idx) = price;
    }

    void MarketData::set_dates(const std::vector<std::string> &dates)
    {
        if (dates.size() != static_cast<size_t>(prices_.rows()))
        {
            throw std::invalid_argument("Dates vector size must match price matrix rows");
        }
        dates_ = dates;
        build_index_maps();
    }

    void MarketData::set_tickers(const std::vector<std::string> &tickers)
    {
        if (tickers.size() != static_cast<size_t>(prices_.cols()))
        {
            throw std::invalid_argument("Tickers vector size must match price matrix columns");
        }
        tickers_ = tickers;
        build_index_maps();
    }

    // ===========================
    // Return Calculation Methods
    // ===========================

    Eigen::MatrixXd MarketData::calculate_returns(ReturnType type) const
    {
        if (prices_.rows() < 2)
        {
            throw std::runtime_error("Need at least 2 price observations to calculate returns");
        }

        Eigen::MatrixXd returns(prices_.rows() - 1, prices_.cols());

        for (int i = 0; i < prices_.rows() - 1; ++i)
        {
            for (int j = 0; j < prices_.cols(); ++j)
            {
                double p_t = prices_(i + 1, j);
                double p_tm1 = prices_(i, j);

                if (std::isnan(p_t) || std::isnan(p_tm1) || p_tm1 == 0.0)
                {
                    returns(i, j) = std::numeric_limits<double>::quiet_NaN();
                }
                else
                {
                    if (type == ReturnType::SIMPLE)
                    {
                        returns(i, j) = (p_t - p_tm1) / p_tm1;
                    }
                    else
                    { // LOG
                        returns(i, j) = std::log(p_t / p_tm1);
                    }
                }
            }
        }

        return returns;
    }

    Eigen::VectorXd MarketData::calculate_returns(const std::string &ticker,
                                                  ReturnType type) const
    {
        Eigen::VectorXd prices = get_prices(ticker);

        if (prices.size() < 2)
        {
            throw std::runtime_error("Need at least 2 price observations");
        }

        Eigen::VectorXd returns(prices.size() - 1);

        for (int i = 0; i < prices.size() - 1; ++i)
        {
            double p_t = prices(i + 1);
            double p_tm1 = prices(i);

            if (std::isnan(p_t) || std::isnan(p_tm1) || p_tm1 == 0.0)
            {
                returns(i) = std::numeric_limits<double>::quiet_NaN();
            }
            else
            {
                if (type == ReturnType::SIMPLE)
                {
                    returns(i) = (p_t - p_tm1) / p_tm1;
                }
                else
                { // LOG
                    returns(i) = std::log(p_t / p_tm1);
                }
            }
        }

        return returns;
    }

    // ================================
    // Data Filtering and Manipulation
    // ================================

    MarketData MarketData::filter_by_date(const std::string &start_date,
                                          const std::string &end_date) const
    {
        int start_idx = find_date_index(start_date);
        int end_idx = find_date_index(end_date);

        if (start_idx < 0)
        {
            throw std::invalid_argument("Start date not found: " + start_date);
        }
        if (end_idx < 0)
        {
            throw std::invalid_argument("End date not found: " + end_date);
        }
        if (start_idx > end_idx)
        {
            throw std::invalid_argument("Start date must be before end date");
        }

        int num_periods = end_idx - start_idx + 1;
        Eigen::MatrixXd filtered_prices = prices_.block(start_idx, 0, num_periods, prices_.cols());

        std::vector<std::string> filtered_dates(dates_.begin() + start_idx,
                                                dates_.begin() + end_idx + 1);

        return MarketData(filtered_prices, filtered_dates, tickers_);
    }

    MarketData MarketData::select_assets(const std::vector<std::string> &selected_tickers) const
    {
        std::vector<int> indices;
        indices.reserve(selected_tickers.size());

        for (const auto &ticker : selected_tickers)
        {
            int idx = find_ticker_index(ticker);
            if (idx < 0)
            {
                throw std::invalid_argument("Ticker not found: " + ticker);
            }
            indices.push_back(idx);
        }

        Eigen::MatrixXd selected_prices(prices_.rows(), indices.size());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            selected_prices.col(i) = prices_.col(indices[i]);
        }

        return MarketData(selected_prices, dates_, selected_tickers);
    }

    MarketData MarketData::drop_missing() const
    {
        // Find rows without NaN values
        std::vector<int> valid_rows;
        for (int i = 0; i < prices_.rows(); ++i)
        {
            bool has_nan = false;
            for (int j = 0; j < prices_.cols(); ++j)
            {
                if (std::isnan(prices_(i, j)))
                {
                    has_nan = true;
                    break;
                }
            }
            if (!has_nan)
            {
                valid_rows.push_back(i);
            }
        }

        if (valid_rows.empty())
        {
            throw std::runtime_error("All rows contain NaN values");
        }

        Eigen::MatrixXd clean_prices(valid_rows.size(), prices_.cols());
        std::vector<std::string> clean_dates;
        clean_dates.reserve(valid_rows.size());

        for (size_t i = 0; i < valid_rows.size(); ++i)
        {
            clean_prices.row(i) = prices_.row(valid_rows[i]);
            clean_dates.push_back(dates_[valid_rows[i]]);
        }

        return MarketData(clean_prices, clean_dates, tickers_);
    }

    MarketData MarketData::forward_fill() const
    {
        Eigen::MatrixXd filled_prices = prices_;

        for (int j = 0; j < filled_prices.cols(); ++j)
        {
            double last_valid = std::numeric_limits<double>::quiet_NaN();

            for (int i = 0; i < filled_prices.rows(); ++i)
            {
                if (!std::isnan(filled_prices(i, j)))
                {
                    last_valid = filled_prices(i, j);
                }
                else if (!std::isnan(last_valid))
                {
                    filled_prices(i, j) = last_valid;
                }
            }
        }

        return MarketData(filled_prices, dates_, tickers_);
    }

    // ====================
    // Statistical Methods
    // ====================

    Eigen::VectorXd MarketData::mean_returns(ReturnType type) const
    {
        Eigen::MatrixXd returns = calculate_returns(type);

        Eigen::VectorXd means(returns.cols());
        for (int j = 0; j < returns.cols(); ++j)
        {
            double sum = 0.0;
            int count = 0;

            for (int i = 0; i < returns.rows(); ++i)
            {
                if (!std::isnan(returns(i, j)))
                {
                    sum += returns(i, j);
                    ++count;
                }
            }

            means(j) = (count > 0) ? sum / count : std::numeric_limits<double>::quiet_NaN();
        }

        return means;
    }

    Eigen::VectorXd MarketData::volatilities(ReturnType type) const
    {
        Eigen::MatrixXd returns = calculate_returns(type);
        Eigen::VectorXd means = mean_returns(type);

        Eigen::VectorXd vols(returns.cols());
        for (int j = 0; j < returns.cols(); ++j)
        {
            double sum_sq_dev = 0.0;
            int count = 0;

            for (int i = 0; i < returns.rows(); ++i)
            {
                if (!std::isnan(returns(i, j)))
                {
                    double dev = returns(i, j) - means(j);
                    sum_sq_dev += dev * dev;
                    ++count;
                }
            }

            vols(j) = (count > 1) ? std::sqrt(sum_sq_dev / (count - 1)) : std::numeric_limits<double>::quiet_NaN();
        }

        return vols;
    }

    Eigen::MatrixXd MarketData::correlation_matrix(ReturnType type) const
    {
        Eigen::MatrixXd returns = calculate_returns(type);

        // Calculate standardized returns
        Eigen::VectorXd means = mean_returns(type);
        Eigen::VectorXd stds = volatilities(type);

        Eigen::MatrixXd std_returns = returns;
        for (int j = 0; j < returns.cols(); ++j)
        {
            for (int i = 0; i < returns.rows(); ++i)
            {
                if (!std::isnan(returns(i, j)) && stds(j) > 0)
                {
                    std_returns(i, j) = (returns(i, j) - means(j)) / stds(j);
                }
            }
        }

        // Calculate correlation
        Eigen::MatrixXd correlation = (std_returns.transpose() * std_returns) / (returns.rows() - 1);

        return correlation;
    }

    // ===================
    // Validation Methods
    // ===================

    bool MarketData::is_valid() const
    {
        if (prices_.rows() == 0 || prices_.cols() == 0)
        {
            return false;
        }
        if (dates_.size() != static_cast<size_t>(prices_.rows()))
        {
            return false;
        }
        if (tickers_.size() != static_cast<size_t>(prices_.cols()))
        {
            return false;
        }
        return true;
    }

    size_t MarketData::count_missing() const
    {
        size_t count = 0;
        for (int i = 0; i < prices_.rows(); ++i)
        {
            for (int j = 0; j < prices_.cols(); ++j)
            {
                if (std::isnan(prices_(i, j)))
                {
                    ++count;
                }
            }
        }
        return count;
    }

    void MarketData::print_summary() const
    {
        std::cout << "\n=== Market Data Summary ===\n";
        std::cout << "Dimensions: " << prices_.rows() << " dates × "
                  << prices_.cols() << " assets\n";
        std::cout << "Date range: " << dates_.front() << " to " << dates_.back() << "\n";
        std::cout << "Assets: ";
        for (const auto &ticker : tickers_)
        {
            std::cout << ticker << " ";
        }
        std::cout << "\nMissing values: " << count_missing() << "\n";
        std::cout << "==========================\n"
                  << std::endl;
    }

    // =========================
    // Private Helper Methods
    // =========================

    int MarketData::find_ticker_index(const std::string &ticker) const
    {
        auto it = ticker_index_.find(ticker);
        if (it != ticker_index_.end())
        {
            return static_cast<int>(it->second);
        }

        // Fallback to linear search if map not built
        auto vec_it = std::find(tickers_.begin(), tickers_.end(), ticker);
        if (vec_it != tickers_.end())
        {
            return static_cast<int>(std::distance(tickers_.begin(), vec_it));
        }

        return -1;
    }

    int MarketData::find_date_index(const std::string &date) const
    {
        auto it = date_index_.find(date);
        if (it != date_index_.end())
        {
            return static_cast<int>(it->second);
        }

        // Fallback to linear search if map not built
        auto vec_it = std::find(dates_.begin(), dates_.end(), date);
        if (vec_it != dates_.end())
        {
            return static_cast<int>(std::distance(dates_.begin(), vec_it));
        }

        return -1;
    }

    void MarketData::build_index_maps()
    {
        date_index_.clear();
        ticker_index_.clear();

        for (size_t i = 0; i < dates_.size(); ++i)
        {
            date_index_[dates_[i]] = i;
        }

        for (size_t i = 0; i < tickers_.size(); ++i)
        {
            ticker_index_[tickers_[i]] = i;
        }
    }

} // namespace portfolio
