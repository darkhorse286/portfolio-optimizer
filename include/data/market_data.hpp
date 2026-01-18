/*
 * @file market_data.hpp
 * @brief Time-series market data storage and manipulation.
 *
 * Provides efficient storage and access to market data using Eigen matrices.
 * Supports multiple assets, data indexing, and return calculations.
 */

#ifndef MARKET_DATA_HPP
#define MARKET_DATA_HPP

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>

namespace portfolio
{
    /**
     *  @enum ReturnType
     *  @brief Type of return calculation.
     */
    enum class ReturnType
    {
        SIMPLE, /**< Simple returns (P_t / P_{t-1} - 1) */
        LOG     /**< Logarithmic returns (log(P_t / P_{t-1})) */
    };

    /**
     * @class MarketData
     * @brief Container for multi-asset time-series price and return data.
     *
     * Stores historical price data as Eigen matrices with associated date indices.
     * Provides methods for return calculations, date alignment, and basic statistics.
     *
     * @note All data is stored in column-major format (time x assets).
     * @note Missing data is represented as NaN values.
     */

    class MarketData
    {
    public:
        /**
         * @brief Constructor with dimensions.
         * @param num_dates Number of time periods.
         * @param num_assets Number of assets.
         */
        MarketData(size_t num_dates, size_t num_assets);

        /**
         * @brief Constructor with data.
         * @param prices Price matrix (dates x assets).
         * @param dates Vector of date strings.
         * @param tickers Vector of asset ticker symbols.
         */
        MarketData(const Eigen::MatrixXd &prices,
                   const std::vector<std::string> &dates,
                   const std::vector<std::string> &tickers);
        
        // Destructor
        ~MarketData() = default;

        /** ===========================================
         *  Data Access Methods
         *  ===========================================
         */

        /**
         * @brief Get the full price matrix.
         * @return Reference to price matrix (dates x assets).
         */
        const Eigen::MatrixXd &get_Prices() const
        {
            return prices_;
        }

        /**
         * @brief Get prices for a specific asset.
         * @param ticker Asset ticker symbol.
         * @return Price vector for the asset.
         */
        Eigen::VectorXd get_prices(const std::string &ticker) const;

        /**
         * @brief Get prices for a specific date.
         * @param date Date string (YYYY-MM-DD).
         * @return Price vector for all assets on that date.
         */
        Eigen::VectorXd get_prices_at_date(const std::string &date) const;

        /**
         * @brief Get price for specific asset and date.
         * @param ticker Asset ticker symbol.
         * @param date Date string (YYYY-MM-DD).
         * @return Price value.
         */
        double get_price(const std::string &ticker, const std::string &date) const;

        /**
         * @brief Get vector of dates.
         */
        const std::vector<std::string> &get_dates() const
        {
            return dates_;
        }

        /**
         * @brief Get vector of tickers.
         */
        const std::vector<std::string> &get_tickers() const
        {
            return tickers_;
        }

        /**
         * @brief Get number of dates (time periods).
         */
        size_t num_dates() const
        {
            return prices_.rows();
        }

        /**
         * @brief Get number of assets.
         */
        size_t num_assets() const
        {
            return prices_.cols();
        }

        /** ===========================================
         *  Data Modification Methods
         *  ===========================================
         */

        /**
         * @brief Set the price matrix
         * @param prices New price matrix
         * @throws std::invalid_argument if dimensions don't match
         */
        void set_prices(const Eigen::MatrixXd &prices);

        /**
         * @brief Set price for specific asset and date
         * @param ticker Asset ticker
         * @param date Date string
         * @param price Price value
         */
        void set_price(const std::string &ticker, const std::string &date, double price);

        /**
         * @brief Set dates vector
         * @param dates Vector of date strings
         * @throws std::invalid_argument if size doesn't match price matrix
         */
        void set_dates(const std::vector<std::string> &dates);

        /**
         * @brief Set tickers vector
         * @param tickers Vector of ticker strings
         * @throws std::invalid_argument if size doesn't match price matrix
         */
        void set_tickers(const std::vector<std::string> &tickers);

        /** ===========================================
         *  Data Calculation Methods
         *  ===========================================
         */

        /**
         * @brief Calculate returns from prices
         * @param type Type of return calculation (SIMPLE or LOG)
         * @return Matrix of returns (dates-1 x assets)
         */
        Eigen::MatrixXd calculate_returns(ReturnType type = ReturnType::SIMPLE) const;

        /**
         * @brief Calculate returns for specific asset
         * @param ticker Asset ticker
         * @param type Type of return calculation
         * @return Vector of returns
         */
        Eigen::VectorXd calculate_returns(const std::string &ticker,
                                          ReturnType type = ReturnType::SIMPLE) const;

        /** ===========================================
         *  Data Filtering and Manipulation Methods
         *  ===========================================
         */

        /**
         * @brief Filter data by date range
         * @param start_date Start date (inclusive)
         * @param end_date End date (inclusive)
         * @return New MarketData object with filtered data
         */
        MarketData filter_by_date(const std::string &start_date,
                                  const std::string &end_date) const;

        /**
         * @brief Select subset of assets
         * @param selected_tickers Vector of tickers to keep
         * @return New MarketData object with selected assets
         */
        MarketData select_assets(const std::vector<std::string> &selected_tickers) const;

        /**
         * @brief Remove rows with missing data
         * @return New MarketData object without NaN values
         */
        MarketData drop_missing() const;

        /**
         * @brief Forward-fill missing data
         * @return New MarketData object with filled values
         */
        MarketData forward_fill() const;

        /** ===========================================
         *  Statistical Methods
         *  ===========================================
         */

        /**
         * @brief Calculate mean returns
         * @param type Type of return calculation
         * @return Vector of mean returns for each asset
         */
        Eigen::VectorXd mean_returns(ReturnType type = ReturnType::SIMPLE) const;

        /**
         * @brief Calculate return volatilities (std dev)
         * @param type Type of return calculation
         * @return Vector of volatilities for each asset
         */
        Eigen::VectorXd volatilities(ReturnType type = ReturnType::SIMPLE) const;

        /**
         * @brief Calculate correlation matrix
         * @param type Type of return calculation
         * @return Correlation matrix (assets x assets)
         */
        Eigen::MatrixXd correlation_matrix(ReturnType type = ReturnType::SIMPLE) const;

        /** ===========================================
         *  Validation Methods
         *  ===========================================
         */
        /**
         * @brief Check if data is valid and consistent
         * @return true if valid, false otherwise
         */
        bool is_valid() const;

        /**
         * @brief Count missing values
         * @return Number of NaN entries in price matrix
         */
        size_t count_missing() const;

        /**
         * @brief Print summary statistics
         */
        void print_summary() const;

    private:
        /** ===========================================
         *  Private Helper Methods
         *  ===========================================
         */

        /**
         * @brief Find index of ticker in tickers_ vector
         * @param ticker Asset ticker
         * @return Index, or -1 if not found
         */
        int find_ticker_index(const std::string &ticker) const;

        /**
         * @brief Find index of date in dates_ vector
         * @param date Date string
         * @return Index, or -1 if not found
         */
        int find_date_index(const std::string &date) const;

        /**
         * @brief Build index maps for fast lookups
         */
        void build_index_maps();

        /** ===========================================
         *  Private Member Variables
         *  ===========================================
         */
        Eigen::MatrixXd prices_;                     ///< Price matrix (dates x assets)
        std::vector<std::string> dates_;             ///< Date strings
        std::vector<std::string> tickers_;           ///< Asset tickers
        std::map<std::string, size_t> date_index_;   ///< Date to index map
        std::map<std::string, size_t> ticker_index_; ///< Ticker to index map
    };

}
#endif // MARKET_DATA_HPP