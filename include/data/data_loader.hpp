/**
 * @file data_loader.hpp
 * @brief Data loading and parsing utilities
 * 
 * Provides functionality to load market data from CSV files and
 * configuration from JSON files.
 */

#ifndef DATA_LOADER_HPP
#define DATA_LOADER_HPP

#include "market_data.hpp"
#include "risk/risk_model_factory.hpp" 
#include <nlohmann/json.hpp>
#include <string>
#include <memory>
#include <optional>


namespace portfolio {

/**
 * @struct DataConfig
 * @brief Configuration parameters for data loading
 */
struct DataConfig {
    std::string start_date;                    ///< Start date filter
    std::string end_date;                      ///< End date filter
    std::vector<std::string> universe;         ///< List of tickers to load
    std::string benchmark;                     ///< Benchmark ticker
    std::string data_file;                     ///< Path to data file
    std::string frequency;                     ///< Data frequency (daily, weekly, etc.)
    
    /**
     * @brief Load from JSON object
     */
    static DataConfig from_json(const nlohmann::json& j);
};

/**
 * @struct OptimizerConfig
 * @brief Configuration for optimizer parameters
 */
struct OptimizerConfig {
    std::string type;                          ///< Optimizer type (mean_variance, qubo)
    std::string objective;                     ///< Objective function
    double risk_aversion;                      ///< Risk aversion parameter
    double max_weight;                         ///< Maximum asset weight
    double min_weight;                         ///< Minimum asset weight
    bool sum_to_one;                           ///< Constraint: weights sum to 1
    bool long_only;                            ///< Constraint: no short positions
    double max_turnover;                       ///< Maximum portfolio turnover
    double risk_free_rate;                     ///< Risk-free rate
    
    static OptimizerConfig from_json(const nlohmann::json& j);
};

/**
 * @struct RiskModelConfig
 * @brief Configuration for risk model estimation
 */
struct RiskModelConfig {
    std::string type;                          ///< Type (sample_covariance, ewma, ledoit_wolf)
    int estimation_window;                     ///< Rolling window size
    int min_periods;                           ///< Minimum required periods
    std::string method;                        ///< Estimation method
    double ewma_lambda;                        ///< EWMA decay parameter
    bool shrinkage_enabled;                    ///< Use shrinkage estimation
    
    static RiskModelConfig from_json(const nlohmann::json& j);
};

/**
 * @struct BacktestConfig
 * @brief Configuration for backtesting parameters
 */
struct BacktestConfig {
    double initial_capital;                    ///< Starting capital
    std::string rebalance_frequency;           ///< Rebalancing period
    int lookback_window;                       ///< Lookback for estimation
    int min_history;                           ///< Minimum history required
    double commission_rate;                    ///< Transaction cost rate
    double slippage_bps;                       ///< Slippage in basis points
    bool allow_cash;                           ///< Allow cash holdings
    double cash_return;                        ///< Return on cash
    
    static BacktestConfig from_json(const nlohmann::json& j);
};

/**
 * @struct PortfolioConfig
 * @brief Complete portfolio configuration
 */
struct PortfolioConfig {
    DataConfig data;
    OptimizerConfig optimizer;
    risk::RiskModelConfig risk_model;
    BacktestConfig backtest;
    
    /**
     * @brief Load complete configuration from JSON file
     */
    static PortfolioConfig load_from_file(const std::string& config_path);
};

/**
 * @class DataLoader
 * @brief Loads and parses market data from various sources
 * 
 * Supports CSV files with standard formats:
 * - Format 1: date, ticker1, ticker2, ... (wide format)
 * - Format 2: date, ticker, price (long format)
 */
class DataLoader {
public:
    /**
     * @brief Constructor
     */
    DataLoader() = default;
    
    /**
     * @brief Destructor
     */
    ~DataLoader() = default;
    
    // ========================================================================
    // CSV Loading Methods
    // ========================================================================
    
    /**
     * @brief Load market data from CSV file (wide format)
     * 
     * Expected format:
     * date,AAPL,MSFT,JPM,...
     * 2020-01-01,150.0,200.0,120.0,...
     * 
     * @param filepath Path to CSV file
     * @param tickers Optional list of tickers to load (loads all if empty)
     * @return MarketData object
     * @throws std::runtime_error if file cannot be loaded
     */
    static MarketData load_csv_wide(const std::string& filepath,
                                    const std::vector<std::string>& tickers = {});
    
    /**
     * @brief Load market data from CSV file (long format)
     * 
     * Expected format:
     * date,ticker,price
     * 2020-01-01,AAPL,150.0
     * 2020-01-01,MSFT,200.0
     * 
     * @param filepath Path to CSV file
     * @param tickers Optional list of tickers to load
     * @return MarketData object
     * @throws std::runtime_error if file cannot be loaded
     */
    static MarketData load_csv_long(const std::string& filepath,
                                    const std::vector<std::string>& tickers = {});
    
    /**
     * @brief Auto-detect CSV format and load
     * @param filepath Path to CSV file
     * @param tickers Optional list of tickers
     * @return MarketData object
     */
    static MarketData load_csv(const std::string& filepath,
                               const std::vector<std::string>& tickers = {});
    
    // ========================================================================
    // Configuration Loading
    // ========================================================================
    
    /**
     * @brief Load JSON configuration file
     * @param filepath Path to JSON config file
     * @return JSON object
     * @throws std::runtime_error if file cannot be loaded
     */
    static nlohmann::json load_json(const std::string& filepath);
    
    /**
     * @brief Load complete portfolio configuration
     * @param config_path Path to config JSON file
     * @return PortfolioConfig struct
     */
    static PortfolioConfig load_config(const std::string& config_path);
    
    // ========================================================================
    // Data Generation (for testing)
    // ========================================================================
    
    /**
     * @brief Generate synthetic market data for testing
     * @param tickers List of ticker symbols
     * @param num_days Number of trading days
     * @param start_date Starting date
     * @param volatility Daily volatility (default 0.02)
     * @param drift Daily drift (default 0.0005)
     * @return MarketData object with synthetic prices
     */
    static MarketData generate_synthetic_data(
        const std::vector<std::string>& tickers,
        size_t num_days,
        const std::string& start_date = "2020-01-01",
        double volatility = 0.02,
        double drift = 0.0005
    );
    
    // ========================================================================
    // Export Methods
    // ========================================================================
    
    /**
     * @brief Save market data to CSV (wide format)
     * @param data MarketData object
     * @param filepath Output file path
     */
    static void save_csv_wide(const MarketData& data, const std::string& filepath);
    
    /**
     * @brief Save market data to CSV (long format)
     * @param data MarketData object
     * @param filepath Output file path
     */
    static void save_csv_long(const MarketData& data, const std::string& filepath);

private:
    // ========================
    // Private Helper Methods
    // ========================
    
    /**
     * @brief Parse CSV line into tokens
     * @param line CSV line string
     * @return Vector of tokens
     */
    static std::vector<std::string> parse_csv_line(const std::string& line);
    
    /**
     * @brief Validate date format (YYYY-MM-DD)
     * @param date Date string
     * @return true if valid
     */
    static bool is_valid_date_format(const std::string& date);
    
    /**
     * @brief Trim whitespace from string
     * @param str Input string
     * @return Trimmed string
     */
    static std::string trim(const std::string& str);
    
    /**
     * @brief Convert string to double safely
     * @param str String representation of number
     * @return Double value, or NaN if conversion fails
     */
    static double safe_stod(const std::string& str);
    
    /**
     * @brief Generate date string for day offset
     * @param start_date Starting date
     * @param days_offset Number of days to add
     * @return Date string
     */
    static std::string add_days(const std::string& start_date, int days_offset);
};

} // namespace portfolio

#endif // DATA_LOADER_HPP
