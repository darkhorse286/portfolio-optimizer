/**
 * @file data_loader.cpp
 * @brief Implementation of DataLoader class and configuration structures
 */

#include "data/data_loader.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <iomanip>
#include <set>
#include <map>

namespace portfolio
{

    // =============================================
    // Configuration Structures - from_json Methods
    // =============================================

    DataConfig DataConfig::from_json(const nlohmann::json &j)
    {
        DataConfig config;
        config.start_date = j.value("start_date", "");
        config.end_date = j.value("end_date", "");
        config.universe = j.value("universe", std::vector<std::string>{});
        config.benchmark = j.value("benchmark", "SPY");
        config.data_file = j.value("data_file", "data/market/historical_prices.csv");
        config.frequency = j.value("frequency", "daily");
        return config;
    }

    OptimizerConfig OptimizerConfig::from_json(const nlohmann::json &j)
    {
        OptimizerConfig config;
        config.type = j.value("type", "mean_variance");
        config.objective = j.value("objective", "max_sharpe");
        config.risk_aversion = j.value("risk_aversion", 1.0);
        config.risk_free_rate = j.value("risk_free_rate", 0.02);

        if (j.contains("constraints"))
        {
            const auto &constraints = j["constraints"];
            config.max_weight = constraints.value("max_weight", 0.30);
            config.min_weight = constraints.value("min_weight", 0.05);
            config.sum_to_one = constraints.value("sum_to_one", true);
            config.long_only = constraints.value("long_only", true);
            config.max_turnover = constraints.value("max_turnover", 0.50);
        }
        else
        {
            config.max_weight = 0.30;
            config.min_weight = 0.05;
            config.sum_to_one = true;
            config.long_only = true;
            config.max_turnover = 0.50;
        }

        return config;
    }

    RiskModelConfig RiskModelConfig::from_json(const nlohmann::json &j)
    {
        RiskModelConfig config;
        config.type = j.value("type", "sample_covariance");
        config.estimation_window = j.value("estimation_window", 252);
        config.min_periods = j.value("min_periods", 60);
        config.method = j.value("method", "sample");
        config.ewma_lambda = j.value("ewma_lambda", 0.94);
        config.shrinkage_enabled = false;

        if (j.contains("shrinkage"))
        {
            config.shrinkage_enabled = j["shrinkage"].value("enabled", false);
        }

        return config;
    }

    BacktestConfig BacktestConfig::from_json(const nlohmann::json &j)
    {
        BacktestConfig config;
        config.initial_capital = j.value("initial_capital", 1000000.0);
        config.rebalance_frequency = j.value("rebalance_frequency", "monthly");
        config.lookback_window = j.value("lookback_window", 252);
        config.min_history = j.value("min_history", 60);

        if (j.contains("transaction_costs"))
        {
            const auto &tc = j["transaction_costs"];
            config.commission_rate = tc.value("commission_rate", 0.001);
            config.slippage_bps = tc.value("slippage_bps", 5.0);
        }
        else
        {
            config.commission_rate = 0.001;
            config.slippage_bps = 5.0;
        }

        if (j.contains("cash_management"))
        {
            const auto &cm = j["cash_management"];
            config.allow_cash = cm.value("allow_cash", true);
            config.cash_return = cm.value("cash_return", 0.0);
        }
        else
        {
            config.allow_cash = true;
            config.cash_return = 0.0;
        }

        return config;
    }

    PortfolioConfig PortfolioConfig::load_from_file(const std::string &config_path)
    {
        return DataLoader::load_config(config_path);
    }

    // ===========================
    // CSV Loading - Wide Format
    // ===========================

    MarketData DataLoader::load_csv_wide(const std::string &filepath,
                                         const std::vector<std::string> &tickers)
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filepath);
        }

        std::string line;
        std::vector<std::string> all_tickers;
        std::vector<std::string> dates;
        std::vector<std::vector<double>> price_data;

        // Read header line
        if (!std::getline(file, line))
        {
            throw std::runtime_error("Empty CSV file");
        }

        auto header = parse_csv_line(line);
        if (header.empty() || header[0] != "date")
        {
            throw std::runtime_error("CSV must start with 'date' column");
        }

        // Extract tickers from header
        for (size_t i = 1; i < header.size(); ++i)
        {
            all_tickers.push_back(trim(header[i]));
        }

        // Determine which columns to load
        std::vector<size_t> column_indices;
        std::vector<std::string> selected_tickers;

        if (tickers.empty())
        {
            // Load all tickers
            for (size_t i = 0; i < all_tickers.size(); ++i)
            {
                column_indices.push_back(i);
                selected_tickers.push_back(all_tickers[i]);
            }
        }
        else
        {
            // Load only specified tickers
            for (const auto &ticker : tickers)
            {
                auto it = std::find(all_tickers.begin(), all_tickers.end(), ticker);
                if (it != all_tickers.end())
                {
                    size_t idx = std::distance(all_tickers.begin(), it);
                    column_indices.push_back(idx);
                    selected_tickers.push_back(ticker);
                }
            }

            if (column_indices.empty())
            {
                throw std::runtime_error("None of the specified tickers found in CSV");
            }
        }

        // Read data rows
        while (std::getline(file, line))
        {
            if (line.empty())
                continue;

            auto fields = parse_csv_line(line);
            if (fields.size() < 2)
                continue;

            std::string date = trim(fields[0]);
            if (!is_valid_date_format(date))
            {
                continue; // Skip invalid dates
            }

            dates.push_back(date);

            std::vector<double> row_prices;
            row_prices.reserve(column_indices.size());

            for (size_t idx : column_indices)
            {
                if (idx + 1 < fields.size())
                {
                    row_prices.push_back(safe_stod(fields[idx + 1]));
                }
                else
                {
                    row_prices.push_back(std::numeric_limits<double>::quiet_NaN());
                }
            }

            price_data.push_back(row_prices);
        }

        file.close();

        if (dates.empty())
        {
            throw std::runtime_error("No valid data found in CSV file");
        }

        // Convert to Eigen matrix
        Eigen::MatrixXd prices(dates.size(), selected_tickers.size());
        for (size_t i = 0; i < dates.size(); ++i)
        {
            for (size_t j = 0; j < selected_tickers.size(); ++j)
            {
                prices(i, j) = price_data[i][j];
            }
        }

        return MarketData(prices, dates, selected_tickers);
    }

    // ===========================
    // CSV Loading - Long Format
    // ===========================

    MarketData DataLoader::load_csv_long(const std::string &filepath,
                                         const std::vector<std::string> &tickers)
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filepath);
        }

        std::string line;
        std::map<std::string, std::map<std::string, double>> data_map; // date -> ticker -> price
        std::set<std::string> all_dates;
        std::set<std::string> all_tickers;

        // Skip header
        std::getline(file, line);

        // Read data
        while (std::getline(file, line))
        {
            if (line.empty())
                continue;

            auto fields = parse_csv_line(line);
            if (fields.size() < 3)
                continue;

            std::string date = trim(fields[0]);
            std::string ticker = trim(fields[1]);
            double price = safe_stod(fields[2]);

            if (!is_valid_date_format(date))
                continue;

            // Filter by tickers if specified
            if (!tickers.empty() &&
                std::find(tickers.begin(), tickers.end(), ticker) == tickers.end())
            {
                continue;
            }

            data_map[date][ticker] = price;
            all_dates.insert(date);
            all_tickers.insert(ticker);
        }

        file.close();

        if (all_dates.empty() || all_tickers.empty())
        {
            throw std::runtime_error("No valid data found in CSV file");
        }

        // Convert sets to vectors
        std::vector<std::string> dates(all_dates.begin(), all_dates.end());
        std::vector<std::string> ticker_vec(all_tickers.begin(), all_tickers.end());

        // Build matrix
        Eigen::MatrixXd prices(dates.size(), ticker_vec.size());
        prices.setConstant(std::numeric_limits<double>::quiet_NaN());

        for (size_t i = 0; i < dates.size(); ++i)
        {
            const auto &date = dates[i];
            if (data_map.count(date))
            {
                for (size_t j = 0; j < ticker_vec.size(); ++j)
                {
                    const auto &ticker = ticker_vec[j];
                    if (data_map[date].count(ticker))
                    {
                        prices(i, j) = data_map[date][ticker];
                    }
                }
            }
        }

        return MarketData(prices, dates, ticker_vec);
    }

    // ========================
    // Auto-detect CSV Format
    // ========================

    MarketData DataLoader::load_csv(const std::string &filepath,
                                    const std::vector<std::string> &tickers)
    {
        // Try to detect format by reading first line
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file: " + filepath);
        }

        std::string line;
        std::getline(file, line);
        file.close();

        auto header = parse_csv_line(line);

        // Wide format: date, ticker1, ticker2, ...
        // Long format: date, ticker, price
        if (header.size() == 3 &&
            (header[1] == "ticker" || header[1] == "symbol"))
        {
            return load_csv_long(filepath, tickers);
        }
        else
        {
            return load_csv_wide(filepath, tickers);
        }
    }

    // ================
    // JSON Loading
    // ================

    nlohmann::json DataLoader::load_json(const std::string &filepath)
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open JSON file: " + filepath);
        }

        nlohmann::json j;
        try
        {
            file >> j;
        }
        catch (const nlohmann::json::exception &e)
        {
            throw std::runtime_error("JSON parsing error: " + std::string(e.what()));
        }

        file.close();
        return j;
    }

    PortfolioConfig DataLoader::load_config(const std::string &config_path)
    {
        auto j = load_json(config_path);

        PortfolioConfig config;

        if (j.contains("data"))
        {
            config.data = DataConfig::from_json(j["data"]);
        }

        if (j.contains("optimizer"))
        {
            config.optimizer = OptimizerConfig::from_json(j["optimizer"]);
        }

        if (j.contains("risk_model"))
        {
            config.risk_model = risk::RiskModelConfig::from_json(j["risk_model"]); // Use risk:: version
        }

        if (j.contains("backtest"))
        {
            config.backtest = BacktestConfig::from_json(j["backtest"]);
        }

        return config;
    }

    // ===========================
    // Synthetic Data Generation
    // ===========================

    MarketData DataLoader::generate_synthetic_data(
        const std::vector<std::string> &tickers,
        size_t num_days,
        const std::string &start_date,
        double volatility,
        double drift)
    {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> dist(drift, volatility);

        Eigen::MatrixXd prices(num_days, tickers.size());
        std::vector<std::string> dates;
        dates.reserve(num_days);

        // Generate dates
        for (size_t i = 0; i < num_days; ++i)
        {
            dates.push_back(add_days(start_date, i));
        }

        // Generate prices (geometric Brownian motion)
        for (size_t j = 0; j < tickers.size(); ++j)
        {
            prices(0, j) = 100.0; // Initial price

            for (size_t i = 1; i < num_days; ++i)
            {
                double return_val = dist(gen);
                prices(i, j) = prices(i - 1, j) * (1.0 + return_val);
            }
        }

        return MarketData(prices, dates, tickers);
    }

    // ==================
    // Export Methods
    // ==================

    void DataLoader::save_csv_wide(const MarketData &data, const std::string &filepath)
    {
        std::ofstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }

        // Write header
        file << "date";
        for (const auto &ticker : data.get_tickers())
        {
            file << "," << ticker;
        }
        file << "\n";

        // Write data
        const auto &prices = data.get_Prices();
        const auto &dates = data.get_dates();

        for (size_t i = 0; i < dates.size(); ++i)
        {
            file << dates[i];
            for (int j = 0; j < prices.cols(); ++j)
            {
                file << "," << std::fixed << std::setprecision(6) << prices(i, j);
            }
            file << "\n";
        }

        file.close();
    }

    void DataLoader::save_csv_long(const MarketData &data, const std::string &filepath)
    {
        std::ofstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }

        // Write header
        file << "date,ticker,price\n";

        // Write data
        const auto &prices = data.get_Prices();
        const auto &dates = data.get_dates();
        const auto &tickers = data.get_tickers();

        for (size_t i = 0; i < dates.size(); ++i)
        {
            for (size_t j = 0; j < tickers.size(); ++j)
            {
                file << dates[i] << "," << tickers[j] << ","
                     << std::fixed << std::setprecision(6) << prices(i, j) << "\n";
            }
        }

        file.close();
    }

    // =======================
    // Private Helper Methods
    // =======================

    std::vector<std::string> DataLoader::parse_csv_line(const std::string &line)
    {
        std::vector<std::string> tokens;
        std::string token;
        bool in_quotes = false;

        for (char c : line)
        {
            if (c == '"')
            {
                in_quotes = !in_quotes;
            }
            else if (c == ',' && !in_quotes)
            {
                tokens.push_back(token);
                token.clear();
            }
            else
            {
                token += c;
            }
        }

        tokens.push_back(token);
        return tokens;
    }

    bool DataLoader::is_valid_date_format(const std::string &date)
    {
        // Simple check for YYYY-MM-DD format
        if (date.length() != 10)
            return false;
        if (date[4] != '-' || date[7] != '-')
            return false;

        for (size_t i = 0; i < date.length(); ++i)
        {
            if (i == 4 || i == 7)
                continue;
            if (!std::isdigit(date[i]))
                return false;
        }

        return true;
    }

    std::string DataLoader::trim(const std::string &str)
    {
        size_t first = str.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return "";

        size_t last = str.find_last_not_of(" \t\r\n");
        return str.substr(first, last - first + 1);
    }

    double DataLoader::safe_stod(const std::string &str)
    {
        try
        {
            std::string trimmed = trim(str);
            if (trimmed.empty() || trimmed == "nan" || trimmed == "NaN")
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return std::stod(trimmed);
        }
        catch (...)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::string DataLoader::add_days(const std::string &start_date, int days_offset)
    {
        // Simple implementation - just for synthetic data generation
        // In production, use proper date library

        struct tm tm = {};
        std::istringstream ss(start_date);
        ss >> std::get_time(&tm, "%Y-%m-%d");

        std::time_t time = std::mktime(&tm);
        time += days_offset * 24 * 3600;

        std::tm *ptm = std::localtime(&time);
        char buffer[11];
        std::strftime(buffer, 11, "%Y-%m-%d", ptm);

        return std::string(buffer);
    }

} // namespace portfolio
