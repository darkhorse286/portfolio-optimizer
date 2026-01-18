/**
 * @file generate_synthetic_data.cpp
 * @brief Generate synthetic market data for portfolio optimizer
 */

#include "data/data_loader.hpp"
#include "data/market_data.hpp"
#include <iostream>
#include <iomanip>

using namespace portfolio;

int main(int argc, char* argv[]) {
    std::cout << "\n=== Synthetic Data Generator ===\n" << std::endl;
    
    // Configuration
    std::vector<std::string> tickers = {
        "AAPL",   // Apple - Tech
        "MSFT",   // Microsoft - Tech
        "JPM",    // JPMorgan - Finance
        "JNJ",    // Johnson & Johnson - Healthcare
        "XOM",    // Exxon Mobil - Energy
        "WMT",    // Walmart - Consumer
        "GOOGL",  // Google - Tech
        "BAC",    // Bank of America - Finance
        "PFE",    // Pfizer - Healthcare
        "CVX"     // Chevron - Energy
    };
    
    // 2 years of trading days (approximately 252 days per year)
    size_t num_days = 504;
    std::string start_date = "2022-01-03";  // First trading day of 2022
    
    // Different volatilities for different sectors
    // Tech stocks: higher volatility
    // Consumer/Healthcare: medium volatility  
    // Finance/Energy: medium-high volatility
    
    std::cout << "Generating data for " << tickers.size() << " assets..." << std::endl;
    std::cout << "Time period: " << start_date << " to ~2 years" << std::endl;
    std::cout << "Trading days: " << num_days << std::endl;
    
    // Parse command line arguments
    std::string output_file = "data/market/historical_prices.csv";
    double base_volatility = 0.015;  // 1.5% daily volatility
    double base_drift = 0.0003;      // ~8% annualized return
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--volatility" && i + 1 < argc) {
            base_volatility = std::stod(argv[++i]);
        } else if (arg == "--drift" && i + 1 < argc) {
            base_drift = std::stod(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n"
                      << "Options:\n"
                      << "  --output FILE      Output CSV file (default: data/market/historical_prices.csv)\n"
                      << "  --volatility VAL   Base daily volatility (default: 0.015)\n"
                      << "  --drift VAL        Base daily drift (default: 0.0003)\n"
                      << "  --help             Show this help\n";
            return 0;
        }
    }
    
    // Generate data with sector-specific characteristics
    std::cout << "\nGenerating prices..." << std::endl;
    auto data = DataLoader::generate_synthetic_data(
        tickers,
        num_days,
        start_date,
        base_volatility,
        base_drift
    );
    
    // Save to CSV
    std::cout << "Saving to " << output_file << "..." << std::endl;
    DataLoader::save_csv_wide(data, output_file);
    
    // Print summary statistics
    std::cout << "\n=== Generated Data Summary ===\n";
    std::cout << "Dates: " << data.num_dates() << " (" 
              << data.get_dates().front() << " to " 
              << data.get_dates().back() << ")\n";
    std::cout << "Assets: " << data.num_assets() << "\n";
    
    auto mean_returns = data.mean_returns(ReturnType::SIMPLE);
    auto volatilities = data.volatilities(ReturnType::SIMPLE);
    
    std::cout << "\nAsset Statistics (Annualized):\n";
    std::cout << std::string(60, '-') << "\n";
    std::cout << std::setw(8) << "Ticker" 
              << std::setw(15) << "Mean Return" 
              << std::setw(15) << "Volatility"
              << std::setw(15) << "Sharpe (0.02)\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (size_t i = 0; i < data.num_assets(); ++i) {
        double ann_return = mean_returns(i) * 252;
        double ann_vol = volatilities(i) * std::sqrt(252);
        double sharpe = (ann_return - 0.02) / ann_vol;
        
        std::cout << std::setw(8) << data.get_tickers()[i]
                  << std::setw(14) << std::fixed << std::setprecision(2) 
                  << (ann_return * 100) << "%"
                  << std::setw(14) << (ann_vol * 100) << "%"
                  << std::setw(15) << std::setprecision(3) << sharpe << "\n";
    }
    std::cout << std::string(60, '-') << "\n";
    
    // Calculate and show correlation matrix
    auto corr = data.correlation_matrix(ReturnType::SIMPLE);
    double avg_corr = 0.0;
    int count = 0;
    for (int i = 0; i < corr.rows(); ++i) {
        for (int j = i + 1; j < corr.cols(); ++j) {
            avg_corr += corr(i, j);
            ++count;
        }
    }
    avg_corr /= count;
    
    std::cout << "\nAverage pairwise correlation: " 
              << std::fixed << std::setprecision(3) << avg_corr << "\n";
    
    std::cout << "\n✓ Data generation complete!\n" << std::endl;
    std::cout << "You can now run:\n";
    std::cout << "  ./build/bin/portfolio_optimizer --config data/config/portfolio_config.json --verbose\n";
    std::cout << std::endl;
    
    return 0;
}