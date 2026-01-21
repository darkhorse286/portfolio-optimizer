/**
 * @file main.cpp
 * @brief Main entry point for Quantum Portfolio Optimizer
 * 
 * Command-line application that loads configuration, runs optimization,
 * performs backtesting, and generates results.
 */

#include "data/data_loader.hpp"
#include "data/market_data.hpp"
#include <iostream>
#include <string>
#include <exception>
#include <iomanip>
#include <chrono>

using namespace portfolio;

/**
 * @brief Print usage information
 */
void print_usage(const char* program_name) {
    std::cout << "Quantum Portfolio Optimizer v1.0.0\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  --config PATH         Path to configuration JSON file (required)\n"
              << "  --output PATH         Path to output directory (default: results/)\n"
              << "  --verbose             Enable verbose logging\n"
              << "  --help, -h            Show this help message\n"
              << "\nExample:\n"
              << "  " << program_name << " --config data/config/portfolio_config.json\n"
              << std::endl;
}

/**
 * @brief Print banner
 */
void print_banner() {
    std::cout << "\n"
              << "╔══════════════════════════════════════════════════════════════╗\n"
              << "║       Quantum Portfolio Optimizer v1.0.0                     ║\n"
              << "║       Modern C++17 Portfolio Optimization                    ║\n"
              << "╚══════════════════════════════════════════════════════════════╝\n"
              << std::endl;
}

/**
 * @brief Parse command-line arguments
 */
struct CommandLineArgs {
    std::string config_path;
    std::string output_dir = "results";
    bool verbose = false;
    bool show_help = false;
    
    static CommandLineArgs parse(int argc, char* argv[]) {
        CommandLineArgs args;
        
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            
            if (arg == "--help" || arg == "-h") {
                args.show_help = true;
            }
            else if (arg == "--config" && i + 1 < argc) {
                args.config_path = argv[++i];
            }
            else if (arg == "--output" && i + 1 < argc) {
                args.output_dir = argv[++i];
            }
            else if (arg == "--verbose") {
                args.verbose = true;
            }
            else {
                std::cerr << "Warning: Unknown argument: " << arg << std::endl;
            }
        }
        
        return args;
    }
    
    bool is_valid() const {
        return !show_help && !config_path.empty();
    }
};

/**
 * @brief Main execution function
 */
int run(const CommandLineArgs& args) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // ====================================================================
        // 1. Load Configuration
        // ====================================================================
        std::cout << "[1/5] Loading configuration..." << std::endl;
        
        auto config = DataLoader::load_config(args.config_path);
        
        if (args.verbose) {
            std::cout << "  - Tickers: ";
            for (const auto& ticker : config.data.universe) {
                std::cout << ticker << " ";
            }
            std::cout << "\n  - Date range: " << config.data.start_date 
                      << " to " << config.data.end_date << std::endl;
        }
        
        // ====================================================================
        // 2. Load Market Data
        // ====================================================================
        std::cout << "[2/5] Loading market data..." << std::endl;
        
        auto data = DataLoader::load_csv(config.data.data_file, config.data.universe);
        
        // Filter by date range
        if (!config.data.start_date.empty() && !config.data.end_date.empty()) {
            data = data.filter_by_date(config.data.start_date, config.data.end_date);
        }
        
        if (args.verbose) {
            data.print_summary();
        }
        
        std::cout << "  - Loaded " << data.num_dates() << " dates, "
                  << data.num_assets() << " assets" << std::endl;
        
        // ====================================================================
        // 3. Calculate Returns
        // ====================================================================
        std::cout << "[3/5] Calculating returns..." << std::endl;
        
        auto returns = data.calculate_returns(ReturnType::SIMPLE);
        auto mean_returns = data.mean_returns(ReturnType::SIMPLE);
        
        if (args.verbose) {
            std::cout << "  - Mean returns:\n" << mean_returns.transpose() << std::endl;
        }
        
        // ====================================================================
        // 4. Risk Model Estimation
        // ====================================================================
        std::cout << "[4/5] Estimating risk model..." << std::endl;
        
        // Calculate sample covariance (placeholder for now)
        Eigen::MatrixXd covariance = (returns.transpose() * returns) / (returns.rows() - 1);
        
        if (args.verbose) {
            std::cout << "  - Covariance matrix size: " << covariance.rows() 
                      << "x" << covariance.cols() << std::endl;
        }
        
        // ====================================================================
        // 5. Portfolio Optimization (Placeholder)
        // ====================================================================
        std::cout << "[5/5] Running optimization..." << std::endl;
        
        // Equal-weight portfolio as placeholder
        Eigen::VectorXd weights = Eigen::VectorXd::Constant(data.num_assets(), 
                                                            1.0 / data.num_assets());
        
        std::cout << "\n  Optimal Weights:\n";
        std::cout << "  " << std::string(50, '-') << "\n";
        for (size_t i = 0; i < data.num_assets(); ++i) {
            std::cout << "  " << std::setw(8) << std::left << data.get_tickers()[i]
                      << ": " << std::fixed << std::setprecision(4) 
                      << weights(i) * 100.0 << "%\n";
        }
        std::cout << "  " << std::string(50, '-') << "\n";
        
        // Calculate portfolio statistics
        double portfolio_return = weights.dot(mean_returns) * 252; // Annualized
        double portfolio_variance = weights.transpose() * covariance * weights;
        double portfolio_vol = std::sqrt(portfolio_variance * 252); // Annualized
        double sharpe_ratio = (portfolio_return - config.optimizer.risk_free_rate) / portfolio_vol;
        
        std::cout << "\n  Portfolio Statistics:\n";
        std::cout << "  " << std::string(50, '-') << "\n";
        std::cout << "  Expected Return:  " << std::fixed << std::setprecision(2) 
                  << portfolio_return * 100 << "%\n";
        std::cout << "  Volatility:       " << portfolio_vol * 100 << "%\n";
        std::cout << "  Sharpe Ratio:     " << std::setprecision(3) << sharpe_ratio << "\n";
        std::cout << "  " << std::string(50, '-') << "\n";
        
        // ====================================================================
        // Timing
        // ====================================================================
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time
        ).count();
        
        std::cout << "\n✓ Optimization completed successfully in " 
                  << duration << " ms\n" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Error: " << e.what() << std::endl;
        return 1;
    }
}

/**
 * @brief Main function
 */
int main(int argc, char* argv[]) {
    // Parse command-line arguments
    auto args = CommandLineArgs::parse(argc, argv);
    
    // Show help if requested or invalid args
    if (args.show_help || !args.is_valid()) {
        print_banner();
        print_usage(argv[0]);
        return args.show_help ? 0 : 1;
    }
    
    // Print banner
    print_banner();
    
    // Run optimization
    return run(args);
}
