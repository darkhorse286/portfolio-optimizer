/**
 * @file main.cpp
 * @brief Main entry point for Quantum Portfolio Optimizer
 *
 * Command-line application that loads configuration, estimates risk models,
 * runs portfolio optimization, and generates results.
 */

#include "data/data_loader.hpp"
#include "data/market_data.hpp"
#include "risk/risk_model_factory.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
#include "optimizer/efficient_frontier.hpp"
#include <iostream>
#include <string>
#include <exception>
#include <iomanip>
#include <chrono>

using namespace portfolio;

/**
 * @brief Print usage information
 */
void print_usage(const char *program_name)
{
    std::cout << "Quantum Portfolio Optimizer v1.1.0\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  --config PATH         Path to configuration JSON file (required)\n"
              << "  --output PATH         Path to output directory (default: results/)\n"
              << "  --frontier            Compute efficient frontier\n"
              << "  --verbose             Enable verbose logging\n"
              << "  --help, -h            Show this help message\n"
              << "\nExample:\n"
              << "  " << program_name << " --config data/config/portfolio_config.json --verbose\n"
              << "  " << program_name << " --config data/config/portfolio_config.json --frontier\n"
              << std::endl;
}

/**
 * @brief Print banner
 */
void print_banner()
{
    std::cout << "\n"
              << "================================================================\n"
              << "       Quantum Portfolio Optimizer v1.1.0                      \n"
              << "       Modern C++17 Portfolio Optimization                     \n"
              << "       Feature Set 2.1: Mean-Variance Optimization             \n"
              << "================================================================\n"
              << std::endl;
}

/**
 * @brief Parse command-line arguments
 */
struct CommandLineArgs
{
    std::string config_path;
    std::string output_dir = "results";
    bool compute_frontier = false;
    bool verbose = false;
    bool show_help = false;

    static CommandLineArgs parse(int argc, char *argv[])
    {
        CommandLineArgs args;

        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];

            if (arg == "--help" || arg == "-h")
            {
                args.show_help = true;
            }
            else if (arg == "--config" && i + 1 < argc)
            {
                args.config_path = argv[++i];
            }
            else if (arg == "--output" && i + 1 < argc)
            {
                args.output_dir = argv[++i];
            }
            else if (arg == "--frontier")
            {
                args.compute_frontier = true;
            }
            else if (arg == "--verbose")
            {
                args.verbose = true;
            }
            else
            {
                std::cerr << "Warning: Unknown argument: " << arg << std::endl;
            }
        }

        return args;
    }

    bool is_valid() const
    {
        return !show_help && !config_path.empty();
    }
};

/**
 * @brief Print optimization result
 */
void print_optimization_result(
    const std::string &title,
    const optimizer::OptimizationResult &result,
    const std::vector<std::string> &tickers)
{
    std::cout << "\n"
              << title << "\n";
    std::cout << std::string(60, '-') << "\n";

    if (!result.success)
    {
        std::cout << "Optimization failed: " << result.message << "\n";
        return;
    }

    std::cout << "Status: SUCCESS\n";
    std::cout << "Iterations: " << result.iterations << "\n\n";

    std::cout << "Portfolio Statistics:\n";
    std::cout << "  Expected Return:  " << std::fixed << std::setprecision(2)
              << result.expected_return * 100 << "% (annualized: "
              << result.expected_return * 252 * 100 << "%)\n";
    std::cout << "  Volatility:       " << result.volatility * 100
              << "% (annualized: " << result.volatility * std::sqrt(252) * 100 << "%)\n";
    std::cout << "  Sharpe Ratio:     " << std::setprecision(3)
              << result.sharpe_ratio << "\n\n";

    std::cout << "Asset Weights:\n";

    // Sort by weight for better display
    std::vector<std::pair<double, std::string>> sorted_weights;
    for (size_t i = 0; i < tickers.size(); ++i)
    {
        if (std::abs(result.weights(i)) > 1e-4)
        {
            sorted_weights.push_back({result.weights(i), tickers[i]});
        }
    }
    std::sort(sorted_weights.begin(), sorted_weights.end(),
              [](const auto &a, const auto &b)
              { return a.first > b.first; });

    for (const auto &pair : sorted_weights)
    {
        std::cout << "  " << std::setw(8) << std::left << pair.second
                  << ": " << std::fixed << std::setprecision(2)
                  << pair.first * 100 << "%\n";
    }

    std::cout << std::string(60, '-') << "\n";
}

/**
 * @brief Main execution function
 */
int run(const CommandLineArgs &args)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    try
    {
        // ====================================================================
        // 1. Load Configuration
        // ====================================================================
        std::cout << "[1/6] Loading configuration..." << std::endl;

        auto config = DataLoader::load_config(args.config_path);

        if (args.verbose)
        {
            std::cout << "  - Tickers: ";
            for (const auto &ticker : config.data.universe)
            {
                std::cout << ticker << " ";
            }
            std::cout << "\n  - Date range: " << config.data.start_date
                      << " to " << config.data.end_date << "\n";
            std::cout << "  - Risk model: " << config.risk_model.type << "\n";
            std::cout << "  - Optimizer: " << config.optimizer.type << "\n";
        }

        // ====================================================================
        // 2. Load Market Data
        // ====================================================================
        std::cout << "[2/6] Loading market data..." << std::endl;

        auto data = DataLoader::load_csv(config.data.data_file, config.data.universe);

        // Filter by date range if specified
        if (!config.data.start_date.empty() && !config.data.end_date.empty())
        {
            data = data.filter_by_date(config.data.start_date, config.data.end_date);
        }

        std::cout << "  - Loaded " << data.num_dates() << " dates, "
                  << data.num_assets() << " assets" << std::endl;

        if (args.verbose)
        {
            data.print_summary();
        }

        // ====================================================================
        // 3. Calculate Returns
        // ====================================================================
        std::cout << "[3/6] Calculating returns..." << std::endl;

        auto returns = data.calculate_returns(ReturnType::SIMPLE);
        auto mean_returns = data.mean_returns(ReturnType::SIMPLE);
        auto volatilities = data.volatilities(ReturnType::SIMPLE);

        if (args.verbose)
        {
            std::cout << "\n  Asset Statistics (Daily):\n";
            std::cout << "  " << std::string(50, '-') << "\n";
            std::cout << "  " << std::setw(8) << "Ticker"
                      << std::setw(12) << "Mean Ret"
                      << std::setw(12) << "Volatility" << "\n";
            std::cout << "  " << std::string(50, '-') << "\n";

            for (size_t i = 0; i < data.num_assets(); ++i)
            {
                std::cout << "  " << std::setw(8) << data.get_tickers()[i]
                          << std::setw(11) << std::fixed << std::setprecision(4)
                          << mean_returns(i) * 100 << "%"
                          << std::setw(11) << volatilities(i) * 100 << "%\n";
            }
            std::cout << "  " << std::string(50, '-') << "\n";
        }

        // ====================================================================
        // 4. Risk Model Estimation
        // ====================================================================
        std::cout << "[4/6] Estimating risk model..." << std::endl;

        auto risk_model = risk::RiskModelFactory::create(config.risk_model);
        auto covariance = risk_model->estimate_covariance(returns);

        if (args.verbose)
        {
            std::cout << "  - Risk model: " << risk_model->get_name() << "\n";
            std::cout << "  - Covariance matrix: " << covariance.rows()
                      << "x" << covariance.cols() << "\n";

            // Check condition number
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(covariance);
            auto sv = svd.singularValues();
            double cond = sv(0) / sv(sv.size() - 1);
            std::cout << "  - Condition number: " << std::scientific
                      << std::setprecision(2) << cond << "\n";
        }

        // ====================================================================
        // 5. Portfolio Optimization
        // ====================================================================
        std::cout << "[5/6] Running portfolio optimization..." << std::endl;

        // Set up constraints from config
        optimizer::OptimizationConstraints constraints;
        constraints.min_weight = config.optimizer.min_weight;
        constraints.max_weight = config.optimizer.max_weight;
        constraints.long_only = config.optimizer.long_only;
        constraints.sum_to_one = config.optimizer.sum_to_one;
        constraints.max_turnover = config.optimizer.max_turnover;

        if (args.verbose)
        {
            std::cout << "\n  Constraints:\n";
            std::cout << "  - Weight range: [" << constraints.min_weight
                      << ", " << constraints.max_weight << "]\n";
            std::cout << "  - Long only: " << (constraints.long_only ? "Yes" : "No") << "\n";
            std::cout << "  - Sum to one: " << (constraints.sum_to_one ? "Yes" : "No") << "\n";
        }

        // 5.1: Minimum Variance Portfolio
        std::cout << "\n  Computing minimum variance portfolio...\n";
        optimizer::MeanVarianceOptimizer min_var_opt(
            optimizer::ObjectiveType::MIN_VARIANCE,
            config.optimizer.risk_free_rate);

        auto min_var_result = min_var_opt.optimize(
            mean_returns, covariance, constraints);

        print_optimization_result(
            "MINIMUM VARIANCE PORTFOLIO",
            min_var_result,
            data.get_tickers());

        // 5.2: Maximum Sharpe Ratio Portfolio
        std::cout << "\n  Computing maximum Sharpe ratio portfolio...\n";
        optimizer::MeanVarianceOptimizer max_sharpe_opt(
            optimizer::ObjectiveType::MAX_SHARPE,
            config.optimizer.risk_free_rate);

        auto max_sharpe_result = max_sharpe_opt.optimize(
            mean_returns, covariance, constraints);

        print_optimization_result(
            "MAXIMUM SHARPE RATIO PORTFOLIO",
            max_sharpe_result,
            data.get_tickers());

        // 5.3: Risk Aversion Portfolio
        std::cout << "\n  Computing risk aversion portfolio...\n";
        optimizer::MeanVarianceOptimizer risk_aversion_opt(
            optimizer::ObjectiveType::RISK_AVERSION,
            config.optimizer.risk_free_rate);

        risk_aversion_opt.set_risk_aversion(config.optimizer.risk_aversion);

        auto risk_aversion_result = risk_aversion_opt.optimize(
            mean_returns, covariance, constraints);

        print_optimization_result(
            "RISK AVERSION PORTFOLIO (lambda=" +
                std::to_string(config.optimizer.risk_aversion) + ")",
            risk_aversion_result,
            data.get_tickers());

        // ====================================================================
        // 6. Efficient Frontier (Optional)
        // ====================================================================
        if (args.compute_frontier)
        {
            std::cout << "\n[6/6] Computing efficient frontier..." << std::endl;

            optimizer::EfficientFrontier frontier;
            frontier.set_num_points(20);

            if (args.verbose)
            {
                std::cout << "  - Computing 20 frontier points...\n";
            }

            auto frontier_result = frontier.compute(
                mean_returns,
                covariance,
                constraints,
                config.optimizer.risk_free_rate);

            if (frontier_result.success)
            {
                frontier_result.print_summary();

                // Export to CSV
                std::string frontier_file = args.output_dir + "/efficient_frontier.csv";
                frontier_result.export_to_csv(frontier_file);
                std::cout << "\n  Frontier data exported to: " << frontier_file << "\n";
            }
            else
            {
                std::cout << "  Frontier computation failed: "
                          << frontier_result.message << "\n";
            }
        }
        else
        {
            std::cout << "\n[6/6] Skipping efficient frontier (use --frontier to enable)\n";
        }

        // ====================================================================
        // Summary
        // ====================================================================
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                            end_time - start_time)
                            .count();

        std::cout << "\n================================================================\n";
        std::cout << "Optimization completed successfully in "
                  << duration << " ms\n";
        std::cout << "================================================================\n"
                  << std::endl;

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
}

/**
 * @brief Main function
 */
int main(int argc, char *argv[])
{
    // Parse command-line arguments
    auto args = CommandLineArgs::parse(argc, argv);

    // Show help if requested or invalid args
    if (args.show_help || !args.is_valid())
    {
        print_banner();
        print_usage(argv[0]);
        return args.show_help ? 0 : 1;
    }

    // Print banner
    print_banner();

    // Run optimization
    return run(args);
}