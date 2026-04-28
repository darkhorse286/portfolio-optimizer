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
#include "backtest/backtest_engine.hpp"
#include "quantum/benchmark_runner.hpp"
#include "quantum/simulated_annealing_solver.hpp"
#include "quantum/qiskit_solver.hpp"
#include <iostream>
#include <string>
#include <exception>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <filesystem>

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
              << "  --benchmark           Run quantum benchmark (MV + SA + Aer QAOA + Aer QAMO)\n"
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
    bool generate_report = false;
    bool quantum_submit = false;
    bool quantum_collect = false;
    bool benchmark = false;

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
            else if (arg == "--report")
            {
                args.generate_report = true;
            }
            else if (arg == "--quantum-submit")
            {
                args.quantum_submit = true;
            }
            else if (arg == "--quantum-collect")
            {
                args.quantum_collect = true;
            }
            else if (arg == "--benchmark")
            {
                args.benchmark = true;
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
 * @brief Return the Python executable to use for subprocesses.
 *
 * Prefers the interpreter inside the active virtualenv (VIRTUAL_ENV env var)
 * so that packages installed via `pip install -r requirements.txt` are found
 * even though std::system() spawns a bare shell without venv activation.
 */
static std::string get_python_executable()
{
    const char* venv = std::getenv("VIRTUAL_ENV");
    if (venv)
    {
        std::string venv_python = std::string(venv) + "/bin/python3";
        if (std::filesystem::exists(venv_python))
            return venv_python;
    }
    return "python3";
}

/**
 * @brief Check if the resolved Python interpreter is available.
 */
static bool python3_available()
{
    std::string cmd = get_python_executable() + " --version >/dev/null 2>&1";
    return std::system(cmd.c_str()) == 0;
}

/**
 * @brief Run the Python report generator script using `std::system`.
 *
 * This function checks for `python3` first; if not available it prints a
 * warning to `std::cerr` and returns without crashing the program.
 *
 * @param output_dir Directory to pass as both input and output to the script.
 */
static void run_report_generator(const std::string &output_dir)
{
    if (!python3_available())
    {
        std::cerr << "Warning: 'python3' not found in PATH; skipping report generation.\n";
        return;
    }

    std::cout << "Starting report generation (this may take a moment)...\n";

    std::string cmd = get_python_executable() + " scripts/generate_report.py --input-dir " + output_dir + " --output-dir " + output_dir;
    int rc = std::system(cmd.c_str());

    if (rc == 0)
    {
        std::cout << "Report generation completed successfully." << std::endl;
    }
    else
    {
        std::cerr << "Warning: report generation command exited with code " << rc << "." << std::endl;
    }
}

/**
 * @brief Run the Qiskit submit script via subprocess.
 *
 * @param output_dir Directory for problem and jobs JSON files.
 */
static void run_quantum_submit(const std::string &output_dir)
{
    if (!python3_available())
    {
        std::cerr << "Warning: 'python3' not found in PATH; skipping quantum submit." << std::endl;
        return;
    }

    std::filesystem::create_directories(output_dir);
    std::string problem_file = output_dir + "/quantum_problem.json";
    std::string jobs_file = output_dir + "/quantum_jobs.json";

    std::cout << "Starting quantum submit...\n";
    std::string cmd = get_python_executable() + " scripts/quantum/qiskit_submit.py"
                      " --problem-file " + std::string("\"") + problem_file + "\""
                      " --jobs-file " + std::string("\"") + jobs_file + "\"";
    int rc = std::system(cmd.c_str());
    if (rc == 0)
    {
        std::cout << "Quantum submit completed successfully." << std::endl;
    }
    else
    {
        std::cerr << "Warning: quantum submit command exited with code " << rc << "." << std::endl;
    }
}

/**
 * @brief Run the Qiskit collect script via subprocess.
 *
 * @param output_dir Directory for jobs and result JSON files.
 */
static void run_quantum_collect(const std::string &output_dir)
{
    if (!python3_available())
    {
        std::cerr << "Warning: 'python3' not found in PATH; skipping quantum collect." << std::endl;
        return;
    }

    std::filesystem::create_directories(output_dir);
    std::string jobs_file = output_dir + "/quantum_jobs.json";

    std::cout << "Starting quantum collect...\n";
    std::string cmd = get_python_executable() + " scripts/quantum/qiskit_collect.py"
                      " --jobs-file " + std::string("\"") + jobs_file + "\""
                      " --results-dir " + std::string("\"") + output_dir + "\"";
    int rc = std::system(cmd.c_str());
    if (rc == 0)
    {
        std::cout << "Quantum collect completed successfully." << std::endl;
    }
    else
    {
        std::cerr << "Warning: quantum collect command exited with code " << rc << "." << std::endl;
    }
}

/**
 * @brief Invoke benchmark_viz.py to render comparison_results.json into an HTML report.
 *
 * @param output_dir Directory containing comparison_results.json and receiving outputs.
 */
static void run_benchmark_viz(const std::string &output_dir)
{
    if (!python3_available())
    {
        std::cerr << "Warning: 'python3' not found in PATH; skipping benchmark visualization.\n";
        return;
    }

    std::cout << "Generating benchmark report...\n";
    std::string comparison_file = output_dir + "/comparison_results.json";
    std::string cmd = get_python_executable() + " scripts/quantum/benchmark_viz.py"
                      " --comparison-file \"" + comparison_file + "\""
                      " --output-dir \"" + output_dir + "\"";
    int rc = std::system(cmd.c_str());
    if (rc == 0)
    {
        std::cout << "Benchmark report written to: " << output_dir << "/quantum_report.html\n";
    }
    else
    {
        std::cerr << "Warning: benchmark_viz.py exited with code " << rc << ".\n";
    }
}

/**
 * @brief Run the four-way quantum benchmark (MV + SA + Aer QAOA + Aer QAMO) and generate the report.
 *
 * @param config      Loaded portfolio configuration.
 * @param data        Market data already loaded and filtered.
 * @param output_dir  Directory to write comparison_results.json and report outputs.
 */
static void run_benchmark(
    const portfolio::PortfolioConfig &config,
    const portfolio::MarketData &data,
    const std::string &output_dir)
{
    std::filesystem::create_directories(output_dir);

    // Remove stale jobs files so this run starts with a clean slate.
    // Without this, every previous failed/successful job accumulates in the file
    // and collect re-processes them all on each rebalancing period.
    std::filesystem::remove(output_dir + "/quantum_jobs_qaoa.json");
    std::filesystem::remove(output_dir + "/quantum_jobs_qaoa_informed.json");
    std::filesystem::remove(output_dir + "/quantum_jobs_qamo.json");
    std::filesystem::remove(output_dir + "/quantum_jobs_qamo_informed.json");

    auto backtest_params = portfolio::backtest::BacktestParams::from_config(config);

    portfolio::quantum::BenchmarkRunner runner(data, backtest_params);

    // Classical baseline — max-Sharpe mean-variance
    auto mv_solver = std::make_shared<portfolio::optimizer::MeanVarianceOptimizer>(
        portfolio::optimizer::ObjectiveType::MAX_SHARPE,
        config.optimizer.risk_free_rate);
    runner.add_classical_solver("markowitz_mv", mv_solver);

    // Quantum-inspired — simulated annealing (use struct defaults)
    portfolio::quantum::SimulatedAnnealingConfig sa_config;
    auto sa_solver = std::make_shared<portfolio::quantum::SimulatedAnnealingSolver>(sa_config);
    runner.add_quantum_solver("sa_classical", "quantum_inspired", sa_solver);

    // Quantum — QAOA on Aer simulator
    // Preflight: verify qiskit is importable before registering the solver.
    // Without this, a missing qiskit produces one error per rebalancing day.
    bool qiskit_ok = std::system(
        (get_python_executable() + " -c \"import qiskit\" >/dev/null 2>&1").c_str()) == 0;
    if (!qiskit_ok)
    {
        std::cerr << "Warning: 'qiskit' not found in Python environment ("
                  << get_python_executable() << "). "
                  << "Skipping qaoa_p1_aer_simulator solver. "
                  << "Ensure the venv is activated and requirements are installed.\n";
    }
    else
    {
        // num_bits_per_asset=2 keeps the circuit at 10×2=20 qubits (16 MB statevector).
        // The default of 4 bits produces 40 qubits (16 TB), which is not simulable locally.

        // QAOA uninformed — raw QUBO Q matrix, no market data augmentation
        portfolio::quantum::QiskitSolverConfig qaoa_cfg;
        qaoa_cfg.backend               = "aer_simulator";
        qaoa_cfg.algorithm_mode        = "qaoa";
        qaoa_cfg.augment_problem_data  = false;
        qaoa_cfg.qaoa_depth            = 1;
        qaoa_cfg.shots                 = 1024;
        qaoa_cfg.problem_file          = output_dir + "/quantum_problem_qaoa.json";
        qaoa_cfg.jobs_file             = output_dir + "/quantum_jobs_qaoa.json";
        qaoa_cfg.results_dir           = output_dir;
        qaoa_cfg.params["num_bits_per_asset"] = 2.0;
        auto qaoa_solver = std::make_shared<portfolio::quantum::QiskitSolver>(qaoa_cfg);
        runner.add_quantum_solver("qaoa_uninformed_aer", "quantum", qaoa_solver);

        // QAOA informed — augmented with EWMA expected_returns and covariance
        portfolio::quantum::QiskitSolverConfig qaoa_informed_cfg;
        qaoa_informed_cfg.backend              = "aer_simulator";
        qaoa_informed_cfg.algorithm_mode       = "qaoa";
        qaoa_informed_cfg.augment_problem_data = true;
        qaoa_informed_cfg.qaoa_depth           = 1;
        qaoa_informed_cfg.shots                = 1024;
        qaoa_informed_cfg.problem_file         = output_dir + "/quantum_problem_qaoa_informed.json";
        qaoa_informed_cfg.jobs_file            = output_dir + "/quantum_jobs_qaoa_informed.json";
        qaoa_informed_cfg.results_dir          = output_dir;
        qaoa_informed_cfg.params["num_bits_per_asset"] = 2.0;
        auto qaoa_informed_solver = std::make_shared<portfolio::quantum::QiskitSolver>(qaoa_informed_cfg);
        runner.add_quantum_solver("qaoa_informed_aer", "quantum", qaoa_informed_solver);

        // QAMO uninformed
        portfolio::quantum::QiskitSolverConfig qamo_cfg;
        qamo_cfg.backend               = "aer_simulator";
        qamo_cfg.algorithm_mode        = "qamo";
        qamo_cfg.augment_problem_data  = false;
        qamo_cfg.qaoa_depth            = 1;
        qamo_cfg.shots                 = 1024;
        qamo_cfg.problem_file          = output_dir + "/quantum_problem_qamo.json";
        qamo_cfg.jobs_file             = output_dir + "/quantum_jobs_qamo.json";
        qamo_cfg.results_dir           = output_dir;
        qamo_cfg.params["num_bits_per_asset"] = 2.0;
        auto qamo_solver = std::make_shared<portfolio::quantum::QiskitSolver>(qamo_cfg);
        runner.add_quantum_solver("qamo_uninformed_aer", "quantum", qamo_solver);

        // QAMO informed — current production behavior
        portfolio::quantum::QiskitSolverConfig qamo_informed_cfg;
        qamo_informed_cfg.backend              = "aer_simulator";
        qamo_informed_cfg.algorithm_mode       = "qamo";
        qamo_informed_cfg.augment_problem_data = true;
        qamo_informed_cfg.qaoa_depth           = 1;
        qamo_informed_cfg.shots                = 1024;
        qamo_informed_cfg.problem_file         = output_dir + "/quantum_problem_qamo_informed.json";
        qamo_informed_cfg.jobs_file            = output_dir + "/quantum_jobs_qamo_informed.json";
        qamo_informed_cfg.results_dir          = output_dir;
        qamo_informed_cfg.params["num_bits_per_asset"] = 2.0;
        auto qamo_informed_solver = std::make_shared<portfolio::quantum::QiskitSolver>(qamo_informed_cfg);
        runner.add_quantum_solver("qamo_informed_aer", "quantum", qamo_informed_solver);
    }

    std::cout << "Running benchmark...\n";
    auto result = runner.run(3);

    result.print_summary();

    std::string comparison_path = output_dir + "/comparison_results.json";
    result.export_comparison_json(comparison_path);
    std::cout << "Comparison results written to: " << comparison_path << "\n";

    run_benchmark_viz(output_dir);
}

/**
 * @brief Main execution function
 */
int run(const CommandLineArgs &args)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    try
    {
        if (args.quantum_submit || args.quantum_collect)
        {
            if (args.quantum_submit)
            {
                run_quantum_submit(args.output_dir);
            }
            if (args.quantum_collect)
            {
                run_quantum_collect(args.output_dir);
            }
            return 0;
        }

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

                // Create output directory if it doesn't exist
                std::filesystem::create_directories(args.output_dir);
                
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

        // Optionally generate the HTML report if requested
        if (args.generate_report)
        {
            run_report_generator(args.output_dir);
        }

        if (args.quantum_submit)
        {
            run_quantum_submit(args.output_dir);
        }

        if (args.quantum_collect)
        {
            run_quantum_collect(args.output_dir);
        }

        if (args.benchmark)
        {
            run_benchmark(config, data, args.output_dir);
        }

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