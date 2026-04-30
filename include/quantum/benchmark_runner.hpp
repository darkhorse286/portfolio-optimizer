/**
 * @file benchmark_runner.hpp
 * @brief Benchmark runner for comparing classical and quantum optimizers
 */

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "data/market_data.hpp"
#include "backtest/backtest_engine.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
#include "quantum/quantum_optimizer.hpp"
#include "analytics/performance_metrics.hpp"

namespace portfolio
{
    namespace quantum
    {

        /**
         * @struct SolverResult
         * @brief Result from running a single solver
         */
        struct SolverResult
        {
            std::string solver_name;
            std::string solver_type;          // "classical" | "quantum_inspired" | "quantum"
            std::string execution_backend;
            int problem_size;                 // number of assets in this run
            double mean_sharpe;
            double std_sharpe;
            double mean_solve_time_ms;
            double mean_circuit_execution_us; // -1.0 for non-hardware runs
            double mean_portfolio_return;
            double mean_portfolio_volatility;
            double solution_quality_vs_classical; // solver_sharpe / markowitz_sharpe
                                        // 1.0 for the classical baseline itself
            int num_runs;
            int successful_runs;
            std::vector<double> nav_series;
            std::vector<std::string> date_series;
            portfolio::backtest::TradeSummary trade_summary;
            std::vector<double> avg_weights;
            std::vector<std::string> tickers;
            std::unordered_map<std::string, double> sector_weights;
            std::unordered_map<std::string, double> benchmark_sector_weights;
        };

        /**
         * @struct BenchmarkResult
         * @brief Aggregate result from benchmark run
         */
        struct BenchmarkResult
        {
            std::string run_id;               // timestamp string, e.g. "20260405_143022"
            std::vector<int> problem_sizes;
            std::vector<SolverResult> results;

            std::string to_json() const;
            void print_summary() const;       // scaling table to stdout
            void export_comparison_json(const std::string& path) const;
        };

        /**
         * @class BenchmarkRunner
         * @brief Orchestrates benchmark runs across multiple solvers
         */
        class BenchmarkRunner
        {
        public:
            BenchmarkRunner(const portfolio::MarketData& data,
                            const portfolio::backtest::BacktestParams& backtest_params);

            void add_classical_solver(
                const std::string& name,
                std::shared_ptr<portfolio::optimizer::MeanVarianceOptimizer> solver);

            void add_quantum_solver(
                const std::string& name,
                const std::string& solver_type,   // "quantum_inspired" | "quantum"
                std::shared_ptr<portfolio::quantum::QuantumOptimizer> solver);

            BenchmarkResult run(int num_runs = 1);

        private:
            const portfolio::MarketData& data_;
            portfolio::backtest::BacktestParams backtest_params_;

            struct RegisteredSolver
            {
                std::string name;
                std::string solver_type;
                std::string execution_backend;
                // one of these is set, the other is nullptr:
                std::shared_ptr<portfolio::optimizer::MeanVarianceOptimizer> classical;
                std::shared_ptr<portfolio::quantum::QuantumOptimizer> quantum;
            };
            std::vector<RegisteredSolver> solvers_;

            SolverResult run_solver(const RegisteredSolver& s, int num_runs) const;
        };

    } // namespace quantum
} // namespace portfolio