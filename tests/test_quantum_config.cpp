#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <data/data_loader.hpp>
#include <quantum/qiskit_solver.hpp>
#include <fstream>
#include <filesystem>

using json = nlohmann::json;

class QuantumConfigTest : public ::testing::Test {
protected:
    std::string temp_path;

    void SetUp() override {
        temp_path = std::filesystem::temp_directory_path() /
                    ("test_quantum_config_" + std::to_string(::getpid()) + ".json");
    }

    void TearDown() override {
        std::filesystem::remove(temp_path);
    }

    void write_config(const json& j) {
        std::ofstream f(temp_path);
        f << j.dump();
    }
};

TEST_F(QuantumConfigTest, ParsesQuantumSolversArray) {
    json cfg = {
        {"data", {{"universe", {"AAPL", "MSFT"}}, {"data_file", "test.csv"}}},
        {"optimizer", {{"risk_free_rate", 0.02}}},
        {"risk_model", {{"type", "ewma"}, {"estimation_window", 252}, {"ewma_lambda", 0.94}}},
        {"backtest", {{"initial_capital", 1000000.0}, {"rebalance_frequency", "monthly"}}},
        {"quantum_solvers", {{
            {"name", "Test QAOA Solver"},
            {"type", "quantum"},
            {"backend", "aer_simulator"},
            {"algorithm_mode", "qaoa"},
            {"augment_problem_data", true},
            {"optimization_level", 2},
            {"qaoa_depth", 2},
            {"shots", 2048},
            {"problem_file", "test_problem.json"},
            {"jobs_file", "test_jobs.json"},
            {"results_dir", "test_results"}
        }}}
    };
    write_config(cfg);

    portfolio::PortfolioConfig config;
    ASSERT_NO_THROW(config = portfolio::DataLoader::load_config(temp_path));
    ASSERT_EQ(config.quantum_solvers.size(), 1u);

    const auto& solver = config.quantum_solvers[0];
    EXPECT_EQ(solver.name, "Test QAOA Solver");
    EXPECT_EQ(solver.type, "quantum");
    EXPECT_EQ(solver.backend, "aer_simulator");
    EXPECT_EQ(solver.algorithm_mode, "qaoa");
    EXPECT_EQ(solver.augment_problem_data, true);
    EXPECT_EQ(solver.optimization_level, 2);
    EXPECT_EQ(solver.qaoa_depth, 2);
    EXPECT_EQ(solver.shots, 2048);
    EXPECT_EQ(solver.problem_file, "test_problem.json");
    EXPECT_EQ(solver.jobs_file, "test_jobs.json");
    EXPECT_EQ(solver.results_dir, "test_results");
}

TEST_F(QuantumConfigTest, DefaultFieldsFallbackToDefaults) {
    json cfg = {
        {"data", {{"universe", {"AAPL"}}, {"data_file", "x.csv"}}},
        {"optimizer", {}},
        {"risk_model", {}},
        {"backtest", {}},
        {"quantum_solvers", {{{"name", "Minimal"}}}}
    };
    write_config(cfg);

    portfolio::PortfolioConfig config;
    ASSERT_NO_THROW(config = portfolio::DataLoader::load_config(temp_path));
    ASSERT_EQ(config.quantum_solvers.size(), 1u);

    const auto& e = config.quantum_solvers[0];
    EXPECT_EQ(e.name, "Minimal");
    EXPECT_EQ(e.type, "quantum");
    EXPECT_EQ(e.backend, "aer_simulator");
    EXPECT_EQ(e.algorithm_mode, "qaoa");
    EXPECT_EQ(e.augment_problem_data, false);
    EXPECT_EQ(e.optimization_level, 1);
    EXPECT_EQ(e.qaoa_depth, 1);
    EXPECT_EQ(e.shots, 1024);
    EXPECT_TRUE(e.problem_file.empty());
    EXPECT_TRUE(e.jobs_file.empty());
    EXPECT_EQ(e.results_dir, "results");
}

TEST_F(QuantumConfigTest, QiskitSolverConfigFromEntry) {
    portfolio::PortfolioConfig::QuantumSolverEntry entry;
    entry.name = "Test Solver";
    entry.backend = "aer_simulator";
    entry.algorithm_mode = "qaoa";
    entry.augment_problem_data = true;
    entry.optimization_level = 3;
    entry.qaoa_depth = 1;
    entry.shots = 1024;
    entry.problem_file = "problem.json";
    entry.jobs_file = "jobs.json";
    entry.results_dir = "results";

    auto cfg = portfolio::quantum::QiskitSolverConfig::from_entry(entry);

    EXPECT_EQ(cfg.name, "Test Solver");
    EXPECT_EQ(cfg.backend, "aer_simulator");
    EXPECT_EQ(cfg.algorithm_mode, "qaoa");
    EXPECT_EQ(cfg.augment_problem_data, true);
    EXPECT_EQ(cfg.optimization_level, 3);
    EXPECT_EQ(cfg.qaoa_depth, 1);
    EXPECT_EQ(cfg.shots, 1024);
    EXPECT_EQ(cfg.problem_file, "problem.json");
    EXPECT_EQ(cfg.jobs_file, "jobs.json");
    EXPECT_EQ(cfg.results_dir, "results");
    EXPECT_EQ(cfg.params.at("num_bits_per_asset"), 2.0);
}
