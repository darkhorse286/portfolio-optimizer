#include <catch2/catch_test_macros.hpp>
#include "data/data_loader.hpp"
#include "quantum/qiskit_solver.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <filesystem>
#include <string>
#include <unistd.h>

namespace fs = std::filesystem;
using json = nlohmann::json;

static std::string write_temp_config(const json& j) {
    auto tmp = fs::temp_directory_path() /
               ("test_quantum_cfg_" + std::to_string(::getpid()) + ".json");
    std::ofstream f(tmp);
    f << j.dump(2);
    return tmp.string();
}

TEST_CASE("quantum_solvers array is parsed into QuantumSolverEntry", "[quantum_config]") {
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
    const std::string path = write_temp_config(cfg);

    portfolio::PortfolioConfig config;
    REQUIRE_NOTHROW(config = portfolio::DataLoader::load_config(path));
    fs::remove(path);

    REQUIRE(config.quantum_solvers.size() == 1u);
    const auto& s = config.quantum_solvers[0];
    CHECK(s.name               == "Test QAOA Solver");
    CHECK(s.type               == "quantum");
    CHECK(s.backend            == "aer_simulator");
    CHECK(s.algorithm_mode     == "qaoa");
    CHECK(s.augment_problem_data == true);
    CHECK(s.optimization_level == 2);
    CHECK(s.qaoa_depth         == 2);
    CHECK(s.shots              == 2048);
    CHECK(s.problem_file       == "test_problem.json");
    CHECK(s.jobs_file          == "test_jobs.json");
    CHECK(s.results_dir        == "test_results");
}

TEST_CASE("missing optional fields fall back to struct defaults", "[quantum_config]") {
    json cfg = {
        {"data", {{"universe", {"AAPL"}}, {"data_file", "x.csv"}}},
        {"optimizer", json::object()},
        {"risk_model", {{"type", "ewma"}}},
        {"backtest", json::object()},
        {"quantum_solvers", {{{"name", "Minimal"}}}}
    };
    const std::string path = write_temp_config(cfg);

    portfolio::PortfolioConfig config;
    REQUIRE_NOTHROW(config = portfolio::DataLoader::load_config(path));
    fs::remove(path);

    REQUIRE(config.quantum_solvers.size() == 1u);
    const auto& e = config.quantum_solvers[0];
    CHECK(e.name               == "Minimal");
    CHECK(e.type               == "quantum");
    CHECK(e.backend            == "aer_simulator");
    CHECK(e.algorithm_mode     == "qaoa");
    CHECK(e.augment_problem_data == false);
    CHECK(e.optimization_level == 1);
    CHECK(e.qaoa_depth         == 1);
    CHECK(e.shots              == 1024);
    CHECK(e.problem_file.empty());
    CHECK(e.jobs_file.empty());
    CHECK(e.results_dir        == "results");
    CHECK(e.ibm_backend_selection == "shortest_queue");
    CHECK(e.ibm_min_qubits     == 20);
    CHECK(e.n_frontier_points  == 20);
}

TEST_CASE("ibm_auto fields are parsed and forwarded to QiskitSolverConfig", "[quantum_config]") {
    json cfg = {
        {"data", {{"universe", {"AAPL"}}, {"data_file", "x.csv"}}},
        {"optimizer", json::object()},
        {"risk_model", {{"type", "ewma"}}},
        {"backtest", json::object()},
        {"quantum_solvers", {
            {
                {"name", "QAOA (IBM)"},
                {"backend", "ibm_auto"},
                {"algorithm_mode", "qaoa"},
                {"augment_problem_data", false},
                {"optimization_level", 3},
                {"qaoa_depth", 1},
                {"shots", 1024},
                {"ibm_min_qubits", 20},
                {"ibm_backend_selection", "shortest_queue"},
                {"problem_file", "quantum_problem_qaoa_ibm.json"},
                {"jobs_file", "quantum_jobs_qaoa_ibm.json"}
            },
            {
                {"name", "QAMOO (IBM)"},
                {"backend", "ibm_auto"},
                {"algorithm_mode", "qamoo"},
                {"augment_problem_data", true},
                {"optimization_level", 3},
                {"qaoa_depth", 1},
                {"shots", 1024},
                {"ibm_min_qubits", 27},
                {"ibm_backend_selection", "shortest_queue"},
                {"n_frontier_points", 15},
                {"problem_file", "quantum_problem_qamoo_ibm.json"},
                {"jobs_file", "quantum_jobs_qamoo_ibm.json"}
            }
        }}
    };
    const std::string path = write_temp_config(cfg);

    portfolio::PortfolioConfig config;
    REQUIRE_NOTHROW(config = portfolio::DataLoader::load_config(path));
    fs::remove(path);

    REQUIRE(config.quantum_solvers.size() == 2u);

    const auto& qaoa = config.quantum_solvers[0];
    CHECK(qaoa.backend               == "ibm_auto");
    CHECK(qaoa.ibm_min_qubits        == 20);
    CHECK(qaoa.ibm_backend_selection == "shortest_queue");
    CHECK(qaoa.n_frontier_points     == 20);  // default

    const auto& qamoo = config.quantum_solvers[1];
    CHECK(qamoo.backend               == "ibm_auto");
    CHECK(qamoo.algorithm_mode        == "qamoo");
    CHECK(qamoo.augment_problem_data  == true);
    CHECK(qamoo.ibm_min_qubits        == 27);
    CHECK(qamoo.ibm_backend_selection == "shortest_queue");
    CHECK(qamoo.n_frontier_points     == 15);

    // Verify all three fields survive QuantumSolverEntry -> QiskitSolverConfig conversion
    const auto solver_cfg = portfolio::quantum::QiskitSolverConfig::from_entry(qamoo);
    CHECK(solver_cfg.backend               == "ibm_auto");
    CHECK(solver_cfg.ibm_min_qubits        == 27);
    CHECK(solver_cfg.ibm_backend_selection == "shortest_queue");
    CHECK(solver_cfg.n_frontier_points     == 15);
}

TEST_CASE("QiskitSolverConfig::from_entry maps all fields correctly", "[quantum_config]") {
    portfolio::PortfolioConfig::QuantumSolverEntry entry;
    entry.name                  = "Test Solver";
    entry.backend               = "aer_simulator";
    entry.algorithm_mode        = "qaoa";
    entry.augment_problem_data  = true;
    entry.optimization_level    = 3;
    entry.qaoa_depth            = 1;
    entry.shots                 = 1024;
    entry.problem_file          = "problem.json";
    entry.jobs_file             = "jobs.json";
    entry.results_dir           = "results";
    entry.ibm_backend_selection = "shortest_queue";
    entry.ibm_min_qubits        = 20;
    entry.n_frontier_points     = 20;

    const auto cfg = portfolio::quantum::QiskitSolverConfig::from_entry(entry);

    CHECK(cfg.name                  == "Test Solver");
    CHECK(cfg.backend               == "aer_simulator");
    CHECK(cfg.algorithm_mode        == "qaoa");
    CHECK(cfg.augment_problem_data  == true);
    CHECK(cfg.optimization_level    == 3);
    CHECK(cfg.qaoa_depth            == 1);
    CHECK(cfg.shots                 == 1024);
    CHECK(cfg.problem_file          == "problem.json");
    CHECK(cfg.jobs_file             == "jobs.json");
    CHECK(cfg.results_dir           == "results");
    CHECK(cfg.ibm_backend_selection == "shortest_queue");
    CHECK(cfg.ibm_min_qubits        == 20);
    CHECK(cfg.n_frontier_points     == 20);
    CHECK(cfg.params.at("num_bits_per_asset") == 2.0);
}
