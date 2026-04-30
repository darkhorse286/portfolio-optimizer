/**
 * @file test_qiskit_solver_integration.cpp
 * @brief Integration tests for QiskitSolver via Aer submit/collect.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <nlohmann/json.hpp>
#include <fstream>
#include <cstdlib>
#include <string>

#include "quantum/qiskit_solver.hpp"

using json = nlohmann::json;
using namespace portfolio::quantum;

static std::pair<Eigen::MatrixXd, Eigen::VectorXd> make_10asset_test_data()
{
    Eigen::MatrixXd cov(10, 10);
    cov.setIdentity();
    cov *= 0.02;

    Eigen::VectorXd mu(10);
    for (int i = 0; i < 10; ++i)
    {
        mu(i) = 0.02 + 0.001 * i;
    }

    return {cov, mu};
}

TEST_CASE("QiskitSolver Aer submit and collect round-trip", "[QiskitSolver][integration]")
{
    if (std::getenv("CI_SKIP_INTEGRATION") != nullptr)
    {
        SKIP("CI_SKIP_INTEGRATION set; skipping Qiskit integration test.");
    }

    if (std::system("command -v python3 >/dev/null 2>&1") != 0)
    {
        WARN("python3 not available; skipping Qiskit integration test");
        SKIP();
    }

    if (std::system("python3 -c 'import qiskit' >/dev/null 2>&1") != 0)
    {
        WARN("qiskit not available; skipping Qiskit integration test");
        SKIP();
    }

    SECTION("10-asset problem produces valid weights via Aer")
    {
        auto [covariance, expected_returns] = make_10asset_test_data();
        portfolio::optimizer::OptimizationConstraints constraints;

        QiskitSolverConfig config;
        config.backend = "aer_simulator";
        config.qaoa_depth = 1;
        config.shots = 512;
        config.problem_file = "/tmp/test_quantum_problem.json";
        config.jobs_file = "/tmp/test_quantum_jobs.json";
        config.results_dir = "/tmp";
        config.params["num_bits_per_asset"] = 2;  // 10 assets * 2 bits = 20 qubits; 4 bits * 10 = 40 qubits exceeds statevector memory

        QiskitSolver solver(config);
        Eigen::VectorXd weights;
        REQUIRE_NOTHROW(weights = solver.optimize(covariance, expected_returns, constraints));

        REQUIRE(weights.size() == 10);
        REQUIRE(weights.sum() == Catch::Approx(1.0).epsilon(0.15));
        for (int i = 0; i < weights.size(); ++i)
        {
            REQUIRE(weights(i) >= -0.01);
        }
        REQUIRE(solver.solve_time_ms() > 0.0);
        REQUIRE(std::abs(solver.circuit_execution_us() + 1.0) < 1e-12);
        REQUIRE_FALSE(solver.convergence_info().empty());
        REQUIRE(solver.execution_backend() == "aer_simulator");
        REQUIRE(solver.solver_name().find("aer_simulator") != std::string::npos);
    }

    SECTION("problem file is written with correct schema")
    {
        // For Aer, the solver communicates synchronously over stdout — no result
        // file or jobs file is written. The only persisted artifact is the QUBO
        // problem file. Test that it is present and well-formed.
        auto [covariance, expected_returns] = make_10asset_test_data();
        portfolio::optimizer::OptimizationConstraints constraints;

        QiskitSolverConfig config;
        config.backend = "aer_simulator";
        config.qaoa_depth = 1;
        config.shots = 512;
        config.problem_file = "/tmp/test_quantum_problem.json";
        config.jobs_file = "/tmp/test_quantum_jobs.json";
        config.results_dir = "/tmp";
        config.params["num_bits_per_asset"] = 2;  // 10 assets * 2 bits = 20 qubits

        QiskitSolver solver(config);
        REQUIRE_NOTHROW(solver.optimize(covariance, expected_returns, constraints));

        std::ifstream ifs(config.problem_file);
        REQUIRE(ifs.good());

        json problem_json;
        ifs >> problem_json;
        REQUIRE(problem_json.contains("num_assets"));
        REQUIRE(problem_json["num_assets"].get<int>() == 10);
        REQUIRE(problem_json.contains("num_bits_per_asset"));
        REQUIRE(problem_json.contains("num_variables"));
        REQUIRE(problem_json.contains("Q"));
        REQUIRE(problem_json["Q"].is_array());
    }

}
