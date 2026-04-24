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

static std::string read_latest_job_id(const std::string &jobs_file)
{
    std::ifstream ifs(jobs_file);
    if (!ifs)
    {
        throw std::runtime_error("Failed to open jobs file: " + jobs_file);
    }

    json jobs_json;
    ifs >> jobs_json;
    auto jobs = jobs_json.at("jobs");
    return jobs.back().at("job_id").get<std::string>();
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

    SECTION("result file is written with correct schema")
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

        QiskitSolver solver(config);
        REQUIRE_NOTHROW(solver.optimize(covariance, expected_returns, constraints));

        std::string job_id = read_latest_job_id(config.jobs_file);
        std::string result_path = config.results_dir + "/quantum_result_" + job_id + ".json";

        std::ifstream ifs(result_path);
        REQUIRE(ifs.good());

        json result_json;
        ifs >> result_json;
        REQUIRE(result_json["schema_version"].get<std::string>() == "1.0");
        REQUIRE(result_json["status"].get<std::string>() == "COMPLETED");
        REQUIRE(result_json["execution_backend"].get<std::string>() == "aer_simulator");
        REQUIRE(result_json["weights"].is_array());
        REQUIRE(result_json["weights"].size() == 10);
        REQUIRE(result_json["bitstring"].is_string());
        REQUIRE(result_json["objective_value"].is_number());
        REQUIRE(result_json["solve_time_ms"].get<double>() >= 0.0);
    }

    SECTION("missing problem file throws runtime_error")
    {
        WARN("Missing-file error path not testable with current QiskitSolver design.");
        SKIP("Skipping missing problem file path due to private result reader design.");
    }
}
