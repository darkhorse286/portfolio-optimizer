/**
 * @file qiskit_solver.cpp
 * @brief Implementation of the Qiskit subprocess wrapper solver.
 */

#include "quantum/qiskit_solver.hpp"
#include <nlohmann/json.hpp>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using json = nlohmann::json;

namespace portfolio
{
    namespace quantum
    {

        QiskitSolverConfig QiskitSolverConfig::from_json(const json &j)
        {
            QiskitSolverConfig cfg;
            if (j.contains("backend")) cfg.backend = j["backend"].get<std::string>();
            if (j.contains("qaoa_depth")) cfg.qaoa_depth = j["qaoa_depth"].get<int>();
            if (j.contains("shots")) cfg.shots = j["shots"].get<int>();
            if (j.contains("problem_file")) cfg.problem_file = j["problem_file"].get<std::string>();
            if (j.contains("jobs_file")) cfg.jobs_file = j["jobs_file"].get<std::string>();
            if (j.contains("results_dir")) cfg.results_dir = j["results_dir"].get<std::string>();
            if (j.contains("timeout_minutes")) cfg.timeout_minutes = j["timeout_minutes"].get<int>();
            if (j.contains("params") && j["params"].is_object()) {
                for (auto it = j["params"].begin(); it != j["params"].end(); ++it) {
                    cfg.params[it.key()] = it.value().get<double>();
                }
            }
            return cfg;
        }

        QiskitSolver::QiskitSolver(const QiskitSolverConfig &config)
            : config_(config)
        {
        }

        Eigen::VectorXd QiskitSolver::optimize(
            const Eigen::MatrixXd &covariance,
            const Eigen::VectorXd &expected_returns,
            const portfolio::optimizer::OptimizationConstraints &constraints)
        {
            if (!python3_available())
            {
                throw std::runtime_error("'python3' not found in PATH; cannot execute Qiskit scripts.");
            }

            int num_bits_per_asset = 4;
            double penalty_budget = 100.0;
            auto bits_it = config_.params.find("num_bits_per_asset");
            if (bits_it != config_.params.end())
            {
                num_bits_per_asset = static_cast<int>(bits_it->second);
            }
            auto penalty_it = config_.params.find("penalty_budget");
            if (penalty_it != config_.params.end())
            {
                penalty_budget = penalty_it->second;
            }

            QUBOFormulation qubo(
                covariance,
                expected_returns,
                num_bits_per_asset,
                1.0,
                penalty_budget);
            qubo.build();
            write_problem_file(qubo);

            auto start_time = std::chrono::high_resolution_clock::now();

            auto quote = [](const std::string &value) {
                return std::string("\"") + value + "\"";
            };

            std::filesystem::path problem_path(config_.problem_file);
            if (problem_path.has_parent_path())
            {
                std::filesystem::create_directories(problem_path.parent_path());
            }

            std::string submit_cmd =
                get_python_executable() + " scripts/quantum/qiskit_submit.py"
                " --problem-file " + quote(config_.problem_file) +
                " --jobs-file " + quote(config_.jobs_file) +
                " --backend " + quote(config_.backend) +
                " --shots " + std::to_string(config_.shots) +
                " --qaoa-depth " + std::to_string(config_.qaoa_depth);
            int rc = std::system(submit_cmd.c_str());
            if (rc != 0)
            {
                throw std::runtime_error("qiskit_submit.py failed with exit code " + std::to_string(rc));
            }

            std::filesystem::path jobs_path(config_.jobs_file);
            std::ifstream jobs_stream(jobs_path);
            if (!jobs_stream)
            {
                throw std::runtime_error("Failed to open jobs file: " + config_.jobs_file);
            }
            json jobs_json;
            jobs_stream >> jobs_json;
            if (!jobs_json.contains("jobs") || !jobs_json["jobs"].is_array() || jobs_json["jobs"].empty())
            {
                throw std::runtime_error("No jobs found in " + config_.jobs_file);
            }
            auto &jobs = jobs_json["jobs"];
            last_job_id_ = jobs.back().value("job_id", "");
            if (last_job_id_.empty())
            {
                throw std::runtime_error("Failed to parse last_job_id from jobs file.");
            }

            std::string collect_cmd =
                get_python_executable() + " scripts/quantum/qiskit_collect.py"
                " --jobs-file " + quote(config_.jobs_file) +
                " --results-dir " + quote(config_.results_dir) +
                " --timeout-min " + std::to_string(config_.timeout_minutes);
            rc = std::system(collect_cmd.c_str());
            if (rc != 0)
            {
                throw std::runtime_error("qiskit_collect.py failed with exit code " + std::to_string(rc));
            }

            auto end_time = std::chrono::high_resolution_clock::now();
            last_solve_time_ms_ = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end_time - start_time).count();

            auto weights = read_result_file();
            return weights;
        }

        double QiskitSolver::solve_time_ms() const
        {
            return last_solve_time_ms_;
        }

        double QiskitSolver::circuit_execution_us() const
        {
            return last_circuit_execution_us_;
        }

        std::string QiskitSolver::convergence_info() const
        {
            return last_convergence_info_;
        }

        std::string QiskitSolver::solver_name() const
        {
            return "qaoa_p" + std::to_string(config_.qaoa_depth) + "_" + config_.backend;
        }

        std::string QiskitSolver::execution_backend() const
        {
            return config_.backend;
        }

        void QiskitSolver::write_problem_file(const QUBOFormulation &qubo) const
        {
            std::filesystem::path path(config_.problem_file);
            if (path.has_parent_path())
            {
                std::filesystem::create_directories(path.parent_path());
            }

            std::ofstream ofs(path);
            if (!ofs)
            {
                throw std::runtime_error("Failed to open problem file: " + config_.problem_file);
            }
            ofs << qubo.to_json();
        }

        Eigen::VectorXd QiskitSolver::read_result_file()
        {
            std::filesystem::path result_path = std::filesystem::path(config_.results_dir) /
                ("quantum_result_" + last_job_id_ + ".json");
            std::ifstream ifs(result_path);
            if (!ifs)
            {
                throw std::runtime_error("Result file not found: " + result_path.string());
            }
            json result_json;
            ifs >> result_json;

            const std::string status = result_json.value("status", "");
            const std::string convergence_info = result_json.value("convergence_info", "");
            last_circuit_execution_us_ = result_json.value("circuit_execution_us", -1.0);
            last_convergence_info_ = convergence_info;

            if (status == "FAILED" || status == "TIMEOUT")
            {
                throw std::runtime_error("Quantum job " + status + ": " + convergence_info);
            }
            if (status != "COMPLETED")
            {
                throw std::runtime_error("Unexpected job status: " + status);
            }

            if (!result_json.contains("weights") || !result_json["weights"].is_array())
            {
                throw std::runtime_error("Invalid result weights in " + result_path.string());
            }

            Eigen::VectorXd weights(result_json["weights"].size());
            for (size_t i = 0; i < result_json["weights"].size(); ++i)
            {
                weights(static_cast<int>(i)) = result_json["weights"][i].get<double>();
            }
            return weights;
        }

        bool QiskitSolver::python3_available()
        {
            std::string cmd = get_python_executable() + " --version >/dev/null 2>&1";
            return std::system(cmd.c_str()) == 0;
        }

        std::string QiskitSolver::get_python_executable()
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

    } // namespace quantum
} // namespace portfolio
