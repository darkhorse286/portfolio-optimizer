/**
 * @file qiskit_solver.cpp
 * @brief Implementation of the Qiskit subprocess wrapper solver.
 */

#include "quantum/qiskit_solver.hpp"
#include <nlohmann/json.hpp>
#include <chrono>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#ifndef _WIN32
#include <poll.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
#include <random>
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
            if (j.contains("timeout_seconds")) cfg.timeout_seconds = j["timeout_seconds"].get<int>();
            if (j.contains("worker_script")) cfg.worker_script = j["worker_script"].get<std::string>();
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

        QiskitSolver::~QiskitSolver()
        {
            std::lock_guard<std::mutex> guard(worker_mutex_);
            shutdown_worker();
        }

        Eigen::VectorXd QiskitSolver::optimize(
            const Eigen::MatrixXd &covariance,
            const Eigen::VectorXd &expected_returns,
            const portfolio::optimizer::OptimizationConstraints &constraints)
        {
            if (!python3_available())
            {
                throw std::runtime_error("'python3' not found in PATH; cannot execute Qiskit worker.");
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

            std::filesystem::path problem_path(config_.problem_file);
            if (problem_path.has_parent_path())
            {
                std::filesystem::create_directories(problem_path.parent_path());
            }
            write_problem_file(qubo);

            auto start_time = std::chrono::high_resolution_clock::now();

            json request = {
                {"request_id", generate_request_id()},
                {"mode", config_.algorithm_mode},
                {"backend", config_.backend},
                {"qaoa_depth", config_.qaoa_depth},
                {"shots", config_.shots},
                {"problem_file", config_.problem_file},
                {"augment_problem_data", config_.augment_problem_data}
            };
            if (config_.params.count("num_bits_per_asset") > 0)
            {
                request["num_bits_per_asset"] = config_.params.at("num_bits_per_asset");
            }
            if (config_.params.count("penalty_budget") > 0)
            {
                request["penalty_budget"] = config_.params.at("penalty_budget");
            }

            json response;
            {
                std::lock_guard<std::mutex> guard(worker_mutex_);
                response = send_request(request);
            }

            if (!response.contains("status"))
            {
                throw std::runtime_error("Invalid response from Qiskit worker: missing status field.");
            }
            if (response["status"].get<std::string>() != "ok")
            {
                std::string error_msg = response.value("error", "unknown error");
                throw std::runtime_error("Qiskit worker failed: " + error_msg);
            }
            if (!response.contains("weights") || !response["weights"].is_array())
            {
                throw std::runtime_error("Invalid response from Qiskit worker: missing weights array.");
            }

            last_solve_time_ms_ = response.value("solve_time_ms", 0.0);
            if (last_solve_time_ms_ <= 0.0)
            {
                last_solve_time_ms_ = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
                    std::chrono::high_resolution_clock::now() - start_time).count();
            }
            last_circuit_execution_us_ = response.value("circuit_execution_us", -1.0);
            last_convergence_info_ = response.value("convergence_info", "");
            last_job_id_ = response.value("job_id", "");

            Eigen::VectorXd weights(response["weights"].size());
            for (size_t i = 0; i < response["weights"].size(); ++i)
            {
                weights(static_cast<int>(i)) = response["weights"][i].get<double>();
            }
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
            return config_.algorithm_mode + "_p" + std::to_string(config_.qaoa_depth) + "_" + config_.backend;
        }

        std::string QiskitSolver::execution_backend() const
        {
            return config_.backend;
        }

        std::string QiskitSolver::generate_request_id()
        {
            static const char *hex_chars = "0123456789abcdef";
            std::random_device rd;
            std::uniform_int_distribution<uint8_t> dist(0, 255);
            std::string result;
            result.reserve(32);
            for (int i = 0; i < 16; ++i)
            {
                uint8_t value = dist(rd);
                result.push_back(hex_chars[value >> 4]);
                result.push_back(hex_chars[value & 0x0F]);
            }
            return result;
        }

        void QiskitSolver::ensure_worker_running()
        {
#ifdef _WIN32
            throw std::runtime_error("Persistent Qiskit worker is not supported on Windows.");
#else
            if (worker_pid_ > 0 && worker_stdin_ != nullptr && worker_stdout_ != nullptr)
            {
                return;
            }

            if (!std::filesystem::exists(config_.worker_script))
            {
                throw std::runtime_error("Qiskit worker script not found: " + config_.worker_script);
            }

            int stdin_pipe[2];
            int stdout_pipe[2];
            if (pipe(stdin_pipe) != 0 || pipe(stdout_pipe) != 0)
            {
                throw std::runtime_error("Failed to create pipes for Qiskit worker: " + std::string(std::strerror(errno)));
            }

            pid_t pid = fork();
            if (pid < 0)
            {
                close(stdin_pipe[0]);
                close(stdin_pipe[1]);
                close(stdout_pipe[0]);
                close(stdout_pipe[1]);
                throw std::runtime_error("Failed to fork Qiskit worker: " + std::string(std::strerror(errno)));
            }

            if (pid == 0)
            {
                if (dup2(stdin_pipe[0], STDIN_FILENO) < 0 || dup2(stdout_pipe[1], STDOUT_FILENO) < 0)
                {
                    _exit(127);
                }
                close(stdin_pipe[0]);
                close(stdin_pipe[1]);
                close(stdout_pipe[0]);
                close(stdout_pipe[1]);
                const std::string python_exe = get_python_executable();
                execlp(python_exe.c_str(), python_exe.c_str(), config_.worker_script.c_str(), (char *)nullptr);
                _exit(127);
            }

            close(stdin_pipe[0]);
            close(stdout_pipe[1]);

            worker_stdin_ = fdopen(stdin_pipe[1], "w");
            worker_stdout_ = fdopen(stdout_pipe[0], "r");
            if (!worker_stdin_ || !worker_stdout_)
            {
                if (worker_stdin_) fclose(worker_stdin_);
                if (worker_stdout_) fclose(worker_stdout_);
                kill(pid, SIGTERM);
                waitpid(pid, nullptr, 0);
                worker_pid_ = -1;
                throw std::runtime_error("Failed to open Qiskit worker streams.");
            }
            setvbuf(worker_stdin_, nullptr, _IOLBF, 0);
            setvbuf(worker_stdout_, nullptr, _IOLBF, 0);
            worker_pid_ = pid;
#endif
        }

        void QiskitSolver::restart_worker()
        {
            shutdown_worker();
            ensure_worker_running();
        }

        nlohmann::json QiskitSolver::send_request(const nlohmann::json &req)
        {
            ensure_worker_running();

            std::string payload = req.dump() + "\n";
            if (std::fputs(payload.c_str(), worker_stdin_) == EOF || std::fflush(worker_stdin_) != 0)
            {
                restart_worker();
                if (std::fputs(payload.c_str(), worker_stdin_) == EOF || std::fflush(worker_stdin_) != 0)
                {
                    throw std::runtime_error("Failed to write request to Qiskit worker.");
                }
            }

            int stdout_fd = fileno(worker_stdout_);
            struct pollfd poll_fd;
            poll_fd.fd = stdout_fd;
            poll_fd.events = POLLIN;
            poll_fd.revents = 0;
            int timeout_ms = config_.timeout_seconds * 1000;
            int ready = poll(&poll_fd, 1, timeout_ms);
            if (ready < 0)
            {
                throw std::runtime_error("Failed to poll Qiskit worker stdout: " + std::string(std::strerror(errno)));
            }
            if (ready == 0)
            {
                restart_worker();
                throw std::runtime_error("Timeout waiting for Qiskit worker response after " + std::to_string(config_.timeout_seconds) + " seconds.");
            }

            std::string response_line;
            char buffer[4096];
            while (true)
            {
                if (!std::fgets(buffer, sizeof(buffer), worker_stdout_))
                {
                    restart_worker();
                    throw std::runtime_error("Failed to read response from Qiskit worker.");
                }
                response_line.append(buffer);
                if (!response_line.empty() && response_line.back() == '\n')
                {
                    break;
                }
            }
            try
            {
                return json::parse(response_line);
            }
            catch (const std::exception &ex)
            {
                throw std::runtime_error(std::string("Failed to parse JSON response from Qiskit worker: ") + ex.what());
            }
        }

        void QiskitSolver::shutdown_worker()
        {
#ifdef _WIN32
            (void)worker_stdin_;
            (void)worker_stdout_;
            (void)worker_pid_;
#else
            if (worker_stdin_)
            {
                fclose(worker_stdin_);
                worker_stdin_ = nullptr;
            }
            if (worker_stdout_)
            {
                fclose(worker_stdout_);
                worker_stdout_ = nullptr;
            }
            if (worker_pid_ > 0)
            {
                kill(worker_pid_, SIGTERM);
                waitpid(worker_pid_, nullptr, 0);
                worker_pid_ = -1;
            }
#endif
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
