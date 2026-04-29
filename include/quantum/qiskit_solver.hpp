/**
 * @file qiskit_solver.hpp
 * @brief Quantum QAOA solver wrapper using Qiskit subprocess scripts.
 */

#pragma once

#include "optimizer/optimizer_interface.hpp"
#include "quantum/quantum_optimizer.hpp"
#include "quantum/qubo_formulation.hpp"
#include <nlohmann/json.hpp>
#include <map>
#include <mutex>
#include <string>
#include <sys/types.h>
#include <cstdio>

namespace portfolio
{
    namespace quantum
    {

        /**
         * @struct QiskitSolverConfig
         * @brief Configuration container for the Qiskit solver.
         */
        struct QiskitSolverConfig
        {
            std::string backend = "aer_simulator";
            std::string algorithm_mode = "qaoa";  // "qaoa" | "qamo"
            int qaoa_depth = 1;
            int shots = 1024;
            std::string problem_file = "results/quantum_problem.json";
            std::string jobs_file = "results/quantum_jobs.json";
            std::string results_dir = "results";
            int timeout_minutes = 30;
            int timeout_seconds = 120;
            std::string worker_script = "scripts/quantum/qiskit_worker.py";
            bool augment_problem_data = false;  // If true, overwrites expected_returns and covariance
                                              // in the problem file with values computed from
                                              // historical prices before circuit submission.
                                              // Enables informed quantum optimization.
            std::map<std::string, double> params;

            /**
             * @brief Construct config from JSON.
             * @param j JSON object
             * @return QiskitSolverConfig
             */
            static QiskitSolverConfig from_json(const nlohmann::json &j);
        };

        /**
         * @class QiskitSolver
         * @brief Wraps Python Qiskit submission and collection scripts.
         */
        class QiskitSolver : public QuantumOptimizer
        {
        public:
            /** @brief Construct with config. */
            explicit QiskitSolver(const QiskitSolverConfig &config);
            ~QiskitSolver() override;

            Eigen::VectorXd optimize(
                const Eigen::MatrixXd &covariance,
                const Eigen::VectorXd &expected_returns,
                const portfolio::optimizer::OptimizationConstraints &constraints) override;

            double solve_time_ms() const override;
            double circuit_execution_us() const override;
            std::string convergence_info() const override;
            std::string solver_name() const override;
            std::string execution_backend() const override;

        private:
            QiskitSolverConfig config_;
            double last_solve_time_ms_ = 0.0;
            double last_circuit_execution_us_ = -1.0;
            std::string last_convergence_info_;
            std::string last_job_id_;

            void write_problem_file(const QUBOFormulation &qubo) const;
            Eigen::VectorXd read_result_file();
            static bool python3_available();
            static std::string get_python_executable();

            // Persistent Python worker bridge
            FILE* worker_stdin_ = nullptr;
            FILE* worker_stdout_ = nullptr;
            pid_t worker_pid_ = -1;
            std::mutex worker_mutex_;

            void ensure_worker_running();
            void restart_worker();
            nlohmann::json send_request(const nlohmann::json& req);
            void shutdown_worker();
            static std::string generate_request_id();
        };

    } // namespace quantum
} // namespace portfolio
