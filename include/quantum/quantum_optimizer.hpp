/**
 * @file quantum_optimizer.hpp
 * @brief Abstract base for quantum optimizers
 * @details Defines the interface for quantum optimization solvers used in
 * portfolio optimization. This interface is intentionally separate from the
 * classical `OptimizerInterface` to allow quantum-specific metadata such as
 * circuit timings and backend identifiers.
 */

#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <map>
#include "optimizer/optimizer_interface.hpp"

namespace portfolio
{
    namespace quantum
    {

        /**
         * @struct QuantumOptimizerConfig
         * @brief Configuration container for quantum optimizers
         */
        struct QuantumOptimizerConfig
        {
            std::string solver_type;                        ///< human-readable solver type
            int max_iterations = 10000;                     ///< algorithm iterations
            int random_seed = 42;                           ///< RNG seed
            int num_bits_per_asset = 4;                     ///< bits per asset
            std::map<std::string, double> params;           ///< solver-specific overrides

            /**
             * @brief Construct from JSON
             * @param j JSON object
             * @return QuantumOptimizerConfig
             */
            static QuantumOptimizerConfig from_json(const nlohmann::json &j);
        };

        inline QuantumOptimizerConfig QuantumOptimizerConfig::from_json(const nlohmann::json &j)
        {
            QuantumOptimizerConfig cfg;
            if (j.contains("solver_type")) cfg.solver_type = j["solver_type"].get<std::string>();
            if (j.contains("max_iterations")) cfg.max_iterations = j["max_iterations"].get<int>();
            if (j.contains("random_seed")) cfg.random_seed = j["random_seed"].get<int>();
            if (j.contains("num_bits_per_asset")) cfg.num_bits_per_asset = j["num_bits_per_asset"].get<int>();
            if (j.contains("params") && j["params"].is_object()) {
                for (auto it = j["params"].begin(); it != j["params"].end(); ++it) {
                    cfg.params[it.key()] = it.value().get<double>();
                }
            }
            return cfg;
        }

        /**
         * @class QuantumOptimizer
         * @brief Abstract base class for quantum-style optimizers
         * @details Exposes an `optimize` method similar to the classical
         * OptimizerInterface but adds quantum-specific metadata accessors.
         */
        class QuantumOptimizer
        {
        public:
            virtual ~QuantumOptimizer() = default;

            /**
             * @brief Run optimization
             * @param covariance Covariance matrix (N x N)
             * @param expected_returns Expected returns (N)
             * @param constraints Optimization constraints (classical)
             * @return Decoded continuous weights vector (N)
             */
            virtual Eigen::VectorXd optimize(
                const Eigen::MatrixXd &covariance,
                const Eigen::VectorXd &expected_returns,
                const portfolio::optimizer::OptimizationConstraints &constraints) = 0;

            /** @brief Wall-clock solve time in milliseconds (includes I/O/queue) */
            virtual double solve_time_ms() const = 0;

            /** @brief Circuit execution time on QPU in microseconds; -1.0 if N/A */
            virtual double circuit_execution_us() const = 0;

            /** @brief Human-readable convergence information */
            virtual std::string convergence_info() const = 0;

            /** @brief Solver name identifier */
            virtual std::string solver_name() const = 0;

            /** @brief Execution backend identifier (e.g. "sa_classical", "qiskit_aer") */
            virtual std::string execution_backend() const = 0;
        };

    } // namespace quantum
} // namespace portfolio
