/**
 * @file simulated_annealing_solver.hpp
 * @brief Simulated annealing QUBO solver (classical)
 * @details Classical simulated annealing solver used to validate QUBO
 * formulation. Produces binary solutions that are decoded to continuous
 * weights via the QUBOFormulation::decode method.
 */

#pragma once

#include "quantum/quantum_optimizer.hpp"
#include "quantum/qubo_formulation.hpp"
#include "optimizer/optimizer_interface.hpp"

#include <random>
#include <string>
#include <map>

namespace portfolio
{
    namespace quantum
    {

        /**
         * @struct SimulatedAnnealingConfig
         * @brief Configuration for simulated annealing solver
         */
        struct SimulatedAnnealingConfig {
            int max_iterations = 50000;
            double initial_temperature = 1000.0;
            double cooling_rate = 0.999;
            int random_seed = 42;
            int num_restarts = 10;
            std::map<std::string, double> params; ///< optional overrides (e.g., num_bits_per_asset, penalty_budget)

            static SimulatedAnnealingConfig from_json(const nlohmann::json &j);
        };

        /**
         * @class SimulatedAnnealingSolver
         * @brief Simulated annealing-based QUBO solver
         */
        class SimulatedAnnealingSolver : public QuantumOptimizer
        {
        public:
            /** @brief Construct with config */
            explicit SimulatedAnnealingSolver(const SimulatedAnnealingConfig &config);

            Eigen::VectorXd optimize(
                const Eigen::MatrixXd &covariance,
                const Eigen::VectorXd &expected_returns,
                const portfolio::optimizer::OptimizationConstraints &constraints) override;

            double solve_time_ms() const override;
            double circuit_execution_us() const override; // always -1.0
            std::string convergence_info() const override;
            std::string solver_name() const override;       // "sa_classical"
            std::string execution_backend() const override; // "sa_classical"

        private:
            SimulatedAnnealingConfig config_;
            double last_solve_time_ms_ = 0.0;
            std::string last_convergence_info_;

            Eigen::VectorXi run_sa(const Eigen::MatrixXd &Q, std::mt19937 &rng) const;
            double compute_delta_energy(const Eigen::MatrixXd &Q,
                                        const Eigen::VectorXi &x,
                                        int pos) const;
        };

    } // namespace quantum
} // namespace portfolio
