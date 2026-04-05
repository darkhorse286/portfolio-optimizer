/**
 * @file simulated_annealing_solver.cpp
 * @brief Implementation of SimulatedAnnealingSolver
 */

#include "quantum/simulated_annealing_solver.hpp"
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace portfolio
{
    namespace quantum
    {

        SimulatedAnnealingConfig SimulatedAnnealingConfig::from_json(const nlohmann::json &j)
        {
            SimulatedAnnealingConfig c;
            if (j.contains("max_iterations")) c.max_iterations = j["max_iterations"].get<int>();
            if (j.contains("initial_temperature")) c.initial_temperature = j["initial_temperature"].get<double>();
            if (j.contains("cooling_rate")) c.cooling_rate = j["cooling_rate"].get<double>();
            if (j.contains("random_seed")) c.random_seed = j["random_seed"].get<int>();
            if (j.contains("num_restarts")) c.num_restarts = j["num_restarts"].get<int>();
            if (j.contains("params") && j["params"].is_object()) {
                for (auto it = j["params"].begin(); it != j["params"].end(); ++it) {
                    c.params[it.key()] = it.value().get<double>();
                }
            }
            return c;
        }

        SimulatedAnnealingSolver::SimulatedAnnealingSolver(const SimulatedAnnealingConfig &config)
            : config_(config)
        {
        }

        Eigen::VectorXd SimulatedAnnealingSolver::optimize(
            const Eigen::MatrixXd &covariance,
            const Eigen::VectorXd &expected_returns,
            const portfolio::optimizer::OptimizationConstraints &constraints)
        {
            // read overrides from config_.params
            int num_bits = 4;
            double penalty_budget = 100.0;
            auto itb = config_.params.find("num_bits_per_asset");
            if (itb != config_.params.end()) num_bits = static_cast<int>(itb->second);
            auto itp = config_.params.find("penalty_budget");
            if (itp != config_.params.end()) penalty_budget = itp->second;

            QUBOFormulation qubo(covariance, expected_returns, num_bits, 1.0, penalty_budget);
            qubo.build();
            const Eigen::MatrixXd &Q = qubo.Q();

            int best_obj_index = -1;
            double best_obj = std::numeric_limits<double>::infinity();
            Eigen::VectorXi best_x;

            auto t0 = std::chrono::high_resolution_clock::now();

            for (int r = 0; r < config_.num_restarts; ++r)
            {
                std::mt19937 rng(config_.random_seed + r);
                Eigen::VectorXi x = run_sa(Q, rng);
                double obj = qubo.objective_value(x);
                if (obj < best_obj)
                {
                    best_obj = obj;
                    best_x = x;
                }
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            last_solve_time_ms_ = std::chrono::duration<double, std::milli>(t1 - t0).count();

            // convergence info
            std::ostringstream ss;
            ss << "restarts: " << config_.num_restarts
               << " | best_objective: " << best_obj
               << " | final_temp: " << config_.initial_temperature * std::pow(config_.cooling_rate, config_.max_iterations)
               << " | iters: " << config_.max_iterations;
            last_convergence_info_ = ss.str();

            Eigen::VectorXd weights = qubo.decode(best_x);
            return weights;
        }

        double SimulatedAnnealingSolver::solve_time_ms() const { return last_solve_time_ms_; }
        double SimulatedAnnealingSolver::circuit_execution_us() const { return -1.0; }
        std::string SimulatedAnnealingSolver::convergence_info() const { return last_convergence_info_; }
        std::string SimulatedAnnealingSolver::solver_name() const { return "sa_classical"; }
        std::string SimulatedAnnealingSolver::execution_backend() const { return "sa_classical"; }

        Eigen::VectorXi SimulatedAnnealingSolver::run_sa(const Eigen::MatrixXd &Q, std::mt19937 &rng) const
        {
            int nvar = static_cast<int>(Q.rows());
            Eigen::VectorXi x(nvar);
            x.setZero();

            std::uniform_int_distribution<int> dist_pos(0, nvar - 1);
            std::uniform_real_distribution<double> dist01(0.0, 1.0);

            double T = config_.initial_temperature;

            for (int iter = 0; iter < config_.max_iterations; ++iter)
            {
                int pos = dist_pos(rng);
                double delta = compute_delta_energy(Q, x, pos);
                if (delta < 0.0 || dist01(rng) < std::exp(-delta / std::max(T, 1e-12)))
                {
                    x(pos) = 1 - x(pos);
                }
                T *= config_.cooling_rate;
                if (T < 1e-10) break;
            }
            return x;
        }

        double SimulatedAnnealingSolver::compute_delta_energy(const Eigen::MatrixXd &Q,
                                                              const Eigen::VectorXi &x,
                                                              int pos) const
        {
            // sign = +1 for 0->1 flip, -1 for 1->0 flip
            int bit = x(pos);
            double sign = (bit == 0) ? +1.0 : -1.0;
            Eigen::VectorXd xd = x.cast<double>();
            double col_dot = Q.col(pos).dot(xd);
            double delta = sign * (2.0 * col_dot - Q(pos, pos));
            return delta;
        }

    } // namespace quantum
} // namespace portfolio
