/**
 * @file qubo_formulation.hpp
 * @brief QUBO formulation for mean-variance portfolio optimization
 * @details
 * Encoding:
 *   w_i = sum_{k=0}^{B-1} 2^k * x_{i,k} / (2^B - 1)
 * Budget constraint sum(w)=1 is encoded as penalty: P * (sum(w) - 1)^2
 * Resulting objective is a quadratic form in binary variables: x^T Q x
 */

#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace portfolio
{
    namespace quantum
    {

        /**
         * @struct QUBOConfig
         * @brief Simple QUBO configuration
         */
        struct QUBOConfig {
            int num_bits_per_asset = 4;
            double penalty_budget = 100.0;
            double penalty_turnover = 0.0;
        };

        /**
         * @class QUBOFormulation
         * @brief Builds a binary quadratic matrix Q from mean-variance inputs
         */
        class QUBOFormulation
        {
        public:
            /**
             * @brief Construct QUBO formulation
             * @param covariance Covariance matrix (N x N)
             * @param expected_returns Expected returns vector (N)
             * @param num_bits_per_asset Bits per asset (B)
             * @param risk_aversion Risk aversion scalar
             * @param penalty_budget Budget penalty coefficient P
             * @throws std::invalid_argument on invalid inputs
             */
            explicit QUBOFormulation(
                const Eigen::MatrixXd &covariance,
                const Eigen::VectorXd &expected_returns,
                int num_bits_per_asset,
                double risk_aversion,
                double penalty_budget);

            /** Build the Q matrix. */
            void build();

            /** @brief Access the Q matrix; throws if build() wasn't called */
            const Eigen::MatrixXd &Q() const;

            /** @brief Decode binary solution to continuous weights */
            Eigen::VectorXd decode(const Eigen::VectorXi &binary_solution) const;

            /** @brief Encode continuous weights to binary solution (quantize) */
            Eigen::VectorXi encode(const Eigen::VectorXd &weights) const;

            /** @brief Number of binary variables (N * B) */
            int num_variables() const;

            /** @brief Objective value of binary vector x (x^T Q x) */
            double objective_value(const Eigen::VectorXi &x) const;

            /** @brief Serialize problem to JSON for Python handoff */
            std::string to_json() const;

            /** Accessors */
            int num_assets() const { return static_cast<int>(expected_returns_.size()); }
            int num_bits() const { return num_bits_per_asset_; }

        private:
            Eigen::MatrixXd covariance_;
            Eigen::VectorXd expected_returns_;
            int num_bits_per_asset_;
            double risk_aversion_;
            double penalty_budget_;
            Eigen::MatrixXd Q_;
            bool built_ = false;
        };

    } // namespace quantum
} // namespace portfolio
