/**
 * @file qubo_formulation.cpp
 * @brief Implementation of QUBOFormulation
 */

#include "quantum/qubo_formulation.hpp"
#include <stdexcept>
#include <sstream>
#include <chrono>

using nlohmann::json;

namespace portfolio
{
    namespace quantum
    {

        static double pow2(int k) { return static_cast<double>(1ULL << k); }

        QUBOFormulation::QUBOFormulation(
            const Eigen::MatrixXd &covariance,
            const Eigen::VectorXd &expected_returns,
            int num_bits_per_asset,
            double risk_aversion,
            double penalty_budget)
            : covariance_(covariance), expected_returns_(expected_returns),
              num_bits_per_asset_(num_bits_per_asset),
              risk_aversion_(risk_aversion), penalty_budget_(penalty_budget)
        {
            // Validation
            if (covariance_.rows() != covariance_.cols())
                throw std::invalid_argument("covariance must be square");
            if ((covariance_ - covariance_.transpose()).cwiseAbs().maxCoeff() > 1e-8)
                throw std::invalid_argument("covariance must be symmetric");
            if (expected_returns_.size() != covariance_.rows())
                throw std::invalid_argument("expected_returns size mismatch");
            if (num_bits_per_asset_ < 1 || num_bits_per_asset_ > 8)
                throw std::invalid_argument("num_bits_per_asset must be in [1,8]");
            if (risk_aversion_ <= 0.0)
                throw std::invalid_argument("risk_aversion must be > 0");
            if (penalty_budget_ <= 0.0)
                throw std::invalid_argument("penalty_budget must be > 0");
        }

        void QUBOFormulation::build()
        {
            int n = covariance_.rows();
            int B = num_bits_per_asset_;
            int Nvar = n * B;
            Q_.setZero(Nvar, Nvar);

            double denom = static_cast<double>((1ULL << B) - 1ULL);

            // 1) Risk term
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    double cov_ij = covariance_(i, j);
                    for (int k = 0; k < B; ++k)
                    {
                        for (int l = 0; l < B; ++l)
                        {
                            int p = i * B + k;
                            int q = j * B + l;
                            double coeff = 0.5 * cov_ij * (pow2(k) * pow2(l)) / (denom * denom);
                            Q_(p, q) += coeff;
                        }
                    }
                }
            }

            // 2) Return term (negated) multiplied by risk_aversion
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < B; ++k)
                {
                    int p = i * B + k;
                    double coeff = risk_aversion_ * expected_returns_(i) * (pow2(k) / denom);
                    Q_(p, p) -= coeff;
                }
            }

            // 3) Budget penalty: P * (sum(w) - 1)^2
            // Expand: P*(sum_i sum_k c_ik x_ik)^2 - 2P*(sum_i sum_k c_ik x_ik) + P
            // where c_ik = 2^k / denom
            std::vector<double> c(Nvar);
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < B; ++k)
                {
                    int p = i * B + k;
                    c[p] = pow2(k) / denom;
                }
            }

            // Quadratic part
            for (int p = 0; p < Nvar; ++p)
            {
                for (int q = 0; q < Nvar; ++q)
                {
                    Q_(p, q) += penalty_budget_ * c[p] * c[q];
                }
            }

            // Linear part contributes to diagonal as -2P*c_p
            for (int p = 0; p < Nvar; ++p)
            {
                Q_(p, p) -= 2.0 * penalty_budget_ * c[p];
            }

            // constant term penalty_budget_ is ignored

            // Symmetrize
            Q_ = 0.5 * (Q_ + Q_.transpose());

            built_ = true;
        }

        const Eigen::MatrixXd &QUBOFormulation::Q() const
        {
            if (!built_)
                throw std::runtime_error("Q accessed before build()");
            return Q_;
        }

        Eigen::VectorXd QUBOFormulation::decode(const Eigen::VectorXi &binary_solution) const
        {
            if (!built_)
                throw std::runtime_error("decode called before build()");
            int n = covariance_.rows();
            int B = num_bits_per_asset_;
            if (binary_solution.size() != n * B)
                throw std::invalid_argument("binary_solution size mismatch");

            Eigen::VectorXd w(n);
            w.setZero();
            double denom = static_cast<double>((1ULL << B) - 1ULL);
            for (int i = 0; i < n; ++i)
            {
                double sum = 0.0;
                for (int k = 0; k < B; ++k)
                {
                    int p = i * B + k;
                    int bit = binary_solution(p);
                    sum += static_cast<double>(bit) * pow2(k);
                }
                w(i) = sum / denom;
                if (w(i) < 0.0) w(i) = 0.0;
                if (w(i) > 1.0) w(i) = 1.0;
            }
            return w;
        }

        Eigen::VectorXi QUBOFormulation::encode(const Eigen::VectorXd &weights) const
        {
            if (!built_)
                throw std::runtime_error("encode called before build()");
            int n = covariance_.rows();
            int B = num_bits_per_asset_;
            if (weights.size() != n)
                throw std::invalid_argument("weights size mismatch");

            Eigen::VectorXi x(n * B);
            x.setZero();
            double denom = static_cast<double>((1ULL << B) - 1ULL);
            for (int i = 0; i < n; ++i)
            {
                double clipped = std::min(1.0, std::max(0.0, weights(i)));
                unsigned int q = static_cast<unsigned int>(std::round(clipped * denom));
                for (int k = 0; k < B; ++k)
                {
                    int p = i * B + k;
                    x(p) = (q >> k) & 1u;
                }
            }
            return x;
        }

        int QUBOFormulation::num_variables() const
        {
            return static_cast<int>(expected_returns_.size()) * num_bits_per_asset_;
        }

        double QUBOFormulation::objective_value(const Eigen::VectorXi &x) const
        {
            if (!built_)
                throw std::runtime_error("objective_value called before build()");
            Eigen::VectorXd xd = x.cast<double>();
            return xd.dot(Q_ * xd);
        }

        std::string QUBOFormulation::to_json() const
        {
            if (!built_)
                throw std::runtime_error("to_json called before build()");
            json j;
            j["schema_version"] = "1.0";
            // simple timestamp id
            auto now = std::chrono::system_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
            j["problem_id"] = std::to_string(ms);
            int n = covariance_.rows();
            j["num_assets"] = n;
            j["num_bits_per_asset"] = num_bits_per_asset_;
            j["num_variables"] = num_variables();
            j["risk_aversion"] = risk_aversion_;
            j["penalty_budget"] = penalty_budget_;

            // Q as nested arrays
            std::vector<std::vector<double>> Qarr(Q_.rows(), std::vector<double>(Q_.cols()));
            for (int i = 0; i < Q_.rows(); ++i)
                for (int k = 0; k < Q_.cols(); ++k)
                    Qarr[i][k] = Q_(i, k);
            j["Q"] = Qarr;

            // expected_returns
            std::vector<double> mu(n);
            for (int i = 0; i < n; ++i) mu[i] = expected_returns_(i);
            j["expected_returns"] = mu;

            // asset names placeholder
            std::vector<std::string> names(n);
            for (int i = 0; i < n; ++i) names[i] = "asset_" + std::to_string(i);
            j["asset_names"] = names;

            return j.dump();
        }

    } // namespace quantum
} // namespace portfolio
