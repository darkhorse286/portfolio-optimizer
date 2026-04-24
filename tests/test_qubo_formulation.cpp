/**
 * @file test_qubo_formulation.cpp
 * @brief Unit tests for QUBOFormulation
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>

#include "quantum/qubo_formulation.hpp"

using namespace portfolio::quantum;

class QUBOTestFixture
{
protected:
    QUBOFormulation *small_; // 2 assets, 2 bits
    QUBOFormulation *large_; // 10 assets, 4 bits

    QUBOTestFixture()
    {
        // 2-asset simple case
        Eigen::MatrixXd cov2(2, 2);
        cov2 << 0.01, 0.002,
                0.002, 0.02;
        Eigen::VectorXd mu2(2);
        mu2 << 0.05, 0.03;
        small_ = new QUBOFormulation(cov2, mu2, 2, 1.0, 10.0);

        // 10-asset random but deterministic
        Eigen::MatrixXd cov10 = Eigen::MatrixXd::Identity(10, 10) * 0.05;
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < i; ++j)
            {
                double v = 0.001 * (i + j + 1);
                cov10(i, j) = v;
                cov10(j, i) = v;
            }
        Eigen::VectorXd mu10(10);
        for (int i = 0; i < 10; ++i) mu10(i) = 0.02 + 0.001 * i;
        large_ = new QUBOFormulation(cov10, mu10, 4, 1.0, 100.0);
    }

    ~QUBOTestFixture()
    {
        delete small_;
        delete large_;
    }
};

TEST_CASE_METHOD(QUBOTestFixture, "Q matrix is square with correct dimensions", "[QUBOFormulation]")
{
    SECTION("small dimensions")
    {
        small_->build();
        auto Q = small_->Q();
        REQUIRE(Q.rows() == 4);
        REQUIRE(Q.cols() == 4);
    }

    SECTION("large dimensions")
    {
        large_->build();
        auto Q = large_->Q();
        REQUIRE(Q.rows() == 40);
        REQUIRE(Q.cols() == 40);
    }
}

TEST_CASE_METHOD(QUBOTestFixture, "Q matrix is symmetric after build", "[QUBOFormulation]")
{
    small_->build();
    auto Q = small_->Q();
    REQUIRE((Q - Q.transpose()).cwiseAbs().maxCoeff() < 1e-10);
}

TEST_CASE_METHOD(QUBOTestFixture, "encode-decode round-trip within quantization tolerance", "[QUBOFormulation]")
{
    small_->build();
    Eigen::VectorXd w(2);
    w << 0.5, 0.5;
    auto x = small_->encode(w);
    auto w2 = small_->decode(x);
    double tol = 1.0 / (static_cast<double>((1ULL << 2) - 1ULL));
    for (int i = 0; i < 2; ++i)
        REQUIRE(std::abs(w(i) - w2(i)) <= tol + 1e-12);

    w << 0.3, 0.7;
    x = small_->encode(w);
    w2 = small_->decode(x);
    for (int i = 0; i < 2; ++i)
        REQUIRE(std::abs(w(i) - w2(i)) <= tol + 1e-12);
}

TEST_CASE_METHOD(QUBOTestFixture, "objective_value matches direct computation", "[QUBOFormulation]")
{
    small_->build();
    auto Q = small_->Q();
    Eigen::VectorXi x(4);
    x << 1, 0, 1, 1;
    double obj = small_->objective_value(x);
    Eigen::VectorXd xd = x.cast<double>();
    double obj2 = xd.dot(Q * xd);
    REQUIRE(std::abs(obj - obj2) < 1e-10);
}

TEST_CASE_METHOD(QUBOTestFixture, "num_variables returns assets times bits", "[QUBOFormulation]")
{
    REQUIRE(small_->num_variables() == 4);
    REQUIRE(large_->num_variables() == 40);
}

TEST_CASE_METHOD(QUBOTestFixture, "to_json contains required schema fields", "[QUBOFormulation]")
{
    large_->build();
    auto s = large_->to_json();
    auto j = nlohmann::json::parse(s);
    REQUIRE(j.contains("schema_version"));
    REQUIRE(j.contains("Q"));
    REQUIRE(j.contains("num_assets"));
    REQUIRE(j.contains("num_bits_per_asset"));
    REQUIRE(j.contains("num_variables"));
    auto Q = j["Q"];
    REQUIRE(Q.size() == static_cast<size_t>(large_->num_variables()));
}

TEST_CASE_METHOD(QUBOTestFixture, "calling Q() before build() throws", "[QUBOFormulation]")
{
    QUBOFormulation tmp = *large_;
    // ensure unbuilt object exists
    REQUIRE_THROWS_AS(tmp.Q(), std::runtime_error);
}

TEST_CASE_METHOD(QUBOTestFixture, "invalid inputs throw std::invalid_argument", "[QUBOFormulation]")
{
    SECTION("non-square covariance")
    {
        Eigen::MatrixXd bad(2, 3);
        Eigen::VectorXd mu(2);
        mu << 0.01, 0.02;
        REQUIRE_THROWS_AS(QUBOFormulation(bad, mu, 2, 1.0, 10.0), std::invalid_argument);
    }

    SECTION("size mismatch")
    {
        Eigen::MatrixXd cov(2, 2);
        cov.setIdentity();
        Eigen::VectorXd mu(3);
        REQUIRE_THROWS_AS(QUBOFormulation(cov, mu, 2, 1.0, 10.0), std::invalid_argument);
    }

    SECTION("num_bits invalid")
    {
        Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(2,2);
        Eigen::VectorXd mu(2);
        mu << 0.01, 0.02;
        REQUIRE_THROWS_AS(QUBOFormulation(cov, mu, 0, 1.0, 10.0), std::invalid_argument);
    }

    SECTION("risk_aversion invalid")
    {
        Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(2,2);
        Eigen::VectorXd mu(2);
        mu << 0.01, 0.02;
        REQUIRE_THROWS_AS(QUBOFormulation(cov, mu, 2, -1.0, 10.0), std::invalid_argument);
    }

    SECTION("penalty_budget invalid")
    {
        Eigen::MatrixXd cov = Eigen::MatrixXd::Identity(2,2);
        Eigen::VectorXd mu(2);
        mu << 0.01, 0.02;
        REQUIRE_THROWS_AS(QUBOFormulation(cov, mu, 2, 1.0, 0.0), std::invalid_argument);
    }
}
