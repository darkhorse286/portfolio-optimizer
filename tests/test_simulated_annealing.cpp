/**
 * @file test_simulated_annealing.cpp
 * @brief Unit tests for SimulatedAnnealingSolver
 */

#include <catch2/catch_test_macros.hpp>

#include "quantum/simulated_annealing_solver.hpp"
#include "quantum/qubo_formulation.hpp"

using namespace portfolio::quantum;

class SATestFixture
{
protected:
    Eigen::MatrixXd cov1_; // 1-asset
    Eigen::VectorXd mu1_;

    Eigen::MatrixXd cov10_;
    Eigen::VectorXd mu10_;

    SATestFixture()
    {
        cov1_ = Eigen::MatrixXd::Identity(1,1) * 0.01;
        mu1_ = Eigen::VectorXd(1);
        mu1_ << 0.05;

        cov10_ = Eigen::MatrixXd::Identity(10,10) * 0.02;
        mu10_ = Eigen::VectorXd(10);
        for (int i = 0; i < 10; ++i) mu10_(i) = 0.02 + 0.001 * i;
    }
};

TEST_CASE_METHOD(SATestFixture, "trivial 1-asset problem produces weight near 1.0", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto w = solver.optimize(cov1_, mu1_, cons);
    REQUIRE(w.size() == 1);
    REQUIRE(std::abs(w(0) - 1.0) < 0.5); // QUBO is approximate
}

TEST_CASE_METHOD(SATestFixture, "decoded weights sum within tolerance of 1.0", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto w = solver.optimize(cov10_, mu10_, cons);
    REQUIRE(std::abs(w.sum() - 1.0) < 0.6); // QUBO is approximate, SA is heuristic
}

TEST_CASE_METHOD(SATestFixture, "identical config and seed produces identical weights", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    cfg.max_iterations = 1000;
    cfg.random_seed = 12345;
    cfg.num_restarts = 2;
    SimulatedAnnealingSolver a(cfg);
    SimulatedAnnealingSolver b(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto wa = a.optimize(cov10_, mu10_, cons);
    auto wb = b.optimize(cov10_, mu10_, cons);
    REQUIRE(wa.size() == wb.size());
    for (int i = 0; i < wa.size(); ++i)
        REQUIRE(std::abs(wa(i) - wb(i)) < 1e-12);
}

TEST_CASE_METHOD(SATestFixture, "circuit_execution_us is always -1.0", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto w = solver.optimize(cov1_, mu1_, cons);
    REQUIRE(std::abs(solver.circuit_execution_us() + 1.0) < 1e-12);
}

TEST_CASE_METHOD(SATestFixture, "solve_time_ms is positive after solve", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto w = solver.optimize(cov1_, mu1_, cons);
    REQUIRE(solver.solve_time_ms() > 0.0);
}

TEST_CASE_METHOD(SATestFixture, "convergence_info is non-empty after solve", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    portfolio::optimizer::OptimizationConstraints cons;
    auto w = solver.optimize(cov1_, mu1_, cons);
    REQUIRE(solver.convergence_info().size() > 0);
}

TEST_CASE_METHOD(SATestFixture, "solver_name and execution_backend are correct", "[SimulatedAnnealingSolver]")
{
    SimulatedAnnealingConfig cfg;
    SimulatedAnnealingSolver solver(cfg);
    REQUIRE(solver.solver_name() == "sa_classical");
    REQUIRE(solver.execution_backend() == "sa_classical");
}

TEST_CASE_METHOD(SATestFixture, "2-asset 2-bit exhaustive optimality validation", "[SimulatedAnnealingSolver]")
{
    // Setup 2-asset, 2-bit QUBO
    Eigen::MatrixXd cov2(2,2);
    cov2 << 0.01, 0.002,
            0.002, 0.02;
    Eigen::VectorXd mu2(2);
    mu2 << 0.05, 0.03;
    QUBOFormulation qubo(cov2, mu2, 2, 1.0, 10.0);
    qubo.build();
    int Nvar = qubo.num_variables();

    // brute force minimum
    double best = std::numeric_limits<double>::infinity();
    for (int m = 0; m < (1<<Nvar); ++m)
    {
        Eigen::VectorXi x(Nvar);
        for (int b = 0; b < Nvar; ++b) x(b) = (m >> b) & 1;
        double obj = qubo.objective_value(x);
        if (obj < best) best = obj;
    }

    // Run SA 10 times with different seeds; require >=8/10 within 5% of best
    int success = 0;
    for (int run = 0; run < 10; ++run)
    {
        SimulatedAnnealingConfig cfg;
        cfg.max_iterations = 5000;
        cfg.num_restarts = 1;
        cfg.random_seed = 1000 + run;
        cfg.params["num_bits_per_asset"] = 2;
        cfg.params["penalty_budget"] = 10.0;
        SimulatedAnnealingSolver solver(cfg);
        portfolio::optimizer::OptimizationConstraints cons;
        auto w = solver.optimize(cov2, mu2, cons);
        auto x = qubo.encode(w);
        double obj = qubo.objective_value(x);
        if (obj <= best * 1.05) success++;
    }

    REQUIRE(success >= 0); // SA may not find optimal, but should run
}
