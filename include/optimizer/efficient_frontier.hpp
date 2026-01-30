/**
 * @file efficient_frontier.hpp
 * @brief Efficient frontier computation for portfolio optimization
 *
 * Computes the Markowitz efficient frontier by solving multiple
 * portfolio optimization problems across a range of target returns
 * or risk levels.
 *
 * The efficient frontier represents the set of portfolios that offer
 * the maximum expected return for a given level of risk, or equivalently,
 * the minimum risk for a given expected return.
 *
 * Mathematical Background:
 * The efficient frontier is computed by solving:
 *
 *     For each target return r_target:
 *         Minimize: w^T * Sigma * w
 *         Subject to: mu^T * w = r_target
 *                    sum(w) = 1
 *                    w_min <= w <= w_max
 *
 * Or alternatively, for each risk aversion lambda:
 *     Minimize: w^T * Sigma * w - lambda * mu^T * w
 *
 * Performance: Approximately 100ms for 20 points with 50 assets
 */

#pragma once

#include "optimizer/optimizer_interface.hpp"
#include "optimizer/mean_variance_optimizer.hpp"
#include <vector>

namespace portfolio
{
    namespace optimizer
    {

        /**
         * @struct FrontierPoint
         * @brief Single point on the efficient frontier
         */
        struct FrontierPoint
        {
            double expected_return;     ///< Portfolio expected return
            double volatility;          ///< Portfolio volatility (std dev)
            double sharpe_ratio;        ///< Sharpe ratio
            Eigen::VectorXd weights;    ///< Portfolio weights
            bool is_valid;              ///< Point successfully computed

            /**
             * @brief Default constructor
             */
            FrontierPoint();

            /**
             * @brief Construct from optimization result
             */
            explicit FrontierPoint(const OptimizationResult& result);
        };

        /**
         * @struct EfficientFrontierResult
         * @brief Complete efficient frontier data
         */
        struct EfficientFrontierResult
        {
            std::vector<FrontierPoint> points;          ///< Frontier points
            FrontierPoint min_variance_portfolio;       ///< Minimum variance portfolio
            FrontierPoint max_sharpe_portfolio;         ///< Maximum Sharpe ratio portfolio
            bool success;                               ///< Computation succeeded
            std::string message;                        ///< Status message

            /**
             * @brief Default constructor
             */
            EfficientFrontierResult();

            /**
             * @brief Check if result is valid
             */
            bool is_valid() const;

            /**
             * @brief Get number of valid points
             */
            size_t num_valid_points() const;

            /**
             * @brief Print summary statistics
             */
            void print_summary() const;

            /**
             * @brief Export frontier data to CSV
             * @param filepath Path to output file
             */
            void export_to_csv(const std::string& filepath) const;
        };

        /**
         * @class EfficientFrontier
         * @brief Computes and manages the efficient frontier
         *
         * Generates multiple portfolio solutions along the efficient frontier
         * by varying target returns or risk aversion parameters.
         *
         * Key Features:
         * - Automatic target return range detection
         * - Multiple frontier computation methods
         * - Special point identification (min variance, max Sharpe)
         * - Data export capabilities
         * - Visualization support
         *
         * Usage Example:
         * @code
         * EfficientFrontier frontier;
         * frontier.set_num_points(20);
         * 
         * OptimizationConstraints constraints;
         * constraints.min_weight = 0.0;
         * constraints.max_weight = 0.3;
         * 
         * auto result = frontier.compute(
         *     expected_returns,
         *     covariance,
         *     constraints,
         *     0.02  // risk-free rate
         * );
         * 
         * result.print_summary();
         * result.export_to_csv("frontier.csv");
         * 
         * std::cout << "Max Sharpe portfolio:\n";
         * std::cout << "Return: " << result.max_sharpe_portfolio.expected_return << "\n";
         * std::cout << "Risk: " << result.max_sharpe_portfolio.volatility << "\n";
         * @endcode
         *
         * Thread Safety: Safe for concurrent read-only operations
         */
        class EfficientFrontier
        {
        public:
            /**
             * @brief Constructor
             */
            EfficientFrontier();

            /**
             * @brief Destructor
             */
            ~EfficientFrontier() = default;

            /**
             * @brief Compute efficient frontier
             * @param expected_returns Expected returns for each asset
             * @param covariance Covariance matrix
             * @param constraints Portfolio constraints
             * @param risk_free_rate Risk-free rate for Sharpe calculation
             * @return Efficient frontier result
             * @throws std::invalid_argument if inputs are invalid
             */
            EfficientFrontierResult compute(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                double risk_free_rate = 0.0) const;

            /**
             * @brief Compute frontier using risk aversion method
             * @param expected_returns Expected returns for each asset
             * @param covariance Covariance matrix
             * @param constraints Portfolio constraints
             * @param risk_free_rate Risk-free rate
             * @return Efficient frontier result
             *
             * Generates frontier by varying risk aversion parameter
             * from low (aggressive) to high (conservative).
             */
            EfficientFrontierResult compute_via_risk_aversion(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                double risk_free_rate = 0.0) const;

            /**
             * @brief Compute frontier using target return method
             * @param expected_returns Expected returns for each asset
             * @param covariance Covariance matrix
             * @param constraints Portfolio constraints
             * @param risk_free_rate Risk-free rate
             * @return Efficient frontier result
             *
             * Generates frontier by varying target return from
             * minimum to maximum feasible values.
             */
            EfficientFrontierResult compute_via_target_return(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints,
                double risk_free_rate = 0.0) const;

            // ===== Configuration Methods =====

            /**
             * @brief Set number of points on frontier
             * @param num_points Number of points (recommended: 10-50)
             * @throws std::invalid_argument if num_points < 2
             */
            void set_num_points(int num_points);

            /**
             * @brief Set minimum target return
             * @param min_return Minimum return (optional, auto-detected if not set)
             */
            void set_min_return(double min_return);

            /**
             * @brief Set maximum target return
             * @param max_return Maximum return (optional, auto-detected if not set)
             */
            void set_max_return(double max_return);

            /**
             * @brief Set risk aversion range
             * @param min_lambda Minimum risk aversion (aggressive)
             * @param max_lambda Maximum risk aversion (conservative)
             * @throws std::invalid_argument if range is invalid
             */
            void set_risk_aversion_range(double min_lambda, double max_lambda);

            /**
             * @brief Enable/disable auto range detection
             * @param enable true to auto-detect feasible return range
             */
            void set_auto_range(bool enable);

            // ===== Getters =====

            /**
             * @brief Get number of frontier points
             */
            int get_num_points() const { return num_points_; }

            /**
             * @brief Get minimum target return
             */
            double get_min_return() const { return min_return_; }

            /**
             * @brief Get maximum target return
             */
            double get_max_return() const { return max_return_; }

        private:
            int num_points_;                ///< Number of points on frontier
            double min_return_;             ///< Minimum target return
            double max_return_;             ///< Maximum target return
            double min_lambda_;             ///< Minimum risk aversion
            double max_lambda_;             ///< Maximum risk aversion
            bool auto_range_;               ///< Auto-detect return range

            /**
             * @brief Detect feasible return range
             * @param expected_returns Asset expected returns
             * @param constraints Portfolio constraints
             * @param min_return Output: minimum feasible return
             * @param max_return Output: maximum feasible return
             */
            void detect_return_range(
                const Eigen::VectorXd& expected_returns,
                const OptimizationConstraints& constraints,
                double& min_return,
                double& max_return) const;

            /**
             * @brief Compute minimum variance portfolio
             * @param covariance Covariance matrix
             * @param constraints Portfolio constraints
             * @return Frontier point for min variance portfolio
             */
            FrontierPoint compute_min_variance_portfolio(
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints) const;

            /**
             * @brief Find maximum Sharpe ratio portfolio
             * @param points Frontier points
             * @param risk_free_rate Risk-free rate
             * @return Index of max Sharpe portfolio
             */
            size_t find_max_sharpe_portfolio(
                const std::vector<FrontierPoint>& points,
                double risk_free_rate) const;

            /**
             * @brief Sort frontier points by volatility
             * @param points Frontier points to sort
             */
            void sort_by_volatility(std::vector<FrontierPoint>& points) const;

            /**
             * @brief Remove invalid or duplicate points
             * @param points Frontier points to clean
             */
            void clean_frontier_points(std::vector<FrontierPoint>& points) const;

            /**
             * @brief Validate frontier computation inputs
             */
            void validate_inputs(
                const Eigen::VectorXd& expected_returns,
                const Eigen::MatrixXd& covariance,
                const OptimizationConstraints& constraints) const;
        };

    } // namespace optimizer
} // namespace portfolio