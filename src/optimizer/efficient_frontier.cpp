/**
 * @file efficient_frontier.cpp
 * @brief Implementation of efficient frontier computation
 */

#include "optimizer/efficient_frontier.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

namespace portfolio
{
    namespace optimizer
    {

        // ============================================================================
        // FrontierPoint Implementation
        // ============================================================================

        FrontierPoint::FrontierPoint()
            : expected_return(0.0),
              volatility(0.0),
              sharpe_ratio(0.0),
              is_valid(false)
        {
        }

        FrontierPoint::FrontierPoint(const OptimizationResult &result)
            : expected_return(result.expected_return),
              volatility(result.volatility),
              sharpe_ratio(result.sharpe_ratio),
              weights(result.weights),
              is_valid(result.success)
        {
        }

        // ============================================================================
        // EfficientFrontierResult Implementation
        // ============================================================================

        EfficientFrontierResult::EfficientFrontierResult()
            : success(false)
        {
        }

        bool EfficientFrontierResult::is_valid() const
        {
            return success && num_valid_points() > 0;
        }

        size_t EfficientFrontierResult::num_valid_points() const
        {
            size_t count = 0;
            for (const auto &point : points)
            {
                if (point.is_valid)
                {
                    ++count;
                }
            }
            return count;
        }

        void EfficientFrontierResult::print_summary() const
        {
            std::cout << "\n=== Efficient Frontier Summary ===\n";
            std::cout << "Status: " << (success ? "SUCCESS" : "FAILED") << "\n";
            std::cout << "Message: " << message << "\n";
            std::cout << "Total points: " << points.size() << "\n";
            std::cout << "Valid points: " << num_valid_points() << "\n";
            std::cout << std::string(60, '-') << "\n";

            if (min_variance_portfolio.is_valid)
            {
                std::cout << "\nMinimum Variance Portfolio:\n";
                std::cout << "  Expected Return:  " << std::fixed << std::setprecision(4)
                          << min_variance_portfolio.expected_return * 100 << "%\n";
                std::cout << "  Volatility:       "
                          << min_variance_portfolio.volatility * 100 << "%\n";
                std::cout << "  Sharpe Ratio:     " << std::setprecision(3)
                          << min_variance_portfolio.sharpe_ratio << "\n";
            }

            if (max_sharpe_portfolio.is_valid)
            {
                std::cout << "\nMaximum Sharpe Ratio Portfolio:\n";
                std::cout << "  Expected Return:  " << std::fixed << std::setprecision(4)
                          << max_sharpe_portfolio.expected_return * 100 << "%\n";
                std::cout << "  Volatility:       "
                          << max_sharpe_portfolio.volatility * 100 << "%\n";
                std::cout << "  Sharpe Ratio:     " << std::setprecision(3)
                          << max_sharpe_portfolio.sharpe_ratio << "\n";
            }

            if (num_valid_points() > 0)
            {
                std::cout << "\nFrontier Range:\n";

                double min_vol = std::numeric_limits<double>::max();
                double max_vol = -std::numeric_limits<double>::max();
                double min_ret = std::numeric_limits<double>::max();
                double max_ret = -std::numeric_limits<double>::max();

                for (const auto &point : points)
                {
                    if (point.is_valid)
                    {
                        min_vol = std::min(min_vol, point.volatility);
                        max_vol = std::max(max_vol, point.volatility);
                        min_ret = std::min(min_ret, point.expected_return);
                        max_ret = std::max(max_ret, point.expected_return);
                    }
                }

                std::cout << "  Return range:     " << std::fixed << std::setprecision(4)
                          << min_ret * 100 << "% to " << max_ret * 100 << "%\n";
                std::cout << "  Volatility range: "
                          << min_vol * 100 << "% to " << max_vol * 100 << "%\n";
            }

            std::cout << "================================\n"
                      << std::endl;
        }

        void EfficientFrontierResult::export_to_csv(const std::string &filepath) const
        {
            std::ofstream file(filepath);
            if (!file.is_open())
            {
                throw std::runtime_error("Could not open file for writing: " + filepath);
            }

            // Write header
            file << "return,volatility,sharpe_ratio,is_valid\n";

            // Write frontier points
            for (const auto &point : points)
            {
                file << std::fixed << std::setprecision(8)
                     << point.expected_return << ","
                     << point.volatility << ","
                     << point.sharpe_ratio << ","
                     << (point.is_valid ? "1" : "0") << "\n";
            }

            file.close();
        }

        // ============================================================================
        // EfficientFrontier Implementation
        // ============================================================================

        EfficientFrontier::EfficientFrontier()
            : num_points_(20),
              min_return_(0.0),
              max_return_(0.0),
              min_lambda_(0.1),
              max_lambda_(10.0),
              auto_range_(true)
        {
        }

        void EfficientFrontier::set_num_points(int num_points)
        {
            if (num_points < 2)
            {
                throw std::invalid_argument(
                    "Number of points must be at least 2, got: " +
                    std::to_string(num_points));
            }
            num_points_ = num_points;
        }

        void EfficientFrontier::set_min_return(double min_return)
        {
            min_return_ = min_return;
            auto_range_ = false;
        }

        void EfficientFrontier::set_max_return(double max_return)
        {
            max_return_ = max_return;
            auto_range_ = false;
        }

        void EfficientFrontier::set_risk_aversion_range(double min_lambda, double max_lambda)
        {
            if (min_lambda <= 0.0 || max_lambda <= 0.0)
            {
                throw std::invalid_argument("Risk aversion parameters must be positive");
            }
            if (min_lambda >= max_lambda)
            {
                throw std::invalid_argument("min_lambda must be less than max_lambda");
            }
            min_lambda_ = min_lambda;
            max_lambda_ = max_lambda;
        }

        void EfficientFrontier::set_auto_range(bool enable)
        {
            auto_range_ = enable;
        }

        void EfficientFrontier::validate_inputs(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints) const
        {
            OptimizerInterface::validate_inputs(expected_returns, covariance);
            constraints.validate();
        }

        EfficientFrontierResult EfficientFrontier::compute(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            double risk_free_rate) const
        {
            // Default to risk aversion method
            return compute_via_risk_aversion(
                expected_returns, covariance, constraints, risk_free_rate);
        }

        EfficientFrontierResult EfficientFrontier::compute_via_risk_aversion(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            double risk_free_rate) const
        {
            validate_inputs(expected_returns, covariance, constraints);

            EfficientFrontierResult result;
            result.points.reserve(num_points_);

            // Compute minimum variance portfolio
            result.min_variance_portfolio = compute_min_variance_portfolio(
                covariance, constraints);

            // Generate points using different risk aversion levels
            std::vector<double> lambdas;
            lambdas.reserve(num_points_);

            // Log-space sampling for better coverage
            double log_min = std::log(min_lambda_);
            double log_max = std::log(max_lambda_);
            double log_step = (log_max - log_min) / (num_points_ - 1);

            for (int i = 0; i < num_points_; ++i)
            {
                double log_lambda = log_min + i * log_step;
                lambdas.push_back(std::exp(log_lambda));
            }

            // Optimize for each risk aversion level
            MeanVarianceOptimizer optimizer(ObjectiveType::RISK_AVERSION, risk_free_rate);

            for (double lambda : lambdas)
            {
                optimizer.set_risk_aversion(lambda);

                try
                {
                    OptimizationResult opt_result = optimizer.optimize(
                        expected_returns, covariance, constraints);

                    FrontierPoint point(opt_result);
                    result.points.push_back(point);
                }
                catch (const std::exception &e)
                {
                    // Add invalid point
                    FrontierPoint invalid_point;
                    invalid_point.is_valid = false;
                    result.points.push_back(invalid_point);
                }
            }

            // Clean and sort points
            clean_frontier_points(result.points);
            sort_by_volatility(result.points);

            // Find maximum Sharpe ratio portfolio
            size_t max_sharpe_idx = find_max_sharpe_portfolio(result.points, risk_free_rate);
            if (max_sharpe_idx < result.points.size())
            {
                result.max_sharpe_portfolio = result.points[max_sharpe_idx];
            }

            // Set result status
            result.success = result.num_valid_points() > 0;
            result.message = result.success ? "Efficient frontier computed successfully" : "Failed to compute efficient frontier";

            return result;
        }

        EfficientFrontierResult EfficientFrontier::compute_via_target_return(
            const Eigen::VectorXd &expected_returns,
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints,
            double risk_free_rate) const
        {
            validate_inputs(expected_returns, covariance, constraints);

            EfficientFrontierResult result;
            result.points.reserve(num_points_);

            // Detect return range if needed
            double min_ret = min_return_;
            double max_ret = max_return_;

            if (auto_range_)
            {
                detect_return_range(expected_returns, constraints, min_ret, max_ret);
            }

            // Compute minimum variance portfolio
            result.min_variance_portfolio = compute_min_variance_portfolio(
                covariance, constraints);

            // Generate target returns
            double return_step = (max_ret - min_ret) / (num_points_ - 1);

            // Optimize for each target return
            MeanVarianceOptimizer optimizer(ObjectiveType::TARGET_RETURN, risk_free_rate);

            for (int i = 0; i < num_points_; ++i)
            {
                double target_return = min_ret + i * return_step;
                optimizer.set_target_return(target_return);

                try
                {
                    OptimizationResult opt_result = optimizer.optimize(
                        expected_returns, covariance, constraints);

                    FrontierPoint point(opt_result);
                    result.points.push_back(point);
                }
                catch (const std::exception &e)
                {
                    // Add invalid point
                    FrontierPoint invalid_point;
                    invalid_point.is_valid = false;
                    result.points.push_back(invalid_point);
                }
            }

            // Clean and sort points
            clean_frontier_points(result.points);
            sort_by_volatility(result.points);

            // Find maximum Sharpe ratio portfolio
            size_t max_sharpe_idx = find_max_sharpe_portfolio(result.points, risk_free_rate);
            if (max_sharpe_idx < result.points.size())
            {
                result.max_sharpe_portfolio = result.points[max_sharpe_idx];
            }

            // Set result status
            result.success = result.num_valid_points() > 0;
            result.message = result.success ? "Efficient frontier computed successfully" : "Failed to compute efficient frontier";

            return result;
        }

        void EfficientFrontier::detect_return_range(
            const Eigen::VectorXd &expected_returns,
            const OptimizationConstraints &constraints,
            double &min_return,
            double &max_return) const
        {
            const int n = expected_returns.size();

            // Minimum return: weighted by min_weight, biased toward lowest returns
            // Maximum return: weighted by max_weight, biased toward highest returns

            // Sort returns
            std::vector<std::pair<double, int>> sorted_returns;
            sorted_returns.reserve(n);
            for (int i = 0; i < n; ++i)
            {
                sorted_returns.push_back({expected_returns(i), i});
            }
            std::sort(sorted_returns.begin(), sorted_returns.end());

            // Minimum return: concentrate in lowest return assets
            min_return = 0.0;
            double remaining_weight = 1.0;
            for (const auto &pair : sorted_returns)
            {
                double weight = std::min(remaining_weight, constraints.max_weight);
                min_return += weight * pair.first;
                remaining_weight -= weight;
                if (remaining_weight < 1e-6)
                    break;
            }

            // Maximum return: concentrate in highest return assets
            max_return = 0.0;
            remaining_weight = 1.0;
            for (auto it = sorted_returns.rbegin(); it != sorted_returns.rend(); ++it)
            {
                double weight = std::min(remaining_weight, constraints.max_weight);
                max_return += weight * it->first;
                remaining_weight -= weight;
                if (remaining_weight < 1e-6)
                    break;
            }

            // Add some buffer
            double buffer = (max_return - min_return) * 0.05;
            min_return -= buffer;
            max_return += buffer;
        }

        FrontierPoint EfficientFrontier::compute_min_variance_portfolio(
            const Eigen::MatrixXd &covariance,
            const OptimizationConstraints &constraints) const
        {
            MeanVarianceOptimizer optimizer(ObjectiveType::MIN_VARIANCE, 0.0);

            try
            {
                // Use zero returns for min variance
                Eigen::VectorXd zero_returns = Eigen::VectorXd::Zero(covariance.rows());
                OptimizationResult result = optimizer.optimize(
                    zero_returns, covariance, constraints);

                return FrontierPoint(result);
            }
            catch (const std::exception &e)
            {
                FrontierPoint invalid;
                invalid.is_valid = false;
                return invalid;
            }
        }

        size_t EfficientFrontier::find_max_sharpe_portfolio(
            const std::vector<FrontierPoint> &points,
            double risk_free_rate) const
        {
            double max_sharpe = -std::numeric_limits<double>::max();
            size_t max_idx = 0;

            for (size_t i = 0; i < points.size(); ++i)
            {
                if (points[i].is_valid && points[i].sharpe_ratio > max_sharpe)
                {
                    max_sharpe = points[i].sharpe_ratio;
                    max_idx = i;
                }
            }

            return max_idx;
        }

        void EfficientFrontier::sort_by_volatility(std::vector<FrontierPoint> &points) const
        {
            std::sort(points.begin(), points.end(),
                      [](const FrontierPoint &a, const FrontierPoint &b)
                      {
                          if (!a.is_valid)
                              return false;
                          if (!b.is_valid)
                              return true;
                          return a.volatility < b.volatility;
                      });
        }

        void EfficientFrontier::clean_frontier_points(std::vector<FrontierPoint> &points) const
        {
            // Remove duplicate points based on volatility
            const double tolerance = 1e-6;

            std::vector<FrontierPoint> cleaned;
            cleaned.reserve(points.size());

            for (const auto &point : points)
            {
                if (!point.is_valid)
                    continue;

                bool is_duplicate = false;
                for (const auto &existing : cleaned)
                {
                    if (std::abs(point.volatility - existing.volatility) < tolerance &&
                        std::abs(point.expected_return - existing.expected_return) < tolerance)
                    {
                        is_duplicate = true;
                        break;
                    }
                }

                if (!is_duplicate)
                {
                    cleaned.push_back(point);
                }
            }

            points = cleaned;
        }

    } // namespace optimizer
} // namespace portfolio