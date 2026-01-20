/**
 * @file risk_model_factory.hpp
 * @brief Factory for creating risk models from configuration
 *
 * Provides a factory pattern for instantiating risk model objects based on
 * configuration parameters. This allows for configuration-driven selection
 * of risk estimation methods without modifying code.
 *
 * Supports JSON-based configuration matching the pattern used in DataLoader:
 *
 * Example configuration:
 * @code{.json}
 * {
 *   "risk_model": {
 *     "type": "ewma",
 *     "ewma_lambda": 0.94,
 *     "bias_correction": true,
 *     "shrinkage_target": "constant_correlation"
 *   }
 * }
 * @endcode
 *
 * This factory integrates seamlessly with the application's configuration
 * system, following the same patterns as the data layer.
 */

#pragma once

#include "risk/risk_model.hpp"
#include "risk/sample_covariance.hpp"
#include "risk/ewma_covariance.hpp"
#include "risk/ledoit_wolf_shrinkage.hpp"
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace portfolio
{
    namespace risk
    {

        /**
         * @struct RiskModelConfig
         * @brief Configuration parameters for risk model creation
         *
         * Holds all parameters needed to configure any supported risk model.
         * Unused parameters are ignored based on the model type.
         *
         * This follows the same pattern as DataLoaderConfig from the data layer.
         */
        struct RiskModelConfig
        {
            /**
             * @brief Type of risk model
             *
             * Supported values:
             * - "sample" or "sample_covariance": SampleCovariance
             * - "ewma" or "ewma_covariance": EWMACovariance
             * - "ledoit_wolf" or "shrinkage": LedoitWolfShrinkage
             *
             * Case-insensitive matching is used.
             */
            std::string type;

            /**
             * @brief Estimation window (number of observations)
             *
             * Used for filtering data before estimation.
             * If -1 (default), uses all available data.
             *
             * Example: estimation_window = 252 uses last 252 days
             */
            int estimation_window = -1;

            /**
             * @brief Apply Bessel's correction (n-1 vs n)
             *
             * Used by: SampleCovariance
             * Default: true (unbiased estimate)
             */
            bool bias_correction = true;

            /**
             * @brief EWMA decay parameter
             *
             * Used by: EWMACovariance
             * Valid range: (0, 1)
             * Default: 0.94 (RiskMetrics standard)
             *
             * Common values:
             * - 0.94: Standard for daily equity data
             * - 0.97: More stable, slower adaptation
             * - 0.99: Very stable
             */
            double ewma_lambda = 0.94;

            /**
             * @brief Shrinkage target type
             *
             * Used by: LedoitWolfShrinkage
             * Valid values:
             * - "identity": Identity matrix target
             * - "constant_variance": Constant variance target
             * - "constant_correlation": Constant correlation (default)
             * - "market_model": Single-factor model
             */
            std::string shrinkage_target = "constant_correlation";

            /**
             * @brief Fixed shrinkage intensity override
             *
             * Used by: LedoitWolfShrinkage
             * Valid range: [0, 1] or -1 for auto
             * Default: -1.0 (compute optimal shrinkage)
             *
             * Set to a value in [0, 1] to use fixed shrinkage.
             * Useful for testing or when you have prior knowledge.
             */
            double shrinkage_intensity = -1.0;

            /**
             * @brief Create configuration from JSON
             * @param j JSON object containing risk model configuration
             * @return RiskModelConfig structure
             * @throws nlohmann::json::exception if JSON is malformed
             *
             * Provides defaults for missing fields.
             *
             * Example JSON:
             * @code{.json}
             * {
             *   "type": "ewma",
             *   "ewma_lambda": 0.94
             * }
             * @endcode
             */
            static RiskModelConfig from_json(const nlohmann::json &j);

            /**
             * @brief Convert configuration to JSON
             * @return JSON object with all configuration parameters
             *
             * Useful for saving configurations or debugging.
             */
            nlohmann::json to_json() const;
        };

        /**
         * @class RiskModelFactory
         * @brief Factory for creating risk model instances
         *
         * Provides static methods to create risk model objects based on
         * configuration parameters. Follows the factory pattern to decouple
         * risk model selection from application logic.
         *
         * Usage Pattern:
         * @code
         * // From JSON configuration
         * auto config = RiskModelConfig::from_json(json_obj);
         * auto model = RiskModelFactory::create(config);
         *
         * // Direct creation
         * auto model = RiskModelFactory::create("ewma", 0.94);
         *
         * // Use the model
         * auto cov = model->estimate_covariance(returns);
         * @endcode
         *
         * Benefits:
         * - Configuration-driven risk model selection
         * - Type-safe parameter handling
         * - Centralized validation
         * - Easy to extend with new risk models
         */
        class RiskModelFactory
        {
        public:
            /**
             * @brief Create risk model from configuration
             * @param config Configuration structure
             * @return Unique pointer to created risk model
             * @throws std::invalid_argument if type is unknown
             * @throws std::invalid_argument if parameters are invalid
             *
             * Creates the appropriate risk model based on config.type and
             * initializes it with the specified parameters.
             *
             * Example:
             * @code
             * RiskModelConfig config;
             * config.type = "ewma";
             * config.ewma_lambda = 0.94;
             * auto model = RiskModelFactory::create(config);
             * @endcode
             */
            static std::unique_ptr<RiskModel> create(const RiskModelConfig &config);

            /**
             * @brief Create risk model from type string and JSON
             * @param type Risk model type
             * @param params JSON object with parameters
             * @return Unique pointer to created risk model
             * @throws std::invalid_argument if type is unknown
             *
             * Convenience method for creating from JSON configuration.
             *
             * Example:
             * @code
             * nlohmann::json params = {{"ewma_lambda", 0.94}};
             * auto model = RiskModelFactory::create("ewma", params);
             * @endcode
             */
            static std::unique_ptr<RiskModel> create(
                const std::string &type,
                const nlohmann::json &params);

            /**
             * @brief Create sample covariance estimator
             * @param bias_correction Apply Bessel's correction
             * @return Unique pointer to SampleCovariance
             *
             * Convenience method for direct creation.
             */
            static std::unique_ptr<RiskModel> create_sample_covariance(
                bool bias_correction = true);

            /**
             * @brief Create EWMA covariance estimator
             * @param lambda Decay parameter
             * @return Unique pointer to EWMACovariance
             * @throws std::invalid_argument if lambda not in (0, 1)
             *
             * Convenience method for direct creation.
             */
            static std::unique_ptr<RiskModel> create_ewma_covariance(
                double lambda = 0.94);

            /**
             * @brief Create Ledoit-Wolf shrinkage estimator
             * @param target Shrinkage target type
             * @param shrinkage_override Fixed shrinkage (-1 for auto)
             * @return Unique pointer to LedoitWolfShrinkage
             *
             * Convenience method for direct creation.
             */
            static std::unique_ptr<RiskModel> create_ledoit_wolf(
                ShrinkageTarget target = ShrinkageTarget::CONSTANT_CORRELATION,
                double shrinkage_override = -1.0);

            /**
             * @brief Get list of supported risk model types
             * @return Vector of supported type strings
             *
             * Useful for validation and help messages.
             */
            static std::vector<std::string> get_supported_types();

        private:
            /**
             * @brief Convert string to shrinkage target enum
             * @param target_str String representation
             * @return ShrinkageTarget enum value
             * @throws std::invalid_argument if string is invalid
             *
             * Supported strings (case-insensitive):
             * - "identity"
             * - "constant_variance"
             * - "constant_correlation"
             * - "market_model"
             */
            static ShrinkageTarget parse_shrinkage_target(const std::string &target_str);

            /**
             * @brief Normalize type string (lowercase, trim)
             * @param type Input type string
             * @return Normalized type string
             */
            static std::string normalize_type(const std::string &type);
        };

    } // namespace risk
} // namespace portfolio