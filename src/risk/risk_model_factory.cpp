/**
 * @file risk_model_factory.cpp
 * @brief Implementation of risk model factory
 */

#include "risk/risk_model_factory.hpp"
#include <algorithm>
#include <cctype>

namespace portfolio
{
    namespace risk
    {

        // RiskmodelConfig implementation
        RiskModelConfig RiskModelConfig::from_json(const nlohmann::json &doc)
        {
            RiskModelConfig config;

            // Type (required)
            if (!doc.contains("type") || !doc["type"].is_string())
            {
                throw std::invalid_argument("Risk model configuration must specify 'type'");
            }

            config.type = doc["type"].get<std::string>();

            // Estimation window (optional)
            if (doc.contains("estimation_window"))
            {
                config.estimation_window = doc["estimation_window"].get<int>();
            }

            // Bias correction (optional)
            if (doc.contains("bias_correction"))
            {
                config.bias_correction = doc["bias_correction"].get<bool>();
            }

            // EWMA lambda (optional)
            if (doc.contains("ewma_lambda"))
            {
                config.ewma_lambda = doc["ewma_lambda"].get<double>();
            }

            // Shrinkage target (optional)
            if (doc.contains("shrinkage_target"))
            {
                config.shrinkage_target = doc["shrinkage_target"].get<std::string>();
            }

            // Shrinkage intensity (optional)
            if (doc.contains("shrinkage_intensity"))
            {
                config.shrinkage_intensity = doc["shrinkage_intensity"].get<double>();
            }

            return config;
        }

        nlohmann::json RiskModelConfig::to_json() const
        {
            nlohmann::json{
                {"type", type},
                {"estimation_window", estimation_window},
                {"bias_correction", bias_correction},
                {"ewma_lambda", ewma_lambda},
                {"shrinkage_target", shrinkage_target},
                {"shrinkage_intensity", shrinkage_intensity}};
        }

        // RiskModelFactory implementation
        std::string RiskModelFactory::normalize_type(const std::string &type)
        {
            std::string normalized = type;

            // Convert to lowercase
            std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char c)
                           { return std::tolower(c); });
            return normalized;
        }

        ShrinkageTarget RiskModelFactory::parse_shrinkage_target(const std::string &target_str)
        {
            std::string normalized = normalize_type(target_str);

            if (normalized == "identity")
            {
                return ShrinkageTarget::IDENTITY;
            }
            else if (normalized == "constant_variance")
            {
                return ShrinkageTarget::CONSTANT_VARIANCE;
            }
            else if (normalized == "constant_correlation")
            {
                return ShrinkageTarget::CONSTANT_CORRELATION;
            }
            else if (normalized == "market_model")
            {
                return ShrinkageTarget::MARKET_MODEL;
            }
            else
            {
                throw std::invalid_argument(
                    "Unknown shrinkage target: '" + target_str + "'. Valid options: identity, constant_variance, constant_correlation, market_model");
            }
        }

        std::unique_ptr<RiskModel> RiskModelFactory::create(const RiskModelConfig &config)
        {
            std::string type = normalize_type(config.type);

            if (type == "sample" || type == "sample_covariance")
            {
                return std::make_unique<SampleCovariance>(config.bias_correction);
            }
            else if (type == "ewma" || type == "ewma_covariance")
            {
                return std::make_unique<EWMACovariance>(config.ewma_lambda);
            }
            else if (type == "ledoit_wolf" || type == "shrinkage")
            {
                ShrinkageTarget target = parse_shrinkage_target(config.shrinkage_target);
                return std::make_unique<LedoitWolfShrinkage>(target, config.shrinkage_intensity);
            }
            else
            {
                throw std::invalid_argument(
                    "Unknown risk model type: '" + config.type + "'. "
                                                                 "Valid options: sample, ewma, ledoit_wolf");
            }
        }

        std::unique_ptr<RiskModel> RiskModelFactory::create(const std::string &type, const nlohmann::json &params)
        {
            // Create config from JSON
            nlohmann::json config_json = params;
            config_json["type"] = type;

            RiskModelConfig config = RiskModelConfig::from_json(config_json);
            return create(config);
        }

        std::unique_ptr<RiskModel> RiskModelFactory::create_sample_covariance(bool bias_correction)
        {
            return std::make_unique<SampleCovariance>(bias_correction);
        }

        std::unique_ptr<RiskModel> RiskModelFactory::create_ewma_covariance(double lambda)
        {
            return std::make_unique<EWMACovariance>(lambda);
        }

        std::unique_ptr<RiskModel> RiskModelFactory::create_ledoit_wolf(ShrinkageTarget target, double shrinkage_override)
        {
            return std::make_unique<LedoitWolfShrinkage>(target, shrinkage_override);
        }

        std::vector<std::string> RiskModelFactory::get_supported_types()
        {
            return {
                "sample",
                "sample_covariance",
                "ewma",
                "ewma_covariance",
                "ledoit_wolf",
                "shrinkage"};
        }

    } // namespace risk
} // namespace portfolio