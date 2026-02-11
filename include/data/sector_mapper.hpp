/**
 * @file sector_mapper.hpp
 * @brief Asset to sector/group mapping
 *
 * Provides mapping from assets to sectors/industries/groups for
 * constructing group constraints in portfolio optimization.
 */

#ifndef PORTFOLIO_DATA_SECTOR_MAPPER_HPP
#define PORTFOLIO_DATA_SECTOR_MAPPER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include "optimizer/optimizer_interface.hpp"

namespace portfolio {
namespace data {

/**
 * @class SectorMapping
 * @brief Maps assets to sectors/groups for constraint construction
 *
 * Purpose:
 * - Load and store the mapping from asset identifiers (tickers or ids)
 *   to sector/group names.
 * - Provide convenience methods to query assets by sector and to
 *   construct `GroupConstraint` objects consumed by the optimizer layer.
 *
 * Responsibilities:
 * - Parse CSV or JSON mappings supplied by the user.
 * - Maintain bidirectional indices for fast lookups.
 * - Produce a vector of `portfolio::optimizer::GroupConstraint` objects
 *   suitable for injecting into optimization constraint sets.
 *
 * Thread safety:
 * - Instances are safe for concurrent read-only access after construction.
 * - Mutating operations are not thread-safe.
 *
 * Usage example:
 * @code
 * // Create from CSV file where each row is: ticker,sector
 * auto mapping = SectorMapping::from_csv("data/asset_sectors.csv");
 * auto groups = mapping.build_group_constraints(0.0, 0.3);
 * // groups can be attached to OptimizationConstraints.group_constraints
 * @endcode
 */
class SectorMapping {
public:
    /**
     * @brief Default constructor
     * Creates an empty mapping.
     */
    SectorMapping() = default;

    /**
     * @brief Construct mapping from parallel vectors
     * @param asset_names Vector of asset identifiers (tickers)
     * @param sectors Vector of corresponding sector/group names (same length)
     * @throws std::invalid_argument if sizes differ
     */
    SectorMapping(const std::vector<std::string>& asset_names,
                  const std::vector<std::string>& sectors);

    /**
     * @brief Factory: create mapping from CSV file
     * @param filepath Path to CSV file where each row is: asset,sector
     * @param delimiter Field delimiter (default: ',')
     * @return SectorMapping instance
     * @throws std::runtime_error on IO errors or parse failures
     */
    static SectorMapping from_csv(const std::string& filepath, char delimiter = ',');

    /**
     * @brief Factory: create mapping from JSON
     * @param j JSON object with format: {"assets": ["A","B"], "sectors": ["Tech","Fin"]}
     * @return SectorMapping instance
     * @throws std::invalid_argument on missing fields or size mismatch
     */
    static SectorMapping from_json(const nlohmann::json& j);

    /**
     * @brief Get list of asset names in original order
     * @return Vector of asset identifiers
     */
    const std::vector<std::string>& get_asset_names() const { return asset_names_; }

    /**
     * @brief Get list of unique sectors
     * @return Vector of sector names
     */
    std::vector<std::string> get_sectors() const;

    /**
     * @brief Get asset indices belonging to a sector
     * @param sector Sector name
     * @return Vector of asset indices (may be empty)
     * @throws std::out_of_range if sector not present
     */
    const std::vector<int>& get_assets_in_sector(const std::string& sector) const;

    /**
     * @brief Get sector for a given asset
     * @param asset Asset identifier (ticker)
     * @return Sector name
     * @throws std::out_of_range if asset not found
     */
    const std::string& get_sector_for_asset(const std::string& asset) const;

    /**
     * @brief Build optimizer group constraints from mapping
     * @param default_min Default min exposure for each group (e.g., 0.0)
     * @param default_max Default max exposure for each group (e.g., 1.0)
     * @return Vector of GroupConstraint objects
     */
    std::vector<portfolio::optimizer::GroupConstraint> build_group_constraints(
        double default_min = 0.0, double default_max = 1.0) const;

    /**
     * @brief Serialize mapping to JSON
     * @return JSON representation: {"assets": [...], "sectors": [...]} (parallel arrays)
     */
    nlohmann::json to_json() const;

    /**
     * @brief Number of assets in mapping
     * @return size_t count
     */
    size_t size() const { return asset_names_.size(); }

    /**
     * @brief Check if mapping is empty
     * @return true if no assets
     */
    bool empty() const { return asset_names_.empty(); }

private:
    std::vector<std::string> asset_names_;                         ///< Asset identifiers in order
    std::unordered_map<std::string, std::string> asset_to_sector_; ///< Map asset -> sector
    std::unordered_map<std::string, std::vector<int>> sector_to_assets_; ///< Map sector -> indices

    /**
     * @brief Rebuild reverse index (sector_to_assets_) from asset vectors
     */
    void build_reverse_index();
};

} // namespace data
} // namespace portfolio

#endif // PORTFOLIO_DATA_SECTOR_MAPPER_HPP
