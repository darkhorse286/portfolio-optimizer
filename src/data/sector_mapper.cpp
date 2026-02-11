/**
 * @file sector_mapper.cpp
 * @brief Implementation of SectorMapping for asset->sector mappings
 */

#include "data/sector_mapper.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <locale>

namespace portfolio {
namespace data {

// Trim helpers (left/right)
static inline std::string ltrim(std::string s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return !std::isspace(ch); }));
    return s;
}

static inline std::string rtrim(std::string s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
    return s;
}

static inline std::string trim(std::string s)
{
    return ltrim(rtrim(std::move(s)));
}

SectorMapping::SectorMapping(const std::vector<std::string>& asset_names,
                             const std::vector<std::string>& sector_names)
    : asset_names_(asset_names)
{
    if (asset_names.size() != sector_names.size())
    {
        throw std::invalid_argument("SectorMapping: asset_names and sector_names size mismatch");
    }

    asset_to_sector_.clear();
    for (size_t i = 0; i < asset_names.size(); ++i)
    {
        const std::string a = trim(asset_names[i]);
        const std::string s = trim(sector_names[i]);
        if (a.empty())
        {
            throw std::invalid_argument("SectorMapping: empty asset name at index " + std::to_string(i));
        }
        if (s.empty())
        {
            throw std::invalid_argument("SectorMapping: empty sector name for asset '" + a + "'");
        }
        asset_to_sector_[a] = s;
    }

    build_reverse_index();
}

SectorMapping SectorMapping::from_csv(const std::string& filepath, char delimiter)
{
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        throw std::runtime_error("SectorMapping::from_csv: could not open file: " + filepath);
    }

    std::string line;
    std::vector<std::string> assets;
    std::vector<std::string> sectors;

    // Read header and ignore it (if present)
    if (!std::getline(file, line))
    {
        throw std::runtime_error("SectorMapping::from_csv: empty file: " + filepath);
    }

    // If header looks like data (no commas), we still treat first line as header per spec

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, delimiter))
        {
            fields.push_back(trim(item));
        }

        if (fields.size() < 2)
        {
            // Skip malformed lines but record as parse error if no valid rows found later
            continue;
        }

        const std::string asset = fields[0];
        const std::string sector = fields[1];

        if (asset.empty() || sector.empty())
        {
            continue; // skip empty entries
        }

        assets.push_back(asset);
        sectors.push_back(sector);
    }

    file.close();

    if (assets.empty())
    {
        throw std::runtime_error("SectorMapping::from_csv: no valid mappings found in: " + filepath);
    }

    return SectorMapping(assets, sectors);
}

SectorMapping SectorMapping::from_json(const nlohmann::json& j)
{
    if (!j.contains("assets") || !j.contains("sectors"))
    {
        throw std::invalid_argument("SectorMapping::from_json: JSON must contain 'assets' and 'sectors' arrays");
    }

    std::vector<std::string> assets = j.at("assets").get<std::vector<std::string>>();
    std::vector<std::string> sectors = j.at("sectors").get<std::vector<std::string>>();

    return SectorMapping(assets, sectors);
}

void SectorMapping::build_reverse_index()
{
    sector_to_assets_.clear();
    for (size_t i = 0; i < asset_names_.size(); ++i)
    {
        const std::string asset = asset_names_[i];
        auto it = asset_to_sector_.find(asset);
        if (it == asset_to_sector_.end())
            continue; // asset has no sector mapping

        const std::string &sector = it->second;
        sector_to_assets_[sector].push_back(static_cast<int>(i));
    }
}

std::vector<std::string> SectorMapping::get_sectors() const
{
    std::vector<std::string> sectors;
    sectors.reserve(sector_to_assets_.size());
    for (const auto &p : sector_to_assets_)
    {
        sectors.push_back(p.first);
    }
    return sectors;
}

const std::vector<int>& SectorMapping::get_assets_in_sector(const std::string& sector) const
{
    auto it = sector_to_assets_.find(sector);
    if (it == sector_to_assets_.end())
    {
        throw std::out_of_range("SectorMapping: sector not found: " + sector);
    }
    return it->second;
}

const std::string& SectorMapping::get_sector_for_asset(const std::string& asset) const
{
    auto it = asset_to_sector_.find(asset);
    if (it == asset_to_sector_.end())
    {
        throw std::out_of_range("SectorMapping: asset not found: " + asset);
    }
    return it->second;
}

std::vector<portfolio::optimizer::GroupConstraint> SectorMapping::build_group_constraints(
    double default_min, double default_max) const
{
    std::vector<portfolio::optimizer::GroupConstraint> groups;
    groups.reserve(sector_to_assets_.size());

    for (const auto &p : sector_to_assets_)
    {
        portfolio::optimizer::GroupConstraint gc;
        gc.name = p.first;
        gc.asset_indices = p.second;
        gc.min_weight = default_min;
        gc.max_weight = default_max;

        // Validate will throw on invalid configuration
        gc.validate();
        groups.push_back(std::move(gc));
    }

    return groups;
}

nlohmann::json SectorMapping::to_json() const
{
    nlohmann::json j;
    j["assets"] = asset_names_;
    std::vector<std::string> sectors(asset_names_.size());
    for (size_t i = 0; i < asset_names_.size(); ++i)
    {
        const auto it = asset_to_sector_.find(asset_names_[i]);
        sectors[i] = (it != asset_to_sector_.end()) ? it->second : std::string();
    }
    j["sectors"] = sectors;
    return j;
}

} // namespace data
} // namespace portfolio
