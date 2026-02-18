#include "simulation_config.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {

Vec2 parse_vec2(const std::string& token)
{
    std::stringstream ss(token);
    std::string item;
    Vec2 v{};
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse x from Vec2 token: " + token);
    }
    v.x = std::stod(item);
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse y from Vec2 token: " + token);
    }
    v.y = std::stod(item);
    return v;
}

DomainBounds parse_domain(const std::string& token)
{
    std::stringstream ss(token);
    std::string item;
    DomainBounds b{};
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse min_x for domain: " + token);
    }
    b.min_x = std::stod(item);
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse max_x for domain: " + token);
    }
    b.max_x = std::stod(item);
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse min_y for domain: " + token);
    }
    b.min_y = std::stod(item);
    if (!std::getline(ss, item, ',')) {
        throw std::runtime_error("Failed to parse max_y for domain: " + token);
    }
    b.max_y = std::stod(item);
    return b;
}

std::pair<std::string, double> parse_param(const std::string& token)
{
    auto eq = token.find('=');
    if (eq == std::string::npos || eq == 0 || eq == token.size() - 1) {
        throw std::runtime_error("Expected --param name=value, got: " + token);
    }
    std::string key = token.substr(0, eq);
    std::string value_str = token.substr(eq + 1);
    double value = std::stod(value_str);
    return {key, value};
}

std::string lower_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

void print_usage(const char* binary)
{
    std::cerr << "Usage: " << binary
              << " [--case {tpb|sent}] [--mesh path] [--initial-tip x,y]\n"
              << "             [--domain xmin,xmax,ymin,ymax] [--Gc value] [--load-scale value]\n"
              << "             [--results-prefix dir] [--param name=value] ...\n";
}

const std::unordered_map<std::string, SimulationConfig>& build_presets()
{
    static const std::unordered_map<std::string, SimulationConfig> presets = {
        {
            "tpb",
            SimulationConfig{
                "tpb",
                "tpb_rect_coarse.msh",
                DomainBounds{0.0, 0.203, 0.0, 0.102},
                Vec2{0.075, 0.051},
                1.0,
                100.0,
                "results/tpb",
                {{"support_span", 0.16}, {"roller_radius", 0.005}}
            }
        },
        {
            "sent",
            SimulationConfig{
                "sent",
                "sent_rect_coarse.msh",
                DomainBounds{0.0, 0.305, 0.0, 0.076},
                Vec2{0.0, 0.038},
                1.0,
                100.0,
                "results/sent",
                {{"applied_displacement", 0.5}}
            }
        }
    };
    return presets;
}

} // namespace

const std::unordered_map<std::string, SimulationConfig>& simulation_presets()
{
    return build_presets();
}

SimulationConfig parse_simulation_config(int argc, char** argv)
{
    std::string requested_case = "tpb";
    std::optional<std::string> mesh_override;
    std::optional<Vec2> init_tip_override;
    std::optional<DomainBounds> domain_override;
    std::optional<double> load_scale_override;
    std::optional<double> gc_override;
    std::optional<std::string> results_prefix_override;
    std::unordered_map<std::string, double> scalar_overrides;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--case" && i + 1 < argc) {
            requested_case = lower_copy(argv[++i]);
        } else if (arg == "--mesh" && i + 1 < argc) {
            mesh_override = argv[++i];
        } else if (arg == "--initial-tip" && i + 1 < argc) {
            init_tip_override = parse_vec2(argv[++i]);
        } else if (arg == "--domain" && i + 1 < argc) {
            domain_override = parse_domain(argv[++i]);
        } else if (arg == "--load-scale" && i + 1 < argc) {
            load_scale_override = std::stod(argv[++i]);
        } else if (arg == "--Gc" && i + 1 < argc) {
            gc_override = std::stod(argv[++i]);
        } else if (arg == "--results-prefix" && i + 1 < argc) {
            results_prefix_override = argv[++i];
        } else if (arg == "--param" && i + 1 < argc) {
            auto kv = parse_param(argv[++i]);
            scalar_overrides.emplace(std::move(kv));
        } else if (arg == "--help") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            print_usage(argv[0]);
            throw std::runtime_error("Unknown CLI argument: " + arg);
        }
    }

    const auto& presets = simulation_presets();
    auto it = presets.find(requested_case);
    if (it == presets.end()) {
        std::stringstream ss;
        ss << "Unknown case preset '" << requested_case << "'. Allowed keys:";
        for (const auto& kv : presets) {
            ss << ' ' << kv.first;
        }
        throw std::runtime_error(ss.str());
    }

    SimulationConfig config = it->second;

    if (mesh_override) {
        config.mesh_path = *mesh_override;
    }
    if (init_tip_override) {
        config.initial_tip = *init_tip_override;
    }
    if (domain_override) {
        config.domain = *domain_override;
    }
    if (load_scale_override) {
        config.load_scale = *load_scale_override;
    }
    if (gc_override) {
        config.fracture_toughness = *gc_override;
    }
    if (results_prefix_override) {
        config.results_prefix = *results_prefix_override;
    }
    config.scalars.insert(scalar_overrides.begin(), scalar_overrides.end());

    return config;
}

