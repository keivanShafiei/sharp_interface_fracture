#pragma once

#include <string>
#include <unordered_map>

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

struct DomainBounds {
    double min_x = 0.0;
    double max_x = 0.0;
    double min_y = 0.0;
    double max_y = 0.0;
};

struct SimulationConfig {
    std::string preset_key;
    std::string mesh_path;
    DomainBounds domain;
    Vec2 initial_tip;
    double load_scale = 1.0;
    double fracture_toughness = 100.0;
    std::string results_prefix = "results";
    std::unordered_map<std::string, double> scalars;
};

const std::unordered_map<std::string, SimulationConfig>& simulation_presets();
SimulationConfig parse_simulation_config(int argc, char** argv);

