#pragma once
#include "mesh.hpp"
#include "crack.hpp"
#include <string>

struct RemeshOptions {
    double radius = 0.0;
    double tip_radius;  // <--- این خط را اضافه کنید
    double theta_min = 22.0;
    bool verbose = false;
};

class Remesher {
public:
    explicit Remesher(RemeshOptions opts = {});
    bool remesh_local(Mesh& mesh, const CrackPath& crack, bool mesh_changed = false);
    void set_options(const RemeshOptions& opts) { opts_ = opts; }
    const RemeshOptions& options() const { return opts_; }

private:
    RemeshOptions opts_;
};

