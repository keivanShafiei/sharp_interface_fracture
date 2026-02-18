#pragma once
#include <vector>
#include <string>
#include <array>
#include <set>
#include <unordered_map>
#include <limits>
#include <cmath>

struct Node {
    double x, y;
};

struct Element {
    std::array<int,3> vid; 
};

class Mesh {
public:
    bool load_gmsh_v2(const std::string& path);
    size_t num_nodes() const { return nodes.size(); }
    size_t num_elements() const { return elems.size(); }
    const Node& node(size_t i) const { return nodes[i]; }
    const Element& element(size_t e) const { return elems[e]; }
    int get_or_add_node(double x, double y, double tol = 1e-6);
    double compute_element_area(int v1, int v2, int v3) const;
    std::set<int> get_boundary_nodes_of_removed_region(std::array<double, 2> center, double radius);
    double average_h() const; // محاسبه میانگین اندازه المان‌ها
    double local_h_at_tip(int tip_node_id) const;
    std::array<double, 2> get_element_center(int e) const;
    double get_dist_to_tip(int e, const std::array<double, 2>& tip_pos) const;
    std::vector<int> nodes_on_y(double y_target, double tol=1e-8) const;
    int node_closest_to(double x0, double y0) const;
    double ymin() const { return _ymin; }
    double ymax() const { return _ymax; }
    void smooth_mesh_around_tip(int tip_node_id, double dx, double dy);

    
    double xmin() const;
    double xmax() const;


    std::vector<Node> nodes;
    std::vector<Element> elems; 

    std::array<double, 2> node_coords(int i) const { 
        return {nodes[i].x, nodes[i].y}; 
    }

    const std::array<int, 3>& element_nodes(int i) const { 
        return elems[i].vid; 
    }

    void remove_elements_in_radius(std::array<double, 2> center, double radius);
    int add_node(double x, double y);
    void add_element(int v1, int v2, int v3);
    int find_closest_node(double x, double y) const {
        int best = -1; double min_d = 1e15;
        for(int i=0; i<num_nodes(); ++i) {
            auto c = node_coords(i);
            double d = std::pow(c[0]-x,2) + std::pow(c[1]-y,2);
            if(d < min_d) { min_d = d; best = i; }
      }
      return best;
    }


private:
    double _ymin = std::numeric_limits<double>::max();
    double _ymax = -std::numeric_limits<double>::max();
    double _xmin = std::numeric_limits<double>::max();
    double _xmax = -std::numeric_limits<double>::max();
}; // اینجا مطمئن شوید که سمیکولن وجود دارد
