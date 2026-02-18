#include "mesh.hpp"
#include <fstream>
#include <sstream>
#include <set>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>


bool Mesh::load_gmsh_v2(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Cannot open Gmsh file: " << path << std::endl;
        return false;
    }

    std::string line;
    double version = 2.2;
    std::map<int, int> gmsh_to_internal;
    nodes.clear();
    elems.clear();

    _xmin = _ymin = std::numeric_limits<double>::max();
    _xmax = _ymax = -std::numeric_limits<double>::max();

    while (std::getline(in, line)) {
        if (line == "$MeshFormat") {
            std::getline(in, line);
            std::istringstream iss(line);
            iss >> version;
        } 
        else if (line == "$Nodes") {
            if (version < 4.0) {
                // منطق قدیمی نسخه 2.2
                std::getline(in, line);
                int num_nodes_file = std::stoi(line);
                for (int i = 0; i < num_nodes_file; ++i) {
                    std::getline(in, line);
                    std::istringstream iss(line);
                    int id; double x, y, z;
                    iss >> id >> x >> y >> z;
                    nodes.push_back({x, y});
                    gmsh_to_internal[id] = (int)nodes.size() - 1;
                    _xmin = std::min(_xmin, x); _xmax = std::max(_xmax, x);
                    _ymin = std::min(_ymin, y); _ymax = std::max(_ymax, y);
                }
            } else {
                // منطق جدید نسخه 4.1 (مطابق فایل جدید شما)
                size_t numEntityBlocks, totalNodes, minTag, maxTag;
                in >> numEntityBlocks >> totalNodes >> minTag >> maxTag;
                nodes.reserve(totalNodes);
                
                for (size_t b = 0; b < numEntityBlocks; ++b) {
                    int entityDim, entityTag, parametric;
                    size_t numNodesInBlock;
                    in >> entityDim >> entityTag >> parametric >> numNodesInBlock;
                    
                    std::vector<int> tags(numNodesInBlock);
                    for (size_t i = 0; i < numNodesInBlock; ++i) in >> tags[i];
                    
                    for (size_t i = 0; i < numNodesInBlock; ++i) {
                        double x, y, z;
                        in >> x >> y >> z;
                        nodes.push_back({x, y});
                        gmsh_to_internal[tags[i]] = (int)nodes.size() - 1;
                        _xmin = std::min(_xmin, x); _xmax = std::max(_xmax, x);
                        _ymin = std::min(_ymin, y); _ymax = std::max(_ymax, y);
                    }
                }
            }
        } 
        else if (line == "$Elements") {
            if (version < 4.0) {
                // نسخه 2.2
                std::getline(in, line);
                int num_elems_file = std::stoi(line);
                for (int i = 0; i < num_elems_file; ++i) {
                    std::getline(in, line);
                    std::istringstream iss(line);
                    int id, type, ntags;
                    iss >> id >> type >> ntags;
                    for (int t = 0; t < ntags; ++t) { int dummy; iss >> dummy; }
                    if (type == 2) { // Triangle
                        int n1, n2, n3;
                        if (iss >> n1 >> n2 >> n3) {
                            if (gmsh_to_internal.count(n1) && gmsh_to_internal.count(n2) && gmsh_to_internal.count(n3)) {
                                Element e;
                                e.vid = {gmsh_to_internal[n1], gmsh_to_internal[n2], gmsh_to_internal[n3]};
                                elems.push_back(e);
                            }
                        }
                    }
                }
            } else {
                // نسخه 4.1
                size_t numEntityBlocks, totalElems, minTag, maxTag;
                in >> numEntityBlocks >> totalElems >> minTag >> maxTag;
                for (size_t b = 0; b < numEntityBlocks; ++b) {
                    int entityDim, entityTag, elementType;
                    size_t numElemsInBlock;
                    in >> entityDim >> entityTag >> elementType >> numElemsInBlock;
                    
                    for (size_t i = 0; i < numElemsInBlock; ++i) {
                        int id;
                        in >> id;
                        if (elementType == 2) { // Triangle
                            int n1, n2, n3;
                            in >> n1 >> n2 >> n3;
                            if (gmsh_to_internal.count(n1) && gmsh_to_internal.count(n2) && gmsh_to_internal.count(n3)) {
                                Element e;
                                e.vid = {gmsh_to_internal[n1], gmsh_to_internal[n2], gmsh_to_internal[n3]};
                                
                                // Validate node indices
                                bool valid = true;
                                for (int k = 0; k < 3; ++k) {
                                    if (e.vid[k] < 0 || e.vid[k] >= (int)nodes.size()) {
                                        std::cerr << "[ERROR] Element " << id << " node " << k 
                                                  << " has invalid index " << e.vid[k] 
                                                  << " (nodes.size=" << nodes.size() << ")\n";
                                        valid = false;
                                    }
                                }
                                
                                if (valid) {
                                    elems.push_back(e);
                                } else {
                                    std::cerr << "[ERROR] Skipping invalid element " << id << "\n";
                                }
                            } else {
                                std::cerr << "[WARNING] Element " << id << " references unknown Gmsh node IDs: "
                                          << n1 << "," << n2 << "," << n3 << "\n";
                            }
                        } else {
                            // Skip non-triangular elements (points, lines, quads, etc.)
                            // Must consume correct number of node IDs from stream
                            int numNodes = 0;
                            switch(elementType) {
                                case 1:  numNodes = 2; break;  // 2-node line
                                case 3:  numNodes = 4; break;  // 4-node quadrangle
                                case 4:  numNodes = 4; break;  // 4-node tetrahedron  
                                case 5:  numNodes = 8; break;  // 8-node hexahedron
                                case 8:  numNodes = 3; break;  // 3-node line
                                case 9:  numNodes = 6; break;  // 6-node triangle
                                case 15: numNodes = 1; break;  // 1-node point
                                case 16: numNodes = 8; break;  // 8-node quadrangle
                                case 20: numNodes = 9; break;  // 9-node triangle
                                case 21: numNodes = 10; break; // 10-node triangle
                                default:
                                    std::cerr << "[WARNING] Unknown element type " << elementType 
                                              << " in element " << id << ", skipping\n";
                                    // Try to skip to next element by reading until newline
                                    std::string remainder;
                                    std::getline(in, remainder);
                                    continue;
                            }
                            // Consume node IDs to keep stream synchronized
                            for(int n=0; n<numNodes; ++n) { 
                                int dummy; 
                                in >> dummy; 
                            }
                        }
                    }
                }
            }
        }
    }
    
    // اعتبارسنجی نهایی
    std::cout << "Successfully loaded Mesh V" << version << " | Nodes: " << nodes.size() << " | Elements: " << elems.size() << std::endl;
    return (!nodes.empty() && !elems.empty());
}

std::vector<int> Mesh::nodes_on_y(double y_target, double tol) const {
    std::vector<int> ids;
    for (int i=0; i<(int)nodes.size(); ++i) {
        if (std::fabs(nodes[i].y - y_target) <= tol) ids.push_back(i);
    }
    return ids;
}

// Change nodes[i][0] -> nodes[i].x
// Change nodes[i][1] -> nodes[i].y
int Mesh::node_closest_to(double x, double y) const {
    int best = -1; 
    double min_d = std::numeric_limits<double>::max();
    for(int i = 0; i < (int)nodes.size(); ++i) {
        double dx = nodes[i].x - x; 
        double dy = nodes[i].y - y;
        double d = dx*dx + dy*dy;
        if(d < min_d) { min_d = d; best = i; }
    }
    return best;
}

double Mesh::xmin() const { return _xmin; }
double Mesh::xmax() const { return _xmax; }
void Mesh::remove_elements_in_radius(std::array<double, 2> center, double radius) {
    elems.erase(std::remove_if(elems.begin(), elems.end(), [&](const Element& e) {
        // اگر مرکز هندسی المان در شعاع بازسازی بود، آن را حذف کن
        double cx = (nodes[e.vid[0]].x + nodes[e.vid[1]].x + nodes[e.vid[2]].x) / 3.0;
        double cy = (nodes[e.vid[0]].y + nodes[e.vid[1]].y + nodes[e.vid[2]].y) / 3.0;
        double dist = std::sqrt(std::pow(cx - center[0], 2) + std::pow(cy - center[1], 2));
        return dist < radius;
    }), elems.end());
}

int Mesh::add_node(double x, double y) {
    nodes.push_back({x, y});
    // آپدیت کردن محدوده‌ها با اضافه شدن نود جدید
    _xmin = std::min(_xmin, x);
    _xmax = std::max(_xmax, x);
    _ymin = std::min(_ymin, y);
    _ymax = std::max(_ymax, y);
    return nodes.size() - 1;
}
void Mesh::add_element(int v1, int v2, int v3) {
    Element e;
    e.vid = {v1, v2, v3};
    elems.push_back(e);
}

int Mesh::get_or_add_node(double x, double y, double tol) {
    for (int i = 0; i < (int)nodes.size(); ++i) {
        if (std::sqrt(std::pow(nodes[i].x - x, 2) + std::pow(nodes[i].y - y, 2)) < tol) {
            return i;
        }
    }
    return add_node(x, y);
}

double Mesh::compute_element_area(int v1, int v2, int v3) const {
    double x1 = nodes[v1].x, y1 = nodes[v1].y;
    double x2 = nodes[v2].x, y2 = nodes[v2].y;
    double x3 = nodes[v3].x, y3 = nodes[v3].y;
    return 0.5 * std::abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
}

std::set<int> Mesh::get_boundary_nodes_of_removed_region(std::array<double, 2> center, double radius) {
    std::set<int> b_nodes;
    for (const auto& e : elems) {
        double cx = (nodes[e.vid[0]].x + nodes[e.vid[1]].x + nodes[e.vid[2]].x) / 3.0;
        double cy = (nodes[e.vid[0]].y + nodes[e.vid[1]].y + nodes[e.vid[2]].y) / 3.0;
        double dist = std::sqrt(std::pow(cx - center[0], 2) + std::pow(cy - center[1], 2));
        
        if (dist < radius) {
            for (int id : e.vid) b_nodes.insert(id);
        }
    }
    return b_nodes;
}



// محاسبه مرکز هندسی المان (میانگین مختصات سه گره)
std::array<double, 2> Mesh::get_element_center(int e) const {
    const auto& v = elems[e].vid;
    return {
        (nodes[v[0]].x + nodes[v[1]].x + nodes[v[2]].x) / 3.0,
        (nodes[v[0]].y + nodes[v[1]].y + nodes[v[2]].y) / 3.0
    };
}

// محاسبه فاصله مرکز المان تا نوک ترک
double Mesh::get_dist_to_tip(int e, const std::array<double, 2>& tip_pos) const {
    std::array<double, 2> center = get_element_center(e);
    return std::sqrt(std::pow(center[0] - tip_pos[0], 2) + 
                     std::pow(center[1] - tip_pos[1], 2));
}

double Mesh::average_h() const {
    if (elems.empty()) return 1.0;
    double total_h = 0.0;
    for (const auto& el : elems) {
        // محاسبه طول لبه‌های المان مثلثی
        double d01 = std::sqrt(std::pow(nodes[el.vid[0]].x - nodes[el.vid[1]].x, 2) + 
                               std::pow(nodes[el.vid[0]].y - nodes[el.vid[1]].y, 2));
        double d12 = std::sqrt(std::pow(nodes[el.vid[1]].x - nodes[el.vid[2]].x, 2) + 
                               std::pow(nodes[el.vid[1]].y - nodes[el.vid[2]].y, 2));
        double d20 = std::sqrt(std::pow(nodes[el.vid[2]].x - nodes[el.vid[0]].x, 2) + 
                               std::pow(nodes[el.vid[2]].y - nodes[el.vid[0]].y, 2));
        total_h += (d01 + d12 + d20) / 3.0;
    }
    return total_h / elems.size();
}

double Mesh::local_h_at_tip(int tip_node_id) const {
    double sum_h = 0.0;
    int count = 0;

    for (size_t e = 0; e < elems.size(); ++e) {
        const auto& elem = elems[e];
        bool involves_tip = false;
        
        // Check if this element is connected to the tip node
        for (int i = 0; i < 3; ++i) {
            if (elem.vid[i] == tip_node_id) {
                involves_tip = true;
                break;
            }
        }
        
        if (involves_tip) {
            // Compute average edge length of this element
            for (int i = 0; i < 3; ++i) {
                int v1 = elem.vid[i];
                int v2 = elem.vid[(i+1)%3];
                
                // FIX: Use .x and .y instead of [0] and [1]
                double dx = nodes[v1].x - nodes[v2].x;
                double dy = nodes[v1].y - nodes[v2].y;
                double edge_len = std::sqrt(dx*dx + dy*dy);
                
                sum_h += edge_len;
                count++;
            }
        }
    }
    
    // Return average h, or a default value if no elements found
    return (count > 0) ? (sum_h / count) : 0.01;
}

void Mesh::smooth_mesh_around_tip(int tip_node_id, double dx, double dy) {
    // ۱. پیدا کردن گره‌های لایه اول (گره‌هایی که با یک المان به نوک ترک وصل هستند)
    std::set<int> layer1;
    for (const auto& el : elems) {
        bool touches_tip = (el.vid[0] == tip_node_id || el.vid[1] == tip_node_id || el.vid[2] == tip_node_id);
        if (touches_tip) {
            for (int vid : el.vid) {
                if (vid != tip_node_id) layer1.insert(vid);
            }
        }
    }

    // ۲. پیدا کردن گره‌های لایه دوم (همسایگان لایه اول)
    std::set<int> layer2;
    for (const auto& el : elems) {
        bool touches_layer1 = false;
        for (int vid : el.vid) {
            if (layer1.count(vid)) { touches_layer1 = true; break; }
        }
        if (touches_layer1) {
            for (int vid : el.vid) {
                if (vid != tip_node_id && !layer1.count(vid)) layer2.insert(vid);
            }
        }
    }

    // ۳. جابجایی نرم (Smoothing) با در نظر گرفتن قید مرزی
    // قانون: اگر گره روی لبه‌های بیرونی تیر است، جابجا نشود.
    auto is_boundary = [&](int id) {
        const auto& n = nodes[id];
        return (std::abs(n.y - _ymin) < 1e-7 || std::abs(n.y - _ymax) < 1e-7 ||
                std::abs(n.x - _xmin) < 1e-7 || std::abs(n.x - _xmax) < 1e-7);
    };

    // جابجایی لایه ۱ (۵۰٪ جابجایی نوک ترک)
    for (int id : layer1) {
        if (!is_boundary(id)) {
            nodes[id].x += 0.5 * dx;
            nodes[id].y += 0.5 * dy;
        }
    }

    // جابجایی لایه ۲ (۲۵٪ جابجایی نوک ترک)
    for (int id : layer2) {
        if (!is_boundary(id)) {
            nodes[id].x += 0.25 * dx;
            nodes[id].y += 0.25 * dy;
        }
    }
    
    std::cout << "  [SMOOTHING] " << layer1.size() << " nodes in Layer 1 and " 
              << layer2.size() << " nodes in Layer 2 adjusted.\n";
}
