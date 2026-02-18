#include "field_transfer.hpp"
#include <iostream>
#include <limits>

// تابعی کمکی برای محاسبه مختصات باریسنتریک و بررسی اینکه آیا نقطه داخل مثلث است یا خیر
bool is_point_in_tri(double px, double py, 
                     double x1, double y1, double x2, double y2, double x3, double y3,
                     double& w1, double& w2, double& w3) {
    double det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
    w1 = ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) / det;
    w2 = ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) / det;
    w3 = 1.0 - w1 - w2;
    
    double tol = 1e-6;
    return (w1 >= -tol && w2 >= -tol && w3 >= -tol);
}

bool FieldTransfer::transfer_displacement(const Mesh& old_mesh,
                                          const Mesh& new_mesh,
                                          const std::vector<double>& u_old,
                                          std::vector<double>& u_new) const {
    
    size_t num_nodes_new = new_mesh.num_nodes();
    u_new.assign(num_nodes_new * 2, 0.0);

    for (size_t i = 0; i < num_nodes_new; ++i) {
        auto node_new = new_mesh.node(i);
        bool found = false;

        // جستجو در المان‌های مش قدیمی (برای بهینه‌سازی می‌توان از Spatial Index استفاده کرد)
        for (size_t e = 0; e < old_mesh.num_elements(); ++e) {
            const auto& el = old_mesh.element(e);
            auto n1 = old_mesh.node(el.vid[0]);
            auto n2 = old_mesh.node(el.vid[1]);
            auto n3 = old_mesh.node(el.vid[2]);

            double w1, w2, w3;
            if (is_point_in_tri(node_new.x, node_new.y, n1.x, n1.y, n2.x, n2.y, n3.x, n3.y, w1, w2, w3)) {
                // درون‌یابی جابجایی در دو راستای x و y
                u_new[2 * i]     = w1 * u_old[2 * el.vid[0]]     + w2 * u_old[2 * el.vid[1]]     + w3 * u_old[2 * el.vid[2]];
                u_new[2 * i + 1] = w1 * u_old[2 * el.vid[0] + 1] + w2 * u_old[2 * el.vid[1] + 1] + w3 * u_old[2 * el.vid[2] + 1];
                found = true;
                break;
            }
        }

        // اگر گره بیرون از مش قبلی بود (به دلیل تغییرات هندسی کوچک)، نزدیک‌ترین گره را پیدا کن
        if (!found) {
            int closest = old_mesh.find_closest_node(node_new.x, node_new.y);
            u_new[2 * i]     = u_old[2 * closest];
            u_new[2 * i + 1] = u_old[2 * closest + 1];
        }
    }

    return true;
}
