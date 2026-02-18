#include "remesher.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef K::Point_2 Point;

Remesher::Remesher(RemeshOptions opts) : opts_(opts) {
    // در اینجا می‌توانید مقداردهی‌های اولیه دیگر را انجام دهید
}

bool Remesher::remesh_local(Mesh& mesh, const CrackPath& crack, bool mesh_changed) {
    auto tip = crack.tip();
    
    // ۱. استخراج گره‌های مرزی قبل از حذف المان‌ها
    std::set<int> boundary_nodes = mesh.get_boundary_nodes_of_removed_region(tip, opts_.radius);
    
    // ۲. حذف المان‌های قدیمی
    mesh.remove_elements_in_radius(tip, opts_.radius);

    CDT cdt;

    // ۳. تزریق سگمنت‌های ترک به عنوان قید (Constraint)
    for (const auto& seg : crack.segments()) {
        cdt.insert_constraint(Point(seg.start[0], seg.start[1]), Point(seg.end[0], seg.end[1]));
    }

// ۴. تزریق گره‌های مرزی به CDT برای حفظ پیوستگی
    for (int node_id : boundary_nodes) {
        auto coord = mesh.node_coords(node_id); // این خروجی std::array<double, 2> است
        
        // به جای coord.x و coord.y از ایندکس‌ها استفاده کنید:
        cdt.insert(Point(coord[0], coord[1])); 
    }

    // ۵. تبدیل مثلث‌های تولید شده به المان‌های سیستم
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        int v1 = mesh.get_or_add_node(fit->vertex(0)->point().x(), fit->vertex(0)->point().y());
        int v2 = mesh.get_or_add_node(fit->vertex(1)->point().x(), fit->vertex(1)->point().y());
        int v3 = mesh.get_or_add_node(fit->vertex(2)->point().x(), fit->vertex(2)->point().y());
        
        if (mesh.compute_element_area(v1, v2, v3) > 1e-9) {
            mesh.add_element(v1, v2, v3);
        }
    }

    return true;
}
