#pragma once
#include "mesh.hpp"
#include <string>
#include <vector>
#include <array>

namespace VTKWriter {

// نوشتن VTU (XML) با داده‌های نقطه‌ای و المانی
// point_data_vecs: نام و بردار داده‌های برداری (مثلاً جابجایی)
// cell_data_scalars: نام و مجموعه اسکالرهای المانی (مثلاً Sxx, Syy, Sxy)
bool write_vtu(const std::string& path,
               const Mesh& mesh,
               const std::vector<std::pair<std::string, std::vector<std::array<double,3>>>>& point_data_vecs,
               const std::vector<std::pair<std::string, std::vector<double>>>& cell_data_scalars);

} // namespace VTKWriter

