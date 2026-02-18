#include "vtk.hpp"
#include <fstream>
#include <iomanip>

namespace VTKWriter {

bool write_vtu(const std::string& path,
               const Mesh& mesh,
               const std::vector<std::pair<std::string, std::vector<std::array<double,3>>>>& point_data_vecs,
               const std::vector<std::pair<std::string, std::vector<double>>>& cell_data_scalars) {
    std::ofstream out(path);
    if (!out) return false;

    out << R"(<?xml version="1.0"?>)"
        << "\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n";
    out << "<Piece NumberOfPoints=\"" << mesh.num_nodes() << "\" NumberOfCells=\"" << mesh.num_elements() << "\">\n";

    // Points
    out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(16);
    for (size_t i=0; i<mesh.num_nodes(); ++i) {
        out << mesh.node(i).x << " " << mesh.node(i).y << " 0\n";
    }
    out << "</DataArray>\n</Points>\n";

    // Cells
    out << "<Cells>\n";
    // connectivity
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t e=0; e<mesh.num_elements(); ++e) {
        const auto& el = mesh.element(e);
        out << el.vid[0] << " " << el.vid[1] << " " << el.vid[2] << "\n";
    }
    out << "</DataArray>\n";
    // offsets
    out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int off = 0;
    for (size_t e=0; e<mesh.num_elements(); ++e) { off += 3; out << off << "\n"; }
    out << "</DataArray>\n";
    // types
    out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t e=0; e<mesh.num_elements(); ++e) out << 5 << "\n"; // VTK_TRIANGLE=5
    out << "</DataArray>\n";
    out << "</Cells>\n";

    // PointData
    if (!point_data_vecs.empty()) {
        out << "<PointData>\n";
        for (const auto& pd : point_data_vecs) {
            out << "<DataArray type=\"Float64\" Name=\"" << pd.first << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const auto& vec : pd.second) {
                out << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
            }
            out << "</DataArray>\n";
        }
        out << "</PointData>\n";
    }

    // CellData
    if (!cell_data_scalars.empty()) {
        out << "<CellData>\n";
        for (const auto& cd : cell_data_scalars) {
            out << "<DataArray type=\"Float64\" Name=\"" << cd.first << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (const auto& val : cd.second) out << val << "\n";
            out << "</DataArray>\n";
        }
        out << "</CellData>\n";
    }

    out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    return true;
}

} // namespace VTKWriter

