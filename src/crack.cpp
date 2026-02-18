#include "crack.hpp"
#include <cmath>
#include <iostream>

CrackPath::CrackPath() : has_tip_(false) {
    tip_[0] = 0.0;
    tip_[1] = 0.0;
}

void CrackPath::reset() {
    has_tip_ = false;
    segments_.clear();
    tip_[0] = 0.0;
    tip_[1] = 0.0;
}

void CrackPath::set_initial_tip(const std::array<double,2>& tip) {
    tip_ = tip;
    has_tip_ = true;
}

void CrackPath::add_segment(const CrackSegment& seg) {
    segments_.push_back(seg);
    tip_ = seg.end;
    has_tip_ = true;
}

std::array<double,2> CrackPath::tip() const {
    if (!has_tip_) {
        std::cerr << "[CrackPath] Warning: Tip requested but not set.\n";
        std::array<double,2> zero;
        zero[0] = 0.0;
        zero[1] = 0.0;
        return zero;
    }
    return tip_;
}

const std::vector<CrackSegment>& CrackPath::segments() const {
    return segments_;
}

double CrackPath::total_length() const {
    double len = 0.0;
    for (size_t i = 0; i < segments_.size(); ++i) {
        const CrackSegment& seg = segments_[i];
        const double dx = seg.end[0] - seg.start[0];
        const double dy = seg.end[1] - seg.start[1];
        len += std::sqrt(dx*dx + dy*dy);
    }
    return len;
}

bool CrackPath::advance_direction(const std::array<double,2>& dir, double step) {
    if (!has_tip_) {
        std::cerr << "[CrackPath] Cannot advance: tip not set\n";
        return false;
    }
    
    // Normalize direction
    double norm = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
    if (norm < 1e-12) {
        std::cerr << "[CrackPath] Invalid direction vector (zero length)\n";
        return false;
    }
    
    std::array<double,2> unit_dir;
    unit_dir[0] = dir[0] / norm;
    unit_dir[1] = dir[1] / norm;
    
    // Create new segment
    CrackSegment seg;
    seg.start = tip_;
    seg.end[0] = tip_[0] + step * unit_dir[0];
    seg.end[1] = tip_[1] + step * unit_dir[1];
    
    add_segment(seg);
    return true;
}

void CrackPath::debug_print() const {
    std::cout << "[CrackPath Debug]\n";
    std::cout << "  Has tip: " << (has_tip_ ? "Yes" : "No") << "\n";
    if (has_tip_) {
        std::cout << "  Tip position: (" << tip_[0] << ", " << tip_[1] << ")\n";
    }
    std::cout << "  Number of segments: " << segments_.size() << "\n";
    std::cout << "  Total length: " << total_length() << " m\n";
    
    for (size_t i = 0; i < segments_.size(); ++i) {
        const CrackSegment& seg = segments_[i];
        std::cout << "  Segment " << i << ": (" 
                  << seg.start[0] << ", " << seg.start[1] << ") -> ("
                  << seg.end[0] << ", " << seg.end[1] << ")\n";
    }
}
