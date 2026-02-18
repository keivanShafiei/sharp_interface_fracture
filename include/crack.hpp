#ifndef CRACK_HPP
#define CRACK_HPP

#include <cmath>
#include <vector>
#include <array>
#include <iostream>

struct CrackSegment {
    std::array<double,2> start;
    std::array<double,2> end;
};

class CrackPath {
public:
    CrackPath();
    
    void reset();
    void set_initial_tip(const std::array<double,2>& tip);
    void add_segment(const CrackSegment& seg);
    
    bool has_tip() const { return has_tip_; }
    std::array<double,2> tip() const;
    
    const std::vector<CrackSegment>& segments() const;
    double total_length() const;
    
    // Additional methods from original crack.cpp
    bool advance_direction(const std::array<double,2>& dir, double step);
    void debug_print() const;

private:
    bool has_tip_;
    std::array<double,2> tip_;
    std::vector<CrackSegment> segments_;
};

#endif // CRACK_HPP
