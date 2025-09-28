#pragma once

#include <vector>
#include <string>
#include <functional>

struct PointData {
    double s;   // arc_length (independent variable)
    double val; // function value
};

// Reads CSV file with columns: time, value, arc_length
std::vector<PointData> read_csv(const std::string& filename);

void read_csv_vector(const std::string& filename, 
                    std::vector<PointData>& data0, 
                    std::vector<PointData>& data1, 
                    std::vector<PointData>& data2);

// Fits polynomial regression of given degree, returns evaluator lambda
std::function<double(double)> poly_fit(const std::vector<PointData>& data, int degree);
