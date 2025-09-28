#include "../include/interpolation.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>

extern "C" {
    void dgelss_(int* m, int* n, int* nrhs,
                double* A, int* lda, double* B, int* ldb,
                double* S, double* rcond, int* rank,
                double* work, int* lwork, int* info);
}

std::vector<PointData> read_csv(const std::string& filename) {
    std::vector<PointData> data;
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open CSV file: " + filename);

    std::string line;
    bool header = true;
    while (std::getline(file, line)) {
        if (header) { header = false; continue; }
        std::stringstream ss(line);
        std::string timeStr, pStr, sStr;
        std::getline(ss, timeStr, ',');
        std::getline(ss, pStr, ',');
        std::getline(ss, sStr, ',');
        if (!pStr.empty() && !sStr.empty()) {
            data.push_back({std::stod(sStr), std::stod(pStr)});
        }
    }
    return data;
}

void read_csv_vector(const std::string& filename, 
                    std::vector<PointData>& data0, 
                    std::vector<PointData>& data1, 
                    std::vector<PointData>& data2) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open CSV file: " + filename);

    std::string line;
    bool header = true;
    while (std::getline(file, line)) {
        if (header) { header = false; continue; }
        std::stringstream ss(line);
        std::string timeStr, u0Str, u1Str, u2Str, sStr;
        
        std::getline(ss, timeStr, ',');
        std::getline(ss, u0Str, ',');
        std::getline(ss, u1Str, ',');
        std::getline(ss, u2Str, ',');
        std::getline(ss, sStr, ',');

        if (!u0Str.empty() && !u1Str.empty() && !u2Str.empty() && !sStr.empty()) {
            double s_val = std::stod(sStr);
            data0.push_back({s_val, std::stod(u0Str)});
            data1.push_back({s_val, std::stod(u1Str)});
            data2.push_back({s_val, std::stod(u2Str)});
        }
    }
}

std::function<double(double)> poly_fit(const std::vector<PointData>& data, int degree) {
    int m = data.size();
    int n = degree + 1;

    // Check if data is nearly constant
    double min_val = data[0].val;
    double max_val = data[0].val;
    for (const auto& point : data) {
        if (point.val < min_val) min_val = point.val;
        if (point.val > max_val) max_val = point.val;
    }
    
    if (std::abs(max_val - min_val) < 1e-10) {
        double avg_val = (min_val + max_val) / 2.0;
        return [avg_val](double) { return avg_val; }; // Remove unused parameter 's'
    }

    // Design matrix
    std::vector<double> A(m * n);
    for (int i = 0; i < m; i++) {
        double s = data[i].s;
        A[i] = 1.0; // Column-major order
        for (int j = 1; j < n; j++) {
            A[i + j * m] = std::pow(s, j); // Correct index for column-major
        }
    }

    // RHS
    std::vector<double> B(m);
    for (int i = 0; i < m; i++) {
        B[i] = data[i].val;
    }

    // Use SVD-based solver (dgelss)
    int nrhs = 1;
    int lda = m;
    int ldb = std::max(m, n);
    B.resize(ldb);
    
    std::vector<double> S(std::min(m, n));
    double rcond = -1;
    int rank;
    int lwork = -1;
    int info;
    double wkopt;

    // First call to determine optimal workspace
    dgelss_(&m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, S.data(), &rcond, &rank, &wkopt, &lwork, &info);
    lwork = static_cast<int>(wkopt);
    std::vector<double> work(lwork);
    
    // Second call to solve the system
    dgelss_(&m, &n, &nrhs, A.data(), &lda, B.data(), &ldb, S.data(), &rcond, &rank, work.data(), &lwork, &info);
    
    if (info != 0) {
        double avg_val = 0.0;
        for (const auto& point : data) avg_val += point.val;
        avg_val /= m;
        return [avg_val](double) { return avg_val; }; // Remove unused parameter 's'
    }

    std::vector<double> coeffs(B.begin(), B.begin() + n);

    std::cout << "Fitted polynomial coefficients: ";
    for (double c : coeffs) std::cout << c << " ";
    std::cout << std::endl;

    for (double coeff : coeffs) {
        if (std::isnan(coeff)) {
            double avg_val = 0.0;
            for (const auto& point : data) avg_val += point.val;
            avg_val /= m;
            return [avg_val](double) { return avg_val; }; // Remove unused parameter 's'
        }
    }

    return [coeffs, n](double s) {
        double val = 0.0;
        double power = 1.0;
        for (int i = 0; i < n; i++) {
            val += coeffs[i] * power;
            power *= s;
        }
        return val;
    };
}