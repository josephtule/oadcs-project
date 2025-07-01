#pragma once 

#include <fstream>

#include "typedefs.h"
#include "Body.h"

struct CelestialBody : public Body {};

struct EGMCoefficients {
    int max_degree, max_order; // max_degree = N, max_order = M
    std::vector<f64> C;
    std::vector<f64> S;

    // Instantiate and fill with zeros
    EGMCoefficients(int N, int M)
        : max_degree(N), max_order(M), C((N + 1) * (M + 1), 0.),
          S((N + 1) * (M + 1), 0.) {};

    void load_egm(std::string filename, int year) {
        int n_ind_s, m_ind_s, C_ind_s, S_ind_s;
        int n_ind_e, m_ind_e, C_ind_e, S_ind_e;

        std::ifstream file(filename);
        std::string line;
        if (!file) {
            throw std::runtime_error(
                "Cannot open Gravity Model file: " + filename
            );
        }

        switch (year) {
        case 1984: {
            n_ind_s = 0;
            n_ind_e = 4;
            m_ind_s = 5;
            m_ind_e = 9;
            C_ind_s = 10;
            C_ind_e = 24;
            S_ind_s = 25;
            S_ind_e = 39;
            break;
        }
        case 1996: {
            n_ind_s = 0;
            n_ind_e = 3;
            m_ind_s = 4;
            m_ind_e = 7;
            C_ind_s = 8;
            C_ind_e = 27;
            S_ind_s = 28;
            S_ind_e = 47;
            break;
        }
        case 2008: {
            n_ind_s = 0;
            n_ind_e = 4;
            m_ind_s = 5;
            m_ind_e = 9;
            C_ind_s = 13;
            C_ind_e = 34;
            S_ind_s = 38;
            S_ind_e = 59;
            break;
        }
        default:
        }

        int n = 0, m = 0;
        f64 Cnm, Snm;

        while (std::getline(file, line)) {
            if (line.empty()) {
                continue;
            }
            std::string temp;

            temp = line.substr(n_ind_s, n_ind_e - n_ind_s + 1);
            std::erase(temp, ' ');
            n = std::stoi(temp);
            if (n > max_degree) {
                break;
            }

            temp = line.substr(m_ind_s, m_ind_e - m_ind_s + 1);
            std::erase(temp, ' ');
            m = std::stoi(temp);
            if (m > max_order) {
                continue;
            }

            temp = line.substr(C_ind_s, C_ind_e - C_ind_s + 1);
            std::erase(temp, ' ');
            Cnm = std::stod(temp);
            temp = line.substr(S_ind_s, S_ind_e - S_ind_s + 1);
            std::erase(temp, ' ');
            Snm = std::stoi(temp);

            if (n >= 0 && n <= max_degree && m >= 0 && m <= max_order) {
                C[idx(n, m)] = Cnm;
                S[idx(n, m)] = Snm;
            }
        }
    };

    inline int idx(int n, int m) const { return n * (max_order + 1) + m; }
    double getC(int n, int m) const { return C[idx(n, m)]; }
    double getS(int n, int m) const { return S[idx(n, m)]; }
    void setC(int n, int m, double v) { C[idx(n, m)] = v; }
    void setS(int n, int m, double v) { S[idx(n, m)] = v; }
};