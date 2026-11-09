/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)
//
// References:
//   Der Kiureghian, A. (1981). A response spectrum method for random
//     vibration analysis of MDF systems. Earthquake Engineering and
//     Structural Dynamics, 9(5), 419-435. doi:10.1002/eqe.4290090503
//   Wilson, E.L., Der Kiureghian, A. and Bayo, E.P. (1981). A replacement
//     for the SRSS method in seismic analysis. Earthquake Engineering and
//     Structural Dynamics, 9(2), 187-192. doi:10.1002/eqe.4290090207
//   U.S. Nuclear Regulatory Commission (1976). Regulatory Guide 1.92,
//     Rev. 1: Combining modal responses and spatial components in seismic
//     response analysis.
//   Chopra, A.K. (2017). Dynamics of Structures, 5th ed. Pearson.
//
// Description: modal combination rules for response-spectrum analysis,
// header only.
//
//   Given the per-mode peak values R_a of ONE response quantity, combine the
//   non-simultaneous modal peaks into a single design value:
//     ABS         R = sum |R_a|                         (upper bound)
//     SRSS        R = sqrt(sum R_a^2)                    (well-separated modes)
//     TenPercent  R = sqrt(sum R_a^2 + 2 sum_{|Ti-Tj|<=0.1 min} |Ri Rj|)
//                                                        (RG 1.92, Rev. 1)
//     CQC         R = sqrt(sum_i sum_j rho_ij R_i R_j)   (closely-spaced modes)
//   with the correlation coefficient of Der Kiureghian (1981)
//     rho_ij = 8 sqrt(xi_i xi_j)(xi_i + r xi_j) r^{3/2}
//              / [ (1-r^2)^2 + 4 xi_i xi_j r (1+r^2) + 4 (xi_i^2+xi_j^2) r^2 ],
//     r = w_j / w_i,  rho_ii = 1.
//
//   Combination is a per-quantity, nonlinear operation: combining nodal
//   displacements and then deriving element forces from them is not the
//   combination of the forces. ResponseSpectrumAnalysis -combine combines
//   the nodal displacement field; other quantities must combine their own
//   per-mode peaks.
#ifndef LadrunoModalCombination_h
#define LadrunoModalCombination_h

#include <vector>
#include <cmath>

namespace LadrunoModalCombination {

enum Rule { RULE_NONE = -1, RULE_SRSS = 0, RULE_CQC = 1, RULE_ABS = 2, RULE_TENPERCENT = 3 };

// parse a rule name (case-insensitive-ish, accepts the common spellings);
// returns RULE_NONE if unrecognized.
inline Rule parseRule(const char* s)
{
    if (s == 0) return RULE_NONE;
    // compare against the accepted spellings
    auto ieq = [](const char* a, const char* b) {
        for (; *a && *b; ++a, ++b) {
            char ca = *a, cb = *b;
            if (ca >= 'A' && ca <= 'Z') ca = char(ca - 'A' + 'a');
            if (cb >= 'A' && cb <= 'Z') cb = char(cb - 'A' + 'a');
            if (ca != cb) return false;
        }
        return *a == 0 && *b == 0;
    };
    if (ieq(s, "SRSS")) return RULE_SRSS;
    if (ieq(s, "CQC")) return RULE_CQC;
    if (ieq(s, "ABS")) return RULE_ABS;
    if (ieq(s, "TenPercent") || ieq(s, "TenPct") || ieq(s, "10pct")) return RULE_TENPERCENT;
    return RULE_NONE;
}

// Der Kiureghian correlation coefficient (0 < rho <= 1, symmetric, rho_ii = 1).
// A non-oscillatory mode (w <= 0: rigid / spurious) is uncorrelated with any
// other by convention: return 0 rather than divide by zero (r = wj/wi = inf).
inline double rhoCQC(double wi, double wj, double xi_i, double xi_j)
{
    if (wi <= 0.0 || wj <= 0.0)
        return 0.0;
    const double r = wj / wi;
    const double num = 8.0 * std::sqrt(xi_i * xi_j) * (xi_i + r * xi_j)
                       * std::pow(r, 1.5);
    const double den = (1.0 - r * r) * (1.0 - r * r)
                       + 4.0 * xi_i * xi_j * r * (1.0 + r * r)
                       + 4.0 * (xi_i * xi_i + xi_j * xi_j) * r * r;
    // den -> 0 only when r == 1 AND both modes are undamped (num -> 0 too): the
    // 0/0 limit of two coincident modes is perfect correlation -> rho = 1.
    if (den <= 0.0)
        return 1.0;
    return num / den;
}

// Combine the per-mode peak values R of ONE response quantity.
//   R, w (rad/s), xi, T (period, s) are all indexed by mode (same length).
//   T is used only by TenPercent; xi/w only by CQC. Returns a non-negative value.
inline double combine(Rule rule,
                      const std::vector<double>& R,
                      const std::vector<double>& w,
                      const std::vector<double>& xi,
                      const std::vector<double>& T)
{
    const int n = static_cast<int>(R.size());
    if (rule == RULE_ABS) {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += std::fabs(R[i]);
        return s;
    }
    if (rule == RULE_SRSS) {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += R[i] * R[i];
        return std::sqrt(s);
    }
    if (rule == RULE_TENPERCENT) {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += R[i] * R[i];
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                const double tmin = (T[i] < T[j]) ? T[i] : T[j];
                if (std::fabs(T[i] - T[j]) <= 0.10 * tmin)
                    s += 2.0 * std::fabs(R[i] * R[j]);
            }
        return std::sqrt(s > 0.0 ? s : 0.0);
    }
    // CQC
    double s = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double rho = (i == j) ? 1.0 : rhoCQC(w[i], w[j], xi[i], xi[j]);
            s += rho * R[i] * R[j];
        }
    return std::sqrt(s > 0.0 ? s : 0.0);
}

} // namespace LadrunoModalCombination

#endif
