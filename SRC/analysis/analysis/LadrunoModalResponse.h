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
//   Nigam, N.C. and Jennings, P.C. (1969). Calculation of response spectra
//     from strong-motion earthquake records. Bulletin of the Seismological
//     Society of America, 59(2), 909-922. doi:10.1785/BSSA0590020909
//   Van Loan, C.F. (1978). Computing integrals involving the matrix
//     exponential. IEEE Transactions on Automatic Control, 23(3), 395-404.
//     doi:10.1109/TAC.1978.1101743
//   Chopra, A.K. (2017). Dynamics of Structures, 5th ed. Pearson
//     (Sec. 5.2, piecewise-exact recurrence).
//   Clough, R.W. and Penzien, J. (1993). Dynamics of Structures, 2nd ed.
//     McGraw-Hill.
//   Rice, S.O. (1944, 1945). Mathematical analysis of random noise. Bell
//     System Technical Journal, 23(3), 282-332 and 24(1), 46-156.
//     doi:10.1002/j.1538-7305.1944.tb00874.x,
//     doi:10.1002/j.1538-7305.1945.tb00453.x
//   Davenport, A.G. (1964). Note on the distribution of the largest value of
//     a random function with application to gust loading. Proceedings of
//     the Institution of Civil Engineers, 28(2), 187-196.
//     doi:10.1680/iicep.1964.10112
//
// Description: modal-superposition response of a linear model on the modes
// of the last `eigen` + `modalProperties` run.
//
//   modalResponseHistory: the piecewise-linear modal-superposition
//   transient.  Each retained mode is a decoupled SDOF
//
//       q_a'' + d_a q_a' + w_a^2 q_a = f_a(t),   d_a = c_a/m_a,
//
//   integrated by a per-mode constant recurrence
//
//       [q_{n+1}; q'_{n+1}] = A_a [q_n; q'_n] + B_a [f_n; f_{n+1}]
//
//   that is exact for a load varying linearly within each step (Nigam and
//   Jennings 1969; first-order hold).  The physical response is recovered by
//   the modal sum u = sum_a (phi_a * Vscale) q_a and one domain step is
//   committed per time station, so ordinary recorders capture displacement,
//   velocity, acceleration, element and reaction histories as in a direct
//   run.
//
//   Excitation channels:
//     -baseAccel: uniform u_g''(t) along a global direction, relative
//         formulation, modal load f_a = -Gamma_a u_g''(t) with Gamma_a the
//         participation factor computed by modalProperties.
//     -load: nodal forces P(t) = s(t) * P.  P is a load pattern's
//         plain-NodalLoad shape (the pattern's own TimeSeries is ignored) and
//         s(t) the -series time series sampled at the stations; modal load
//         f_a(t) = (psi_a^T P / m~_a) s(t) with psi = evec*Vscale and
//         m~ = generalizedMasses(), independent of the eigenvector
//         normalization.  The response is absolute for -load (the relative
//         formulation is a base-motion concept).
//
//   The (w, d) parametrization (rather than (w, xi)) covers all damping
//   branches through the discriminant Delta = w^2 - (d/2)^2 (under,
//   critical, over) and a rigid-body mode carrying Rayleigh mass-proportional
//   damping (w=0, d=a0 > 0), which a pure q''=f rigid branch would miss.
//
//   frequencyResponse / steadyStateDynamics: the steady harmonic response
//   to a uniform base acceleration u_g''(t) = amp e^{+iOm t} (relative
//   formulation) is the modal sum
//
//       u_hat(Om) = amp * sum_a (phi_a Vscale_a) (-Gamma_a) H_a(Om),
//       H_a(Om)   = 1 / (w_a^2 - Om^2 + i Om d_a),   d_a = c_a/m_a,
//
//   the frequency-domain image of the transient recovery above (same
//   psi_a = phi_a Vscale_a and Gamma_a).  Velocity = i Om u_hat; relative
//   acceleration = -Om^2 u_hat.  The time convention is e^{+iOm t} (the
//   +i Om d_a term gives the -90 deg response lag at resonance); it agrees
//   with the direct complex solve (K - Om^2 M + i Om C)^{-1}(-M R).
//   Frequencies are in Hz (Om = 2 pi f), in and out.  With -load the
//   excitation is amp P e^{+iOm t} and f_a = psi_a^T P / m~_a; -baseAccel is
//   the special case P = -M R (then f_a = -Gamma_a).
//
//   randomResponse: for a stationary excitation given as a one-sided power
//   spectral density G(f) in Hz (sigma^2 = int_0^inf G df, the convention of
//   wind, floor-vibration and equipment specifications), the response
//   auto-PSD is G_xx(f) = |H_x(f)|^2 G(f) with H_x the frequencyResponse
//   transfer function, and the band moments over the sweep grid are
//   m_k = int f^k G_xx df (k = 0, 2, trapezoid rule): RMS = sqrt(m_0), mean
//   zero-upcrossing rate nu_0 = sqrt(m_2/m_0) [Hz] (Rice 1944), and expected
//   peak over a duration T, g * RMS with the Davenport (1964) peak factor
//   g = sqrt(2 ln(nu_0 T)) + 0.5772/sqrt(2 ln(nu_0 T)), for nu_0 T > 1.
//   Against the two-sided rad/s PSD S(Om) of random-vibration texts,
//   G(f) = 4 pi S(Om); the white-noise SDOF variance is
//   sigma_x^2 = G0/(8 xi w^3).

#ifndef LadrunoModalResponse_h
#define LadrunoModalResponse_h

#include <vector>
#include <string>

class AnalysisModel;
class Domain;
class TimeSeries;
class DomainModalProperties;

// ---------------------------------------------------------------------------
// Shared modal damping specification: how the per-mode d_a = c_a/m_a (= 2 xi w
// for w>0) is obtained.  Used by the frequency-domain commands; the transient
// carries its own equivalent copy.
// ---------------------------------------------------------------------------
struct LadrunoModalDamping
{
    enum Kind { UNIFORM, RAYLEIGH, MODAL_LIST };
    Kind                kind = UNIFORM;
    double              xi   = 0.0;   // uniform ratio
    double              a0   = 0.0;   // Rayleigh mass factor
    double              a1   = 0.0;   // Rayleigh stiffness factor
    std::vector<double> xis;          // per absolute-mode ratios (MODAL_LIST)

    // d_a = c_a/m_a for absolute 0-based mode index at circular frequency w.
    // Rayleigh: a0 + a1 w^2 (rigid mode w=0 -> a0). Ratio channels: 2 xi w
    // (rigid mode -> 0, stiffness-proportional damping vanishing there).
    double coeff(int mode0, double w) const;
};

class LadrunoModalResponse
{
public:
    // damping specification (how the per-mode d_a = c_a/m_a is obtained)
    enum DampingKind { DAMP_UNIFORM, DAMP_RAYLEIGH, DAMP_MODAL_LIST };

    LadrunoModalResponse(AnalysisModel* theModel,
                         TimeSeries*    baseAccel,
                         int            direction,
                         double         dt,
                         int            nsteps,
                         double         t0);
    ~LadrunoModalResponse();

    // damping channels (mutually exclusive; last one set wins)
    void setDampingUniform(double xi);
    void setDampingRayleigh(double a0, double a1);
    void setDampingModalList(const std::vector<double>& xi_per_mode); // absolute mode order

    // optional 1-based mode subset (default: all modes from the last eigen run)
    void setModeSubset(const std::vector<int>& modes_1based);

    // switch the excitation to a nodal-force pattern P(t) = s(t)*P, where P is
    // the pattern's plain-NodalLoad shape and s(t) is the ctor's time series
    // (the `baseAccel` ctor argument doubles as s; `direction` is unused)
    void setLoadPattern(int patternTag)
    { m_useLoad = true; m_patternTag = patternTag; }

    // run the whole transient, committing one domain step per time station
    int analyze();

    // exposed for the unit oracle: exact PWL recurrence blocks for one SDOF.
    //   A[2][2], B[2][2] s.t. [q1;qd1] = A[q0;qd0] + B[f0;f1].
    // w  = undamped circular frequency (>= 0), d = c/m damping coefficient (>= 0).
    static void recurrenceBlocks(double w, double d, double dt,
                                 double A[2][2], double B[2][2]);

private:
    int  check(DomainModalProperties& mp);
    // per-mode damping coefficient d_a = c_a/m_a for absolute mode index (0-based)
    double dampingCoeff(int mode0, double w) const;

private:
    AnalysisModel* m_model;
    TimeSeries*    m_baseAccel;   // excitation series: u_g''(t), or s(t) w/ -load
    int            m_direction;   // 1-based global excitation direction (base)
    bool           m_useLoad = false; // excitation = nodal-force pattern
    int            m_patternTag = 0;  // load-pattern tag (m_useLoad)
    double         m_dt;
    int            m_nsteps;
    double         m_t0;          // start time (series sampled at t0 + n*dt)

    DampingKind         m_damp_kind;
    double              m_xi;      // uniform ratio
    double              m_a0, m_a1; // Rayleigh factors
    std::vector<double> m_xi_list; // per absolute-mode ratios

    std::vector<int>    m_modes;   // absolute 0-based mode indices to include
};


// ===========================================================================
// frequencyResponse / steadyStateDynamics: modal frequency response and
// steady-state amplitude (see the description at the top of this file).
//
//   EX_BASE_ACCEL  uniform base acceleration, modal load -Gamma_a
//   EX_LOAD        nodal-force pattern amp P e^{+iOm t}; P is assembled from
//                  the pattern's plain NodalLoads (reference values; the
//                  pattern's own TimeSeries is ignored, the sweep replaces
//                  it) and the modal load is
//
//                      f_a = psi_a^T P / m~_a,   m~_a = psi_a^T M psi_a,
//
//                  with psi_a = evec_a*Vscale_a, the basis in which
//                  DomainModalProperties computes generalizedMasses() (the
//                  -unorm scaling is applied before the generalized-mass and
//                  participation products), so m~_a = generalizedMasses()(a).
// ===========================================================================
class LadrunoFrequencyResponse
{
public:
    enum SweepKind { SWEEP_LIN, SWEEP_LOG, SWEEP_BIASED };
    enum RespKind  { RESP_DISP, RESP_VEL, RESP_ACCEL };
    enum ExKind    { EX_BASE_ACCEL, EX_LOAD };

    explicit LadrunoFrequencyResponse(AnalysisModel* theModel);
    ~LadrunoFrequencyResponse();

    // excitation: uniform base acceleration of amplitude amp along a global dir
    void setBaseAccel(int direction, double amp)
    { m_ex = EX_BASE_ACCEL; m_direction = direction; m_amp = amp; }
    // excitation: harmonic nodal-force pattern, amp * P * e^{+iOm t}
    void setLoadPattern(int patternTag, double amp)
    { m_ex = EX_LOAD; m_patternTag = patternTag; m_amp = amp; }
    // frequency sweep in Hz
    void setSweep(double fmin, double fmax, int nf, SweepKind k)
    { m_fmin = fmin; m_fmax = fmax; m_nf = nf; m_sweep = k; }
    void setDamping(const LadrunoModalDamping& d) { m_damp = d; }
    void setModeSubset(const std::vector<int>& modes_1based);
    void setResponse(int nodeTag, int dof, RespKind r)
    { m_nodeTag = nodeTag; m_dof = dof; m_resp = r; }
    void setOutputFile(const char* fname) { m_outfile = fname ? fname : ""; }

    // Evaluate the sweep.  Fills `rows`, one per frequency:
    //   ampOnly==false -> {f, Re, Im}   (full complex FRF)
    //   ampOnly==true  -> {f, |H|}      (steady-state amplitude, SSD)
    // Also writes `rows` to the output file if one was set. Returns 0 / -1.
    int run(std::vector<std::vector<double> >& rows, bool ampOnly);

private:
    int  buildGrid(std::vector<double>& freqs) const; // Hz sample points

    AnalysisModel*      m_model;
    ExKind              m_ex;        // excitation channel
    int                 m_direction; // 1-based global excitation dir (EX_BASE_ACCEL)
    int                 m_patternTag;// load pattern tag (EX_LOAD)
    double              m_amp;
    double              m_fmin, m_fmax;
    int                 m_nf;
    SweepKind           m_sweep;
    LadrunoModalDamping m_damp;
    std::vector<int>    m_modes;     // absolute 0-based mode subset (empty=all)
    int                 m_nodeTag;
    int                 m_dof;       // 1-based response DOF
    RespKind            m_resp;
    std::string         m_outfile;
};


// ===========================================================================
// randomResponse: stationary random response, PSD -> RMS (see the
// description at the top of this file).
//
//   A zero-damped mode whose resonance lies inside the sweep band is refused:
//   its variance integral diverges, so any finite number the grid returned
//   would be meaningless (and a -biased grid samples exactly on the
//   resonance).  Zero-damped modes strictly outside the band keep a finite
//   in-band integrand and are accepted.
//
//   -load: the excitation is P(t) = s(t) * P with P the pattern's nodal-load
//   shape and s(t) a stationary scalar of one-sided PSD G_s(f) [Hz] (all
//   loads fully correlated); H_x is then the -load transfer function and
//   G_xx = |H_x|^2 G_s.  For a white force PSD G_F on an SDOF,
//   sigma_u^2 = G_F/(8 xi w^3 m^2).
// ===========================================================================
class LadrunoRandomResponse
{
public:
    struct Stats {
        double rms  = 0.0;   // sqrt(m0)
        double m0   = 0.0;   // int   G_xx df   over the band
        double m2   = 0.0;   // int f^2 G_xx df over the band
        double nu0  = 0.0;   // sqrt(m2/m0), Hz
        double peak = 0.0;   // Davenport expected peak (valid iff hasPeak)
        bool   hasPeak = false;
    };

    explicit LadrunoRandomResponse(AnalysisModel* theModel);
    ~LadrunoRandomResponse();

    // excitation: uniform base acceleration along a global dir, with one-sided
    // PSD G(f) [Hz] sampled from the time series (getFactor(f), f in Hz)
    void setBaseAccel(int direction)
    { m_ex = LadrunoFrequencyResponse::EX_BASE_ACCEL; m_direction = direction; }
    // excitation: nodal-force pattern shape P scaled by a stationary scalar
    // s(t) with one-sided PSD G_s(f) [Hz]: P(t) = s(t)*P (fully correlated)
    void setLoadPattern(int patternTag)
    { m_ex = LadrunoFrequencyResponse::EX_LOAD; m_patternTag = patternTag; }
    void setInputPSD(TimeSeries* psd)    { m_psd = psd; }
    void setSweep(double fmin, double fmax, int nf,
                  LadrunoFrequencyResponse::SweepKind k)
    { m_fmin = fmin; m_fmax = fmax; m_nf = nf; m_sweep = k; }
    void setDamping(const LadrunoModalDamping& d) { m_damp = d; }
    void setModeSubset(const std::vector<int>& modes_1based)
    { m_modes_1based = modes_1based; }
    void setResponse(int nodeTag, int dof, LadrunoFrequencyResponse::RespKind r)
    { m_nodeTag = nodeTag; m_dof = dof; m_resp = r; }
    void setOutputFile(const char* fname) { m_outfile = fname ? fname : ""; }
    void setDuration(double T)            { m_duration = T; } // for E[peak]

    // Evaluate: FRF sweep -> G_xx rows {f, G_in, G_xx} -> band moments/stats.
    // Also writes `rows` to the output file if one was set. Returns 0 / -1.
    int run(std::vector<std::vector<double> >& rows, Stats& out);

private:
    AnalysisModel*      m_model;
    TimeSeries*         m_psd;
    LadrunoFrequencyResponse::ExKind m_ex;
    int                 m_direction;
    int                 m_patternTag;
    double              m_fmin, m_fmax;
    int                 m_nf;
    LadrunoFrequencyResponse::SweepKind m_sweep;
    LadrunoModalDamping m_damp;
    std::vector<int>    m_modes_1based;  // pass-through to the FRF (empty=all)
    int                 m_nodeTag;
    int                 m_dof;
    LadrunoFrequencyResponse::RespKind  m_resp;
    std::string         m_outfile;
    double              m_duration;      // <=0: no peak estimate
};

#endif
