/*---------------------------------------------------------------------------*\
    udrgm_gerg.c - Fluent 2026 R1 UDRGM wrapper over the teqp GERG-2008
    kernel (gerg_kernel.cpp). Route mirrors RFAM: GERG-2008 primary, PR
    fallback reserved (ideal-gas guard below the GERG failure path).

    30-species mixture (sandiego30, Fluent-safe names):
      N2 H O2 OH O H2 H2O HO2 H2O2 CO CO2 HCO CH3 CH4 CH2O CH2S
      C2H4 CH3O C2H5 C2H6 C2H2 C2H3 CH2CHO CH2CO CH3CHO C3H4 C3H6 C3H8
      IC3H7 NC3H7
\*---------------------------------------------------------------------------*/
#include "udf.h"  /* pulls in real_gas.h (RGAS_Functions) */
#include "GergUdrgmNasa7.H"
#include <math.h>
#include <string.h>
#include <stdio.h>

#define GERG_NSPECIES 30
#define GERG_RU 8.314462618  /* J/mol/K */

/* ---- kernel (gerg_kernel.cpp, C++ TU) ---- */
extern int gerg_init(const char* fluiddata, const char* const* names, int n);
extern int gerg_last_error(void);
extern int gerg_rho_pT(double T, double p, const double* y, int n, double* rho);
extern int gerg_p_rhoT(double T, double rho, const double* y, int n, double* p);
extern int gerg_e_rhoT(double T, double rho, const double* y, int n, double* e);
extern int gerg_h_rhoT(double T, double rho, const double* y, int n, double* h);
extern int gerg_a_rhoT(double T, double rho, const double* y, int n, double* a);
extern int gerg_cp_cv(double T, double rho, const double* y, int n,
                      double* cp, double* cv);
extern int gerg_mu_k(double T, double rho, const double* y, int n,
                     double* mu, double* k);

static int gerg_initialized = 0;
static const char* gerg_fluiddata = "/root/teqp/dev/fluiddata";

/* Fluent mixture species names (must match Fluent materials panel order) */
static const char* gerg_names[GERG_NSPECIES] = {
    "N2", "H", "O2", "OH", "O", "H2", "H2O", "HO2", "H2O2", "CO",
    "CO2", "HCO", "CH3", "CH4", "CH2O", "CH2S", "C2H4", "CH3O", "C2H5",
    "C2H6", "C2H2", "C2H3", "CH2CHO", "CH2CO", "CH3CHO", "C3H4", "C3H6",
    "C3H8", "IC3H7", "NC3H7"
};

static int (*gerg_msg)(const char*, ...) = NULL;
static void (*gerg_err)(const int, const char*, const char*) = NULL;
static int gerg_silent = 0;
static int gerg_fail_count = 0;

static void gerg_ensure_init(void)
{
    if (gerg_initialized) return;
    const char* fd = getenv("RFAM_GERG_FLUIDDATA");
    if (!fd || !*fd) fd = gerg_fluiddata;
    if (gerg_init(fd, gerg_names, GERG_NSPECIES) == 0)
        gerg_initialized = 1;
    else if (gerg_msg)
        gerg_msg("GERG UDRGM: kernel init FAILED (fluiddata=%s), ideal-gas "
                 "fallback only\n", fd);
}

/* ---- mass-fraction hygiene (mirrors PR UDRGM) ---- */
static void gerg_ensure_yi(cell_t c, Thread* t, double yi[])
{
    int i;
    double sum = 0.0;
    int ok = 1;
    if (yi)
    {
        for (i = 0; i < GERG_NSPECIES; ++i)
        {
            double v = yi[i];
            if (!isfinite(v) || v < -1e-12 || v > 1.2) { ok = 0; break; }
            if (v > 0.0) sum += v;
        }
        if (ok && sum > 0.5 && sum < 1.5)
        {
            for (i = 0; i < GERG_NSPECIES; ++i)
                if (yi[i] < 0.0) yi[i] = 0.0;
            return;
        }
    }
    if (c)
    {
        double csum = 0.0;
        for (i = 0; i < GERG_NSPECIES; ++i)
        {
            double v = C_YI(c, t, i);
            if (isfinite(v) && v > 0.0) csum += v;
            yi[i] = (isfinite(v) && v > 0.0) ? v : 0.0;
        }
        if (csum > 0.5 && csum < 1.5) return;
    }
    for (i = 0; i < GERG_NSPECIES; ++i) yi[i] = 0.0;
    yi[2] = 1.0; /* pure O2 */
}

/* ---- NASA7 ideal helpers (from generated header) ---- */
static double gerg_nasa_cp(int idx, double T)
{
    const double* a = (T >= GERG_NASA7[idx].Tmid)
        ? GERG_NASA7[idx].hi : GERG_NASA7[idx].lo;
    return GERG_RU * (a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T)
        / GERG_NASA7[idx].mw;
}
static double gerg_nasa_h(int idx, double T)
{
    const double* a = (T >= GERG_NASA7[idx].Tmid)
        ? GERG_NASA7[idx].hi : GERG_NASA7[idx].lo;
    return GERG_RU * T * (a[0] + a[1]*T/2 + a[2]*T*T/3 + a[3]*T*T*T/4
        + a[4]*T*T*T*T/5 + a[5]/T) / GERG_NASA7[idx].mw;
}
static double gerg_nasa_s(int idx, double T)
{
    const double* a = (T >= GERG_NASA7[idx].Tmid)
        ? GERG_NASA7[idx].hi : GERG_NASA7[idx].lo;
    return GERG_RU * (a[0]*log(T) + a[1]*T + a[2]*T*T/2 + a[3]*T*T*T/3
        + a[4]*T*T*T*T/4 + a[6]) / GERG_NASA7[idx].mw;
}
static double gerg_mw_from_y(const double yi[])
{
    double m = 0.0;
    int i;
    for (i = 0; i < GERG_NSPECIES; ++i) m += yi[i] / GERG_NASA7[i].mw;
    return m > 1e-300 ? 1.0 / m : 0.0;
}
static double gerg_cp_ideal(const double yi[], double T)
{
    double s = 0.0;
    int i;
    for (i = 0; i < GERG_NSPECIES; ++i) s += yi[i] * gerg_nasa_cp(i, T);
    return s;
}
static double gerg_h_ideal(const double yi[], double T)
{
    double s = 0.0;
    int i;
    for (i = 0; i < GERG_NSPECIES; ++i) s += yi[i] * gerg_nasa_h(i, T);
    return s;
}
static double gerg_s_ideal(const double yi[], double T, double p)
{
    double s = 0.0;
    int i;
    for (i = 0; i < GERG_NSPECIES; ++i)
        s += yi[i] * gerg_nasa_s(i, T);
    /* pressure term: -R/MW * ln(p/p0), p0 = 101325 Pa */
    s -= GERG_RU / gerg_mw_from_y(yi) * log(p / 101325.0);
    return s;
}

/* ---- RGAS callbacks ---- */
static void gerg_initialize(Domain* domain, cxboolean vapor_phase,
                            char* name,
                            int (*msg)(const char*, ...),
                            void (*err)(const int, const char*, const char*))
{
    (void)domain; (void)vapor_phase; (void)name;
    gerg_msg = msg;
    gerg_err = err;
    gerg_ensure_init();
    if (gerg_msg)
        gerg_msg("GERG UDRGM: teqp GERG-2008 %d-species kernel %s\n",
                 GERG_NSPECIES,
                 gerg_initialized ? "initialized" : "UNAVAILABLE (ideal-gas)");
}

static double gerg_density(cell_t c, Thread* t, cxboolean vapor_phase,
                           double T, double P, double yi[])
{
    double rho = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized || T < 60.0 || P <= 0.0
        || gerg_rho_pT(T, P, yi, GERG_NSPECIES, &rho) != 0
        || !isfinite(rho) || rho <= 0.0)
    {
        /* ideal-gas guard (PR fallback hook point) */
        double mw = gerg_mw_from_y(yi);
        if (mw > 0.0 && T > 0.0)
            rho = P * mw / (GERG_RU * T);
        else
            rho = 1.0;
        gerg_fail_count++;
    }
    if (!(isfinite(rho) && rho > 0.0)) rho = 1.0;
    if (rho > 1.0e5) rho = 1.0e5;
    if (rho < 1.0e-3) rho = 1.0e-3;
    return rho;
}

static double gerg_enthalpy(cell_t c, Thread* t, double T, double P,
                            double rho, double yi[])
{
    double h = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized || !(isfinite(rho) && rho > 0.0)
        || gerg_h_rhoT(T, rho, yi, GERG_NSPECIES, &h) != 0 || !isfinite(h))
    {
        h = gerg_h_ideal(yi, T);
        if (!isfinite(h)) h = 0.0;
    }
    return h;
}

static double gerg_entropy(cell_t c, Thread* t, double T, double P,
                           double rho, double yi[])
{
    (void)rho;
    double s;
    gerg_ensure_yi(c, t, yi);
    s = gerg_s_ideal(yi, T, P);
    if (!isfinite(s)) s = 0.0;
    return s;
}

static double gerg_specific_heat(cell_t c, Thread* t, double T, double P,
                                 double rho, double yi[])
{
    double cp = 0.0, cv = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized || !(isfinite(rho) && rho > 0.0)
        || gerg_cp_cv(T, rho, yi, GERG_NSPECIES, &cp, &cv) != 0
        || !(isfinite(cp) && cp > 0.0))
    {
        cp = gerg_cp_ideal(yi, T);
        if (!(isfinite(cp) && cp > 0.0)) cp = 1000.0;
    }
    return cp;
}

static double gerg_molecular_weight(double yi[])
{
    double mw = gerg_mw_from_y(yi);
    return (isfinite(mw) && mw > 0.0) ? mw : 32.0;
}

static double gerg_speed_of_sound(cell_t c, Thread* t, double T, double P,
                                  double rho, double yi[])
{
    double a = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized || !(isfinite(rho) && rho > 0.0)
        || gerg_a_rhoT(T, rho, yi, GERG_NSPECIES, &a) != 0
        || !(isfinite(a) && a > 1.0))
    {
        double cp = gerg_cp_ideal(yi, T);
        double mw = gerg_mw_from_y(yi);
        double gamma = 1.4;
        if (mw > 0.0)
        {
            double cv = cp - GERG_RU / mw;
            if (cv > 1.0) gamma = cp / cv;
        }
        a = sqrt(gamma * GERG_RU * T / mw);
        if (!isfinite(a) || a < 1.0) a = 340.0;
    }
    return a;
}

static double gerg_viscosity(cell_t c, Thread* t, double T, double P,
                             double rho, double yi[])
{
    double mu = 0.0, k = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized
        || gerg_mu_k(T, rho, yi, GERG_NSPECIES, &mu, &k) != 0
        || !(isfinite(mu) && mu > 0.0))
        mu = 1.8e-5 * pow(T / 300.0, 0.7);
    return mu;
}

static double gerg_thermal_conductivity(cell_t c, Thread* t, double T,
                                        double P, double rho, double yi[])
{
    double mu = 0.0, k = 0.0;
    gerg_ensure_yi(c, t, yi);
    gerg_ensure_init();
    if (!gerg_initialized
        || gerg_mu_k(T, rho, yi, GERG_NSPECIES, &mu, &k) != 0
        || !(isfinite(k) && k > 0.0))
        k = gerg_viscosity(c, t, T, P, rho, yi) * gerg_cp_ideal(yi, T) * 1.3;
    return k;
}

static double gerg_rho_t(cell_t c, Thread* t, double T, double P,
                         double rho, double yi[])
{
    (void)rho;
    double dT = T * 1e-4 + 1.0;
    double r0 = gerg_density(c, t, TRUE, T, P, yi);
    double r1 = gerg_density(c, t, TRUE, T + dT, P, yi);
    return (r1 - r0) / dT;
}

static double gerg_rho_p(cell_t c, Thread* t, double T, double P,
                         double rho, double yi[])
{
    (void)rho;
    double dP = P * 1e-4 + 1.0;
    double r0 = gerg_density(c, t, TRUE, T, P, yi);
    double r1 = gerg_density(c, t, TRUE, T, P + dP, yi);
    return (r1 - r0) / dP;
}

static double gerg_enthalpy_t(cell_t c, Thread* t, double T, double P,
                              double rho, double yi[])
{
    return gerg_specific_heat(c, t, T, P, rho, yi);
}

static double gerg_enthalpy_p(cell_t c, Thread* t, double T, double P,
                              double rho, double yi[])
{
    (void)c; (void)t; (void)T; (void)P; (void)rho; (void)yi;
    return 0.0;
}

static double gerg_total_enthalpy(cell_t c, Thread* t, double T, double P,
                                  double rho, double yi[], double vel2[])
{
    return gerg_enthalpy(c, t, T, P, rho, yi);
}

static void gerg_flash(double t0, double p0, double z[],
                       double z0[], double xliq[], double xvap[])
{
    /* Full RGAS flash contract (13 fields), single-phase all-vapor:
       z0[0]=rho, [1]=rho_liq, [2]=rho_vap, [3]=quality(mole), [4]=e, [5]=h,
       [6]=s, [7]=cv, [8]=cp, [9]=a, [10]=mw_liq, [11]=mw_vap,
       [12]=vapor vol frac */
    double yi[GERG_NSPECIES];
    int i;
    for (i = 0; i < GERG_NSPECIES; ++i) yi[i] = z[i];
    if (z0)
    {
        double rho = gerg_density((cell_t)0, (Thread*)0, TRUE, t0, p0, yi);
        double h = gerg_h_ideal(yi, t0);
        double cp = gerg_cp_ideal(yi, t0);
        double mw = gerg_mw_from_y(yi);
        z0[0] = rho;
        z0[1] = rho;
        z0[2] = rho;
        z0[3] = 1.0;   /* quality: all vapor */
        z0[4] = h;
        z0[5] = h;
        z0[6] = 0.0;
        z0[7] = (cp > 0.0 && mw > 0.0) ? cp - GERG_RU / mw : 1000.0;
        z0[8] = cp;
        z0[9] = gerg_speed_of_sound((cell_t)0, (Thread*)0, t0, p0, rho, yi);
        z0[10] = mw;
        z0[11] = mw;
        z0[12] = 1.0;
    }
    if (xliq) for (i = 0; i < GERG_NSPECIES; ++i) xliq[i] = 0.0;
    if (xvap) for (i = 0; i < GERG_NSPECIES; ++i) xvap[i] = yi[i];
}

static void gerg_sat(double tp0, double z[], int index0, double z0[])
{
    (void)tp0; (void)z; (void)index0; (void)z0;
}

static void gerg_mass_mole_fractions(int index, double y[], double x[])
{
    int i;
    double moles = 0.0;
    for (i = 0; i < GERG_NSPECIES; ++i)
    {
        double yi = y[i] > 0.0 ? y[i] : 0.0;
        moles += yi / GERG_NASA7[i].mw;
    }
    for (i = 0; i < GERG_NSPECIES; ++i)
        x[i] = moles > 0.0 ? (y[i] / GERG_NASA7[i].mw) / moles : 0.0;
    (void)index;
}

static void gerg_set_sat_clip(const cxboolean sat_clip_flag) { (void)sat_clip_flag; }

static void gerg_set_error_handling(const cxboolean silent)
{
    gerg_silent = silent ? 1 : 0;
}

static int gerg_get_error(void) { return 0; }
static int gerg_get_non_fatal_error_count(void) { return gerg_fail_count; }
static int gerg_get_warning_count(int* codes, char* msg)
{
    (void)codes; (void)msg;
    return 0;
}

/* ---- material info: same numeric layout as the working PR UDRGM
   (first double is small so the string fields are NUL-terminated) ---- */
static const RGAS_Material_Info gerg_rgas_info = {
    50.0, 4000.0,      /* min/max temperature [K] */
    2000.0, 2.0e7,     /* min/max pressure [Pa] */
    436.0,             /* critical density [kg/m3] (O2) */
    5.043e6,           /* critical pressure [Pa] (O2) */
    154.6,             /* critical temperature [K] (O2) */
    50.0, 4000.0, 2000.0, 2.0e7, 2000.0,  /* eos range */
    50.0, 4000.0, 2000.0, 2.0e7,          /* tcond range */
    50.0, 4000.0, 2000.0, 2.0e7,          /* visc range */
};

UDF_EXPORT RGAS_Functions RealGasFunctionList = {
    gerg_initialize,
    gerg_density,
    gerg_enthalpy,
    gerg_entropy,
    gerg_specific_heat,
    gerg_molecular_weight,
    gerg_speed_of_sound,
    gerg_viscosity,
    gerg_thermal_conductivity,
    gerg_rho_t,
    gerg_rho_p,
    gerg_enthalpy_t,
    gerg_enthalpy_p,
    gerg_total_enthalpy,
    gerg_flash,
    gerg_sat,
    gerg_mass_mole_fractions,
    gerg_set_sat_clip,
    gerg_set_error_handling,
    gerg_get_error,
    gerg_get_non_fatal_error_count,
    gerg_get_warning_count,
    &gerg_rgas_info,        /* rgas_info */
};
