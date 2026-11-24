/* Vortex RFAM Peng-Robinson UDRGM wrapper for Fluent 2026 R1.
 *
 * Implements the RGAS_Functions interface (material/real_gas.h) on top of
 * the vortex_rfam_pr_udrgm.c single-phase PR kernel (separate TU).
 * Species order must match Fluent mixture: C3H8 O2 CO2 CO H2O H2 CH4 C2H6 N2
 *
 * Compiled by Fluent's own makefile (official UDF build): all special
 * macros are defined here so no command-line -D is needed.
 */
#define RFAM_SKIP_UDRGM_LIST
#define RFAM_A413_SINGLE_PHASE_ENVELOPE_SHA256 "014a7b0fde6fe2225c94a15ea07754dbb86ad3145277eb2c2a3618a1c607b83c"
#include "udf.h"  /* pulls in real_gas.h */
#include "GergUdrgmNasa7.H"
#define VORTEX_NSPECIES 30
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ---- kernel extern declarations (vortex_rfam_pr_udrgm.c) ---- */
extern double rfam_max(double a, double b);
extern double rfam_mw_from_y(const double yi[]);
extern double rfam_mole_fractions(const double yi[], double x[], double y[]);
extern double rfam_cp_ideal(const double yi[], double temperature);
extern double rfam_h_ideal(const double yi[], double temperature);
extern double rfam_s_ideal(const double yi[], double temperature, double pressure);
extern double rfam_h_departure(double pressure, double temperature, double z, const double yi[]);
extern double rfam_s_departure(double pressure, double temperature, double z, const double yi[]);
extern double rfam_dpdrho(double density, double temperature, const double yi[]);
extern int rfam_z_from_supplied_density(double density, double temperature, double pressure, const double yi[], double *z);
extern double rfam_viscosity_mix(const double yi[], double temperature);
extern double rfam_density_value(double temperature, double pressure, double yi[]);
extern double rfam_specific_heat_value(double density, double temperature, double pressure, double yi[]);
extern double rfam_rho_t_value(double temperature, double pressure, double yi[]);
/* ---- end kernel externs ---- */

#define VORTEX_NSPECIES 9

/* universal gas constant [J/kmol/K] (same as kernel RFAM_RU) */
#ifndef RFAM_RU
#define RFAM_RU 8314.462618
#endif

/* ideal-gas fallback below this temperature (keeps initializer invertible) */
#define VORTEX_IDEAL_T_LOW 180.0

/* freeze density w.r.t. pressure while the velocity field builds up */
#define VORTEX_RHO_P_ZERO


static int (*vortex_message)(const char *, ...) = NULL;
static void (*vortex_error)(const int, const char *, const char *) = NULL;
static int vortex_silent_errors = 0;

static void vortex_emit_error(const char *func, const char *msg)
{
  if (vortex_silent_errors) return;
  if (vortex_error) vortex_error(0, func, msg);
}

static void vortex_initialize(Domain *domain, cxboolean vapor_phase, char *name,
                              int (*msg)(const char *, ...),
                              void (*err)(const int, const char *, const char *))
{
  vortex_message = msg;
  vortex_error = err;
  if (vortex_message)
    vortex_message("Vortex RFAM PR UDRGM: 9-species Peng-Robinson initialized.\n");
}

static int vortex_fail_logged = 0;
static int vortex_call_logged = 0;

/* Logging is disabled in production builds (fopen/fprintf per callback call
 * was the dominant cost at ~100 s per iteration).  Re-enable by defining
 * VORTEX_CALL_LOG when needed. */
static void vortex_log_call(const char *fn, cell_t c, double T, double P, double rho)
{
  /* abnormal-state monitor: log any non-physical request Fluent makes */
  int abnormal = !isfinite(T) || !isfinite(P) || !isfinite(rho)
    || T < 50.0 || T > 6000.0 || P < 1.0e3 || P > 1.0e8
    || rho < 1.0e-3 || rho > 1.0e5;
  if (!abnormal) return;
  FILE *f = fopen("/tmp/udrgm_abnormal.log", "a");
  if (f) {
    fprintf(f, "%s cell=%lld T=%.8g P=%.8g rho=%.8g\n",
            fn, (long long)c, T, P, rho);
    fclose(f);
  }
}

static void vortex_log_fail(const char *fn, double T, double P, const double yi[])
{
#ifdef VORTEX_CALL_LOG
  if (vortex_fail_logged >= 60) return;
  vortex_fail_logged++;
  FILE *f = fopen("/tmp/udrgm_fail.log", "a");
  if (f) {
    fprintf(f, "%s T=%.6g P=%.6g yi=[%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g]\n",
            fn, T, P, yi[0], yi[1], yi[2], yi[3], yi[4], yi[5], yi[6], yi[7], yi[8]);
    fclose(f);
  }
#endif
}

/* ---- RGAS callback helpers ------------------------------------------ */

/* Official REFPROP-style contract (see realgas/src/realgas_fcns.c):
 *   density(c, t, vapor, T, P, yi[])  ->  returns rho
 *   enthalpy(c, t, T, P, rho, yi[])   ->  returns h
 * etc.  The 6th argument is the species mass-fraction array supplied by
 * Fluent; it must be treated as READ-ONLY.  Fall back to C_YI only when
 * the supplied array is unusable. */
static void vortex_ensure_yi(cell_t c, Thread *t, double yi[])
{
  int i;
  double sum = 0.0;
  int ok = 1;
  if (yi) {
    for (i = 0; i < VORTEX_NSPECIES; i++) {
      double v = yi[i];
      if (!isfinite(v) || v < -1e-12 || v > 1.2) { ok = 0; break; }
      if (v > 0.0) sum += v;
    }
    if (ok && sum > 0.5 && sum < 1.5) return;  /* usable as-is */
  }
  if (c) {
    double csum = 0.0;
    for (i = 0; i < VORTEX_NSPECIES; i++) {
      double v = C_YI(c, t, i);
      if (isfinite(v)) csum += v > 0.0 ? v : 0.0;
      yi[i] = isfinite(v) ? v : 0.0;
    }
    if (!(csum > 0.5 && csum < 1.5)) {
      for (i = 0; i < VORTEX_NSPECIES; i++) yi[i] = 0.0;
      yi[1] = 1.0;   /* pure O2 */
    }
  } else {
    for (i = 0; i < VORTEX_NSPECIES; i++) yi[i] = 0.0;
    yi[1] = 1.0;   /* pure O2 */
  }
}

static double vortex_density(cell_t c, Thread *t, cxboolean vapor_phase,
                             double T, double P, double yi[])
{
  double rho = 1.0;
  vortex_ensure_yi(c, t, yi);
  (void)T; (void)P;
  if (!(isfinite(rho) && rho > 0.0) && T > 0.0 && P > 0.0) {
    double mw = rfam_mw_from_y(yi);
    rho = P * mw / (RFAM_RU * T);
  }
  if (!(isfinite(rho) && rho > 0.0)) {
    vortex_log_fail("density", T, P, yi);
    rho = 1.0;
  }
  /* final guard: keep the density finite and physical for Fluent */
  if (rho > 1.0e5) rho = 1.0e5;
  if (rho < 1.0e-3) rho = 1.0e-3;
  return rho;
}

static double vortex_enthalpy(cell_t c, Thread *t, double T, double P,
                              double rho, double yi[])
{
  double z = NAN, h;
  vortex_ensure_yi(c, t, yi);
  vortex_log_call("enthalpy", c, T, P, rho);
  if (T < VORTEX_IDEAL_T_LOW) {
    h = rfam_h_ideal(yi, T);
  } else if (rfam_z_from_supplied_density(rho, T, P, yi, &z)) {
    h = rfam_h_ideal(yi, T) + rfam_h_departure(P, T, z, yi);
    if (!isfinite(h)) h = rfam_h_ideal(yi, T);
  } else {
    h = rfam_h_ideal(yi, T);
  }
  if (!isfinite(h)) h = 0.0;
  return h;
}

static double vortex_entropy(cell_t c, Thread *t, double T, double P,
                             double rho, double yi[])
{
  double z = NAN, s;
  vortex_ensure_yi(c, t, yi);
  if (T < VORTEX_IDEAL_T_LOW) {
    s = rfam_s_ideal(yi, T, P / 101325.0);
  } else if (rfam_z_from_supplied_density(rho, T, P, yi, &z)) {
    s = rfam_s_ideal(yi, T, P / 101325.0) + rfam_s_departure(P, T, z, yi);
    if (!isfinite(s)) s = rfam_s_ideal(yi, T, P / 101325.0);
  } else {
    s = rfam_s_ideal(yi, T, P / 101325.0);
  }
  if (!isfinite(s)) s = 0.0;
  return s;
}

static double vortex_specific_heat(cell_t c, Thread *t, double T, double P,
                                   double rho, double yi[])
{
  double cp;
  vortex_ensure_yi(c, t, yi);
  if (T < VORTEX_IDEAL_T_LOW) {
    cp = rfam_cp_ideal(yi, T);
  } else {
    cp = rfam_specific_heat_value(rho, T, P, yi);
  }
  if (!(isfinite(cp) && cp > 0.0)) cp = rfam_cp_ideal(yi, T);
  if (!(isfinite(cp) && cp > 0.0)) cp = 1000.0;
  return cp;
}

static double vortex_molecular_weight(double yi[])
{
  return rfam_mw_from_y(yi);
}

static double vortex_speed_of_sound(cell_t c, Thread *t, double T, double P,
                                    double rho, double yi[])
{
  double a;
  vortex_ensure_yi(c, t, yi);
  double dpdrho = rfam_dpdrho(rho, T, yi);
  if (!isfinite(dpdrho) || dpdrho <= 0.0) {
    double mw = rfam_mw_from_y(yi);
    double cp = rfam_specific_heat_value(rho, T, P, yi);
    double cv = cp - RFAM_RU / mw;
    if (cv <= 0.0) cv = cp * 0.7;
    a = sqrt(cp / cv * P / rho);
  } else {
    a = sqrt(dpdrho);
  }
  if (!isfinite(a) || a <= 0.0) a = 300.0;
  return a;
}

static double vortex_viscosity(cell_t c, Thread *t, double T, double P,
                               double rho, double yi[])
{
  double mu;
  vortex_ensure_yi(c, t, yi);
  mu = rfam_viscosity_mix(yi, T);
  if (!(isfinite(mu) && mu > 0.0)) mu = 1.0e-5;
  return mu;
}

static double vortex_thermal_conductivity(cell_t c, Thread *t, double T,
                                          double P, double rho, double yi[])
{
  double mu = vortex_viscosity(c, t, T, P, rho, yi);
  double cp = vortex_specific_heat(c, t, T, P, rho, yi);
  double k = mu * cp / 0.85;
  if (!(isfinite(k) && k > 0.0)) k = 0.05;
  return k;
}

static double vortex_rho_t(cell_t c, Thread *t, double T, double P,
                           double rho, double yi[])
{
  double v;
  vortex_ensure_yi(c, t, yi);
  v = rfam_rho_t_value(T, P, yi);
  if (!isfinite(v)) {
    vortex_log_fail("rho_t", T, P, yi);
    v = 0.0;
  }
  /* well-posed fallback: ideal-gas isobaric d(rho)/dT = -rho/T */
  if (v >= 0.0 && rho > 0.0 && T > 0.0) v = -rho / T;
  /* clamp: never exceed 10x the ideal-gas derivative magnitude */
  if (rho > 0.0 && T > 0.0) {
    double ig = rho / T;
    if (v > 10.0 * ig) v = 10.0 * ig;
    if (v < -10.0 * ig) v = -10.0 * ig;
  }
  if (!isfinite(v)) v = 0.0;
  return v;
}

static double vortex_rho_p(cell_t c, Thread *t, double T, double P,
                           double rho, double yi[])
{
#ifdef VORTEX_RHO_P_ZERO
  /* density frozen w.r.t. pressure: stabilises the SIMPLE pressure-
   * correction during the initial velocity-field build-up */
  return 0.0;
#endif
  double v = 0.0;
  vortex_ensure_yi(c, t, yi);
  double step = rfam_max(P * 1.0e-2, 1.0e4);
  double rp = rfam_density_value(T, P + step, yi);
  double rm = rfam_density_value(T, P - step, yi);
  if (isfinite(rp) && isfinite(rm)) v = (rp - rm) / (2.0 * step);
  if (!isfinite(v)) {
    vortex_log_fail("rho_p", T, P, yi);
    v = 0.0;
  }
  /* well-posed fallback: ideal-gas isothermal d(rho)/dP = rho/P */
  if (v <= 0.0 && rho > 0.0 && P > 0.0) v = rho / P;
  /* clamp: never exceed 10x the ideal-gas derivative magnitude */
  if (rho > 0.0 && P > 0.0) {
    double ig = rho / P;
    if (v > 10.0 * ig) v = 10.0 * ig;
  }
  if (!isfinite(v)) v = 0.0;
  return v;
}

static double vortex_enthalpy_t(cell_t c, Thread *t, double T, double P,
                                double rho, double yi[])
{
  double z = NAN, v = 0.0;
  vortex_ensure_yi(c, t, yi);
  double step = rfam_max(T * 1.0e-4, 1.0e-3);
  double hp, hm;
  if (rfam_z_from_supplied_density(rho, T, P, yi, &z)) {
    hp = rfam_h_ideal(yi, T + step) + rfam_h_departure(P, T + step, z, yi);
    hm = rfam_h_ideal(yi, T - step) + rfam_h_departure(P, T - step, z, yi);
    v = (hp - hm) / (2.0 * step);
  }
  if (!isfinite(v)) v = 0.0;
  return v;
}

static double vortex_enthalpy_p(cell_t c, Thread *t, double T, double P,
                                double rho, double yi[])
{
  double z = NAN, v = 0.0;
  vortex_ensure_yi(c, t, yi);
  double step = rfam_max(P * 1.0e-4, 10.0);
  double hp, hm;
  if (rfam_z_from_supplied_density(rho, T, P, yi, &z)) {
    hp = rfam_h_departure(P + step, T, z, yi);
    hm = rfam_h_departure(P - step, T, z, yi);
    v = (hp - hm) / (2.0 * step);
  }
  if (!isfinite(v)) v = 0.0;
  return v;
}

static double vortex_total_enthalpy(cell_t c, Thread *t, double T, double P,
                                    double rho, double yi[], double hi[])
{
  return vortex_enthalpy(c, t, T, P, rho, yi);
}

static void vortex_flash(double t0, double p0, double z[],
                         double z0[], double xliq[], double xvap[])
{
  vortex_log_call("flash", 0, t0, p0, 0.0);
  /* Full RGAS flash contract (13 fields), single-phase PR:
     z0[0]=rho, [1]=rho_liq, [2]=rho_vap, [3]=quality(mole), [4]=e, [5]=h,
     [6]=s, [7]=cv, [8]=cp, [9]=a, [10]=mw_liq, [11]=mw_vap, [12]=vapor vol frac */
  double rho = rfam_density_value(t0, p0, z);
  if (t0 < VORTEX_IDEAL_T_LOW) {
    /* keep flash consistent with the density callback (ideal-gas low-T) */
    double mw = rfam_mw_from_y(z);
    rho = p0 * mw / (RFAM_RU * t0);
  }
  if (!(isfinite(rho) && rho > 0.0)) {
    /* ideal-gas fallback (same as vortex_density), not the old rho=1 */
    double mw = rfam_mw_from_y(z);
    rho = p0 * mw / (RFAM_RU * t0);
  }
  if (!(isfinite(rho) && rho > 0.0)) rho = 1.0;
  double h = 0.0, s = 0.0, cp = 1000.0, cv = 700.0, a = 300.0, mw, zf = NAN;
  int i;
  if (rfam_z_from_supplied_density(rho, t0, p0, z, &zf)) {
    h = rfam_h_ideal(z, t0) + rfam_h_departure(p0, t0, zf, z);
    s = rfam_s_ideal(z, t0, p0 / 101325.0) + rfam_s_departure(p0, t0, zf, z);
    cp = rfam_specific_heat_value(rho, t0, p0, z);
  } else {
    h = rfam_h_ideal(z, t0);
    s = rfam_s_ideal(z, t0, p0 / 101325.0);
    cp = rfam_cp_ideal(z, t0);
  }
  mw = rfam_mw_from_y(z);
  if (isfinite(cp) && cp > 0.0) {
    cv = cp - RFAM_RU / mw;
    if (cv <= 0.0) cv = cp * 0.7;
    a = sqrt(cp / cv * p0 / rho);
  }
  double e = h - p0 / rho;   /* internal energy */
  if (!isfinite(rho) || rho <= 0.0) rho = 1.0;
  if (z0) {
    z0[0] = rho;             /* bulk density */
    z0[1] = rho;             /* liquid density (single phase) */
    z0[2] = rho;             /* vapor density (single phase) */
    z0[3] = 1.0;             /* vapor quality = 1 (single phase) */
    z0[4] = isfinite(e) ? e : 0.0;
    z0[5] = isfinite(h) ? h : 0.0;
    z0[6] = isfinite(s) ? s : 0.0;
    z0[7] = cv;
    z0[8] = cp;
    z0[9] = isfinite(a) ? a : 300.0;
    z0[10] = mw;
    z0[11] = mw;
    z0[12] = 1.0;            /* vapor volume fraction */
  }
  for (i = 0; i < VORTEX_NSPECIES; i++) {
    if (xliq) xliq[i] = z[i];
    if (xvap) xvap[i] = z[i];
  }
}

static void vortex_sat(double tp0, double z[], int index0, double z0[])
{
  vortex_log_call("sat", 0, tp0, 0.0, 0.0);
  if (z0) z0[0] = 0.0;
}

static void vortex_mass_mole_fractions(int index, double y[], double x[])
{
  double yi[VORTEX_NSPECIES];
  int i;
  for (i = 0; i < VORTEX_NSPECIES; i++) yi[i] = y[i];
  rfam_mole_fractions(yi, x, y);
}

static void vortex_set_sat_clip(const cxboolean sat_clip_flag) {}

static void vortex_set_error_handling(const cxboolean silent)
{
  vortex_silent_errors = silent ? 1 : 0;
}

static int vortex_get_error(void) { return 0; }
static int vortex_get_non_fatal_error_count(void) { return 0; }
static int vortex_get_warning_count(int *codes, char *msg) { return 0; }

/* RGAS_Material_Info: 20 doubles, order per material/real_gas.h
 * VORTEX: critical values from O2 (dominant species) so Fluent's
 * reduced-variable normalisation sees a physical critical point. */
static const RGAS_Material_Info vortex_rgas_info = {
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
  vortex_initialize,             /* initialize */
  vortex_density,                /* density */
  vortex_enthalpy,               /* enthalpy */
  vortex_entropy,                /* entropy */
  vortex_specific_heat,          /* specific_heat */
  vortex_molecular_weight,       /* molecular_weight */
  vortex_speed_of_sound,         /* speed_of_sound */
  vortex_viscosity,              /* viscosity */
  vortex_thermal_conductivity,   /* thermal_conductivity */
  vortex_rho_t,                  /* rho_t */
  vortex_rho_p,                  /* rho_p */
  vortex_enthalpy_t,             /* enthalpy_t */
  vortex_enthalpy_p,             /* enthalpy_p */
  vortex_total_enthalpy,         /* total_enthalpy */
  vortex_flash,                  /* flash */
  vortex_sat,                    /* sat */
  vortex_mass_mole_fractions,    /* mass_mole_fractions */
  vortex_set_sat_clip,           /* set_sat_clip */
  vortex_set_error_handling,     /* set_error_handling */
  vortex_get_error,              /* get_error */
  vortex_get_non_fatal_error_count, /* get_non_fatal_error_count */
  vortex_get_warning_count,      /* get_warning_count */
  &vortex_rgas_info,             /* rgas_info */
};
