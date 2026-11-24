/* test_gerg_sanity.c - standalone kernel sanity (no Fluent headers). */
#include <stdio.h>
#include <math.h>

extern int gerg_init(const char* fluiddata, const char* const* names, int n);
extern int gerg_rho_pT(double T, double p, const double* y, int n, double* rho);
extern int gerg_p_rhoT(double T, double rho, const double* y, int n, double* p);
extern int gerg_h_rhoT(double T, double rho, const double* y, int n, double* h);
extern int gerg_a_rhoT(double T, double rho, const double* y, int n, double* a);
extern int gerg_cp_cv(double T, double rho, const double* y, int n,
                      double* cp, double* cv);

static const char* NAMES[] = {
    "N2","H","O2","OH","O","H2","H2O","HO2","H2O2","CO","CO2","HCO","CH3",
    "CH4","CH2O","CH2S","C2H4","CH3O","C2H5","C2H6","C2H2","C2H3","CH2CHO",
    "CH2CO","CH3CHO","C3H4","C3H6","C3H8","IC3H7","NC3H7"
};
#define N 30

static double y[N];

static void pure(int idx) { for (int i = 0; i < N; ++i) y[i] = 0.0; y[idx] = 1.0; }
static int name_idx(const char* s) {
    for (int i = 0; i < N; ++i) if (strcmp(NAMES[i], s) == 0) return i;
    return -1;
}
#include <string.h>

int main(void)
{
    if (gerg_init("/root/teqp/dev/fluiddata", NAMES, N) != 0)
    {
        printf("GERG_SANITY_FAIL init\n");
        return 1;
    }
    int fails = 0;
    double rho, p, h, a, cp, cv;

    /* 1) pure O2 130 K / 6 MPa: subcooled liquid, rho ~ 909 kg/m3 */
    pure(name_idx("O2"));
    if (gerg_rho_pT(130.0, 6e6, y, N, &rho) != 0) { printf("FAIL O2 rho\n"); fails++; }
    else printf("O2 130K/6MPa: rho=%.2f kg/m3 (ref ~909) %s\n", rho,
                fabs(rho - 909.0) / 909.0 < 0.05 ? "OK" : "CHECK");

    /* 2) pure C3H8 450 K / 4.5 MPa: supercritical, rho ~ 65.5 */
    pure(name_idx("C3H8"));
    if (gerg_rho_pT(450.0, 4.5e6, y, N, &rho) != 0) { printf("FAIL C3H8 rho\n"); fails++; }
    else printf("C3H8 450K/4.5MPa: rho=%.2f kg/m3 (ref ~65.5) %s\n", rho,
                fabs(rho - 65.5) / 65.5 < 0.10 ? "OK" : "CHECK");

    /* 3) pure O2 130K sound speed + enthalpy sanity */
    pure(name_idx("O2"));
    gerg_rho_pT(130.0, 6e6, y, N, &rho);
    gerg_a_rhoT(130.0, rho, y, N, &a);
    gerg_h_rhoT(130.0, rho, y, N, &h);
    printf("O2 130K liquid: a=%.1f m/s h=%.0f J/kg (liquid a ~800-1100)\n", a, h);

    /* 4) combustion products at 3000 K / 4 MPa (equilibrium-ish mix) */
    for (int i = 0; i < N; ++i) y[i] = 0.0;
    y[name_idx("H2O")] = 0.35; y[name_idx("CO2")] = 0.40; y[name_idx("CO")] = 0.10;
    y[name_idx("H2")] = 0.03; y[name_idx("OH")] = 0.03; y[name_idx("O2")] = 0.05;
    y[name_idx("H")] = 0.01; y[name_idx("N2")] = 0.03;
    if (gerg_rho_pT(3000.0, 4e6, y, N, &rho) != 0) { printf("FAIL prod rho\n"); fails++; }
    else {
        gerg_a_rhoT(3000.0, rho, y, N, &a);
        gerg_cp_cv(3000.0, rho, y, N, &cp, &cv);
        printf("products 3000K/4MPa: rho=%.3f kg/m3 a=%.0f m/s cp=%.1f cv=%.1f\n",
               rho, a, cp, cv);
    }
    printf(fails == 0 ? "GERG_SANITY_PASS\n" : "GERG_SANITY_FAIL\n");
    return fails == 0 ? 0 : 1;
}
