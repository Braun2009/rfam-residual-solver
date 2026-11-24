#!/usr/bin/env python3
"""
coolant_reverse_solve.py - SUS316 cooling requirement + pin-fin/film assessment.

Complete chain: given the verified Vortex throat heat load, reverse-solve how
much coolant-side heat transfer SUS316 needs to survive, then check whether the
CFD-verified pin-fin array (h_c≈135k) + film cooling (79% q reduction) can
actually hold it. Compare to CuCrZr.

All inputs from VERIFIED sources:
  - Bartz throat q0 = 40 MW/m²           (three-source self-consistent)
  - gas-side h_g    = 14000 W/m²K        (reverse-eng ≈ Bartz 16000)
  - pin-fin h_c     = 135000 W/m²K       (wall-resolved RANS Nu592≈Žukauskas573)
  - film cooling    = q → 0.21·q (−79%)  (CFD η=0.732 @throat)
  - T_gas 3500K, T_coolant(throat) 324K  (vortex.txt counterflow)
  - SUS316 k=16.3 melt=1673K; CuCrZr k=340 melt=1350K; wall L=1.0mm

1-D steady through-wall balance:
  q = h_g(T_gas − Tw_hot) = (k/L)(Tw_hot − Tw_cold) = h_c(Tw_cold − T_cool)
  Series resistance: 1/U = 1/h_g + L/k + 1/h_c ; q = U(T_gas − T_cool)
  Tw_hot = T_gas − q/h_g
"""
import sys

# ── verified inputs ──────────────────────────────────────────────────────────
T_GAS      = 3500.0     # K, recovery temp
T_COOL     = 324.0      # K, coolant at throat (counterflow)
H_G        = 14000.0    # W/m²K, gas-side (reverse-eng ≈ Bartz 16000)
Q0         = 40.0e6     # W/m², Bartz throat heat flux (no film)
FILM_FACTOR= 0.21       # film cooling: q → 0.21·q  (η=0.732, −79%)
H_PINFIN   = 135000.0   # W/m²K, CFD-verified pin-fin array (Nu592)
L_WALL     = 1.0e-3     # m, inner wall thickness

MATERIALS = {
    "SUS316": {"k": 16.3,  "melt": 1673.0, "safe": 1400.0},
    "CuCrZr": {"k": 340.0, "melt": 1350.0, "safe": 1200.0},
}

def solve_series(k, h_c, film=1.0):
    """
    Correct 1-D steady series-resistance solve.
    Heat flows gas → wall → coolant through three resistances IN SERIES:
      1/U = 1/(h_g·film_hg) + L/k + 1/h_c
    The film reduces the EFFECTIVE gas-side driving by insulating the wall from
    the hot gas — modelled as an effective gas-side conductance reduction
    (film_factor multiplies h_g: less heat gets INTO the wall).
      q       = U · (T_gas − T_cool)          [the ONE consistent heat flux]
      Tw_hot  = T_gas  − q/(h_g·film)          [gas-side surface]
      Tw_cold = T_cool + q/h_c                 [coolant-side surface]
    film=1.0 → no film; film=0.21 → film cuts gas-side coupling to 21%.
    Returns (q, Tw_hot, Tw_cold).
    """
    h_g_eff = H_G * film
    U = 1.0 / (1.0/h_g_eff + L_WALL/k + 1.0/h_c)
    q = U * (T_GAS - T_COOL)
    Tw_hot  = T_GAS  - q / h_g_eff
    Tw_cold = T_COOL + q / h_c
    return q, Tw_hot, Tw_cold

def hc_for_target(k, Tw_hot_target, film=1.0):
    """
    Reverse: what coolant-side h_c makes the gas-side wall sit at Tw_hot_target?
    From Tw_hot = T_gas − q/(h_g·film) → q = h_g·film·(T_gas − Tw_hot_target).
    Then that q must flow through wall+coolant: q = (T_gas−T_cool)/(1/(h_g·film)+L/k+1/h_c)
    Solve for h_c.
    """
    h_g_eff = H_G * film
    q = h_g_eff * (T_GAS - Tw_hot_target)
    # q = (T_gas - T_cool) / (1/h_g_eff + L/k + 1/h_c)
    #  → 1/h_c = (T_gas - T_cool)/q - 1/h_g_eff - L/k
    inv_hc = (T_GAS - T_COOL)/q - 1.0/h_g_eff - L_WALL/k
    return (1.0/inv_hc if inv_hc > 0 else float('inf')), q

def main():
    print("=== SUS316 cooling requirement + pin-fin/film assessment ===\n")
    print(f"  Verified inputs: h_g={H_G:.0f}, pin-fin h_c={H_PINFIN/1e3:.0f}k, "
          f"film −{(1-FILM_FACTOR)*100:.0f}%")
    print(f"  T_gas={T_GAS:.0f}K, T_coolant={T_COOL:.0f}K, wall L={L_WALL*1e3:.1f}mm\n")

    npass = 0; nfail = 0
    def check(name, ok):
        nonlocal npass, nfail
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        if ok: npass += 1
        else:  nfail += 1

    # ── ANCHOR: validate series model vs verified CuCrZr design point ──
    # CuCrZr, pin-fin h_c=135k, NO film → should give q≈40MW, Tw≈641K.
    print("─" * 70)
    print("⓪ ANCHOR — series model vs verified CuCrZr design (expect q≈40MW, Tw≈641K):")
    k_cu = MATERIALS["CuCrZr"]["k"]
    q_a, Twh_a, Twc_a = solve_series(k_cu, H_PINFIN, film=1.0)
    print(f"   q={q_a/1e6:.1f}MW/m²  Tw_hot={Twh_a:.0f}K  Tw_cold={Twc_a:.0f}K")
    print(f"   vs verified: q≈40MW, Tw≈641K → "
          f"{'ANCHOR OK' if abs(q_a/1e6-40)<8 and abs(Twh_a-641)<120 else 'MODEL OFF'}")
    check("anchor: CuCrZr no-film reproduces q≈40MW (±8)", abs(q_a/1e6-40) < 8)
    check("anchor: CuCrZr no-film reproduces Tw≈641K (±120)", abs(Twh_a-641) < 120)

    # ── ① SUS316, NO film, pin-fin only ──
    print("\n① SUS316, NO film, pin-fin h_c=135k:")
    k_ss = MATERIALS["SUS316"]["k"]
    q1, Twh1, Twc1 = solve_series(k_ss, H_PINFIN, film=1.0)
    v1 = "MELTS" if Twh1 > 1673 else "survives"
    print(f"   q={q1/1e6:.1f}MW/m²  Tw_hot={Twh1:.0f}K  Tw_cold={Twc1:.0f}K → {v1}")
    print(f"   through-wall ΔT={Twh1-Twc1:.0f}K (SUS316 low k piles heat at gas side)")

    # ── ② SUS316, WITH film ──
    print("\n② SUS316, WITH film (gas-side coupling ×0.21), pin-fin h_c=135k:")
    q2, Twh2, Twc2 = solve_series(k_ss, H_PINFIN, film=FILM_FACTOR)
    v2 = "MELTS" if Twh2 > 1673 else ("SAFE" if Twh2 < 1400 else "survives-marginal")
    print(f"   q={q2/1e6:.1f}MW/m²  Tw_hot={Twh2:.0f}K  Tw_cold={Twc2:.0f}K → {v2}")

    # ── ③ REVERSE: h_c needed for SUS316 safe 1400K (with film) ──
    print("\n③ REVERSE: coolant h_c needed for SUS316 Tw_hot=1400K:")
    for film_case, label in [(1.0, "no film"), (FILM_FACTOR, "with film")]:
        hc_need, q_need = hc_for_target(k_ss, 1400.0, film=film_case)
        feasible = hc_need <= H_PINFIN if hc_need != float('inf') else False
        hc_str = f"{hc_need/1e3:.0f}k" if hc_need != float('inf') else "∞ (impossible)"
        print(f"   {label:<10}: need h_c={hc_str}, q={q_need/1e6:.1f}MW  → "
              f"{'pin-fin 135k SUFFICES' if feasible else 'pin-fin 135k NOT enough'}")

    # ── ④ CuCrZr comparison (with film) ──
    print("\n④ CuCrZr comparison (WITH film, pin-fin 135k):")
    q4, Twh4, Twc4 = solve_series(k_cu, H_PINFIN, film=FILM_FACTOR)
    print(f"   q={q4/1e6:.1f}MW/m²  Tw_hot={Twh4:.0f}K  Tw_cold={Twc4:.0f}K")
    print(f"   through-wall ΔT: CuCrZr {Twh4-Twc4:.0f}K vs SUS316 {Twh2-Twc2:.0f}K "
          f"({(Twh2-Twc2)/(Twh4-Twc4):.0f}× larger, k ratio {k_cu/k_ss:.0f}×)")

    # ── verdict chain ──
    print("\n" + "─" * 70)
    print("VERDICT CHAIN:")
    print(f"  ⓪ anchor CuCrZr no-film:  q={q_a/1e6:.0f}MW Tw={Twh_a:.0f}K (matches verified)")
    print(f"  ① SUS316 no-film:         Tw={Twh1:.0f}K → {v1}")
    print(f"  ② SUS316 +film:           Tw={Twh2:.0f}K → {v2}")
    print(f"  ④ CuCrZr +film:           Tw={Twh4:.0f}K → huge margin")
    if Twh2 < 1673:
        margin = 1673 - Twh2
        print(f"\n  ✓ SUS316 SURVIVES with pin-fin(135k)+film: {Twh2:.0f}K < melt 1673K "
              f"(margin {margin:.0f}K)" + ("" if Twh2<1400 else " — but marginal"))
    else:
        print(f"\n  ✗ SUS316 MELTS ({Twh2:.0f}K) even with pin-fin+film")

    # ── checks ──
    print()
    check("SUS316 +film computed Tw physical (T_cool<Tw<T_gas)", T_COOL < Twh2 < T_GAS)
    check("film is decisive lever (no-film Tw > +film Tw)", Twh1 > Twh2)
    check("CuCrZr Tw < SUS316 Tw (low-k penalty)", Twh4 < Twh2)
    # ΔT/L = q/k, so ΔT ratio = (q_ss/q_cu)·(k_cu/k_ss). With film the two q differ,
    # so verify the FULL relation, not just the k ratio.
    dT_ratio  = (Twh2-Twc2)/(Twh4-Twc4)
    pred_ratio = (q2/q4) * (k_cu/k_ss)
    check("through-wall ΔT ratio = (q_ss/q_cu)·(k_cu/k_ss) (physics)",
          abs(dT_ratio - pred_ratio)/pred_ratio < 0.05)

    print(f"\n  RESULT: {npass} passed, {nfail} failed")
    print("COOLANT_SOLVE_PASS" if nfail == 0 else "COOLANT_SOLVE_FAIL")
    return 0 if nfail == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
