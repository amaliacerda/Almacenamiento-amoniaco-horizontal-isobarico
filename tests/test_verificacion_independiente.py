"""Verificación independiente del modelo (no depende de que los balances cierren).

error_masa y error_energia cierran casi por construcción (el esquema conservativo
telescopa), así que no detectan errores en P_j, A_f, k, U ni en el gradiente de la
interfase. Estos tests comparan contra fórmulas cerradas de geometria.py, cuadraturas
propias y CoolProp directo (sin el spline), mediante tank.diagnostico(t).
Geometría y condiciones: caso 3 (NH3, d_o = d_i, tapas planas) salvo indicación.
"""
import numpy as np
import pytest
import CoolProp.CoolProp as CP
from cryoevap.cryogens import Cryogen

from conftest import CASO3, crear_horizontal
from horizontal import GeometriaHorizontal, TablaVapor

T_DIAG = 20 * 3600.0  # s: vapor ya estratificado, T_j > T_sat en todas las celdas


def _crear(fluido, t_f=T_DIAG, **kw):
    c = CASO3
    kw.setdefault("eta_w", c["eta_w"])
    tk = crear_horizontal(fluido, c["P"], c["V"], c["d"], c["L"], c["LF"], c["T_air"], q=c["q"],
                          time_interval=3600.0, **kw)
    tk.evaporate(t_f)
    return tk


@pytest.fixture(scope="module")
def nh3():
    return _crear("Ammonia")


@pytest.fixture(scope="module")
def n2():
    return _crear("Nitrogen")


# ---------------------------------------------------------------- verificaciones
def verif_geometria(tk, d):
    """Devuelve dict de discrepancias relativas de la geometría por celda."""
    g = tk.geo
    h = d["h"]
    zf = d["z_f"]
    return dict(
        sum_P=abs(d["P_j"].sum() - g.A_dry(h)) / g.A_dry(h),
        sum_V=abs(d["V_j"].sum() - (tk.V - g.V_L(h))) / (tk.V - g.V_L(h)),
        # primitiva cerrada de P (no pasa por integral_P): P_j de cada celda
        P_celda=np.max(np.abs(d["P_j"] - (g._prim_P(zf[1:]) - g._prim_P(zf[:-1]))) / d["P_j"]),
    )


def Q_V_independiente(tk, d, n=40001):
    """U_V int_h^H P(z) (T_air - T(z)) dz con T lineal entre (h, T_sat), (z_c, T_j), techo = T_N.

    Integral en el ángulo z = R(1 - cos t): P dz = (2 L R_o + 4 R^2 sin^2 t) dt (planas),
    regular; regla trapezoidal fina propia.
    """
    g = tk.geo
    R = g.R
    t0, t1 = np.arccos((R - d["h"]) / R), np.pi
    t = np.linspace(t0, t1, n)
    z = R * (1 - np.cos(t))
    T = np.interp(z, np.concatenate([[d["h"]], d["z_c"], [g.H]]),
                  np.concatenate([[tk.cryogen.T_sat], d["T_j"], [d["T_j"][-1]]]))
    f = (2 * g.L * g.R_o + 4 * R ** 2 * np.sin(t) ** 2) * (tk.T_air - T)
    return tk.U_V * np.trapezoid(f, t)


def verif_Q(tk, d):
    g = tk.geo
    Q_L = tk.U_L * g.A_wet(d["h"]) * (tk.T_air - tk.cryogen.T_sat)
    Q_V = Q_V_independiente(tk, d)
    return dict(Q_L=abs(d["Q_L"] - Q_L) / Q_L, Q_V=abs(d["Q_V"] - Q_V) / Q_V,
                Q_in=abs(d["Q_in"] - (Q_L + Q_V)) / (Q_L + Q_V))


def verif_props(tk, d):
    """rho, h, k del vapor contra CoolProp directo en (P, T_j); dH_LV, rho_L, h_L en saturación."""
    P, fl = tk.cryogen.P, tk.cryogen.name
    rho = np.array([CP.PropsSI("D", "T", T, "P", P, fl) for T in d["T_j"]])
    hh = np.array([CP.PropsSI("H", "T", T, "P", P, fl) for T in d["T_j"]])
    k = np.array([CP.PropsSI("L", "T", T, "P", P, fl) for T in d["T_j"]])
    hV = CP.PropsSI("H", "P", P, "Q", 1, fl)
    hL = CP.PropsSI("H", "P", P, "Q", 0, fl)
    return dict(
        rho=np.max(np.abs(d["rho_j"] - rho) / rho),
        h=np.max(np.abs(d["h_j"] - hh) / np.abs(hh - hL)),  # relativo al calor latente (h puede ser ~0 o < 0)
        k=np.max(np.abs(tk.tabla.k(d["T_j"]) - k) / k),
        rho_L=abs(tk.cryogen.rho_L - CP.PropsSI("D", "P", P, "Q", 0, fl)) / tk.cryogen.rho_L,
        dH_LV=abs(tk.dH_LV - (hV - hL)) / (hV - hL),
    )


# ------------------------------------------------------------------------ tests
def test_geometria_por_celda(nh3):
    """sum(P_j) = A_dry(h), sum(V_j) = V - V_L(h), P_j = integral_P de cada celda (rel 1e-10)."""
    for t in (0.0, 5 * 3600.0, T_DIAG):
        r = verif_geometria(nh3, nh3.diagnostico(t))
        assert max(r.values()) < 1e-10, (t, r)


def test_Q_in_independiente(nh3):
    """Q_L cerrado (rel 1e-10) y Q_V contra cuadratura propia con T lineal (rel 1e-4).

    Medido: 3e-6 (NH3). La tolerancia 1e-4 refleja la discretización: el modelo toma T constante en cada celda
    (40 celdas, T_j en el centro) y la cuadratura propia usa T lineal; la diferencia es
    O(dz^2 T''), mayor cerca de la interfase donde el perfil se curva, y relativa a
    U_V*A_dry*(T_air - T) con T_air - T ~ 54 K y salto de T en el vapor de pocos K.
    """
    d = nh3.diagnostico(T_DIAG)
    r = verif_Q(nh3, d)
    assert r["Q_L"] < 1e-10, r
    assert r["Q_V"] < 1e-4, r
    assert r["Q_in"] < 1e-4, r
    assert d["Q_in"] == pytest.approx(d["Q_L"] + d["Q_V"] + d["Q_b"] + d["Q_roof"], rel=1e-12)


def test_Q_VL_independiente(nh3):
    """Q_VL = k(T_sat) A(h) dT/dz con gradiente de 2.º orden propio y k de CoolProp (rel 2e-3)."""
    d = nh3.diagnostico(T_DIAG)
    P = nh3.cryogen.P
    k_sat = CP.PropsSI("L", "P", P, "Q", 1, "Ammonia")
    d1, d2 = d["z_c"][0] - d["h"], d["z_c"][1] - d["h"]
    Ts = nh3.cryogen.T_sat
    dTdz = (d["T_j"][0] - Ts) * d2 / (d1 * (d2 - d1)) - (d["T_j"][1] - Ts) * d1 / (d2 * (d2 - d1))
    esperado = k_sat * nh3.geo.A(d["h"]) * dTdz
    assert d["Q_VL"] == pytest.approx(esperado, rel=2e-3)


@pytest.mark.parametrize("fx", ["nh3", "n2"])
def test_propiedades_vs_coolprop(fx, request):
    """rho_j y h_j (rel 1e-6), k (rel 1e-4), rho_L y dH_LV (rel 1e-6) contra CoolProp directo."""
    tk = request.getfixturevalue(fx)
    d = tk.diagnostico(T_DIAG)
    assert np.all(d["T_j"] > tk.cryogen.T_sat)
    r = verif_props(tk, d)
    assert r["rho"] < 1e-6, r
    assert r["h"] < 1e-6, r
    assert r["k"] < 1e-4, r
    assert r["rho_L"] < 1e-6, r
    assert r["dH_LV"] < 1e-6, r


def test_eta_w_1_todo_al_liquido():
    """eta_w = 1: el vapor queda a T_sat y BOG = Q_in/dH_LV (1 - rho_V_sat/rho_L) (rel 1e-3)."""
    tk = _crear("Ammonia", eta_w=1.0)
    d = tk.diagnostico(T_DIAG)
    Ts = tk.cryogen.T_sat
    assert np.max(np.abs(d["T_j"] - Ts)) < 1e-4  # K
    assert d["Q_Vv"] == pytest.approx(0.0, abs=1e-9)
    A_dry = tk.geo.A_dry(d["h"])
    assert d["Q_wi"] == pytest.approx(tk.U_V * A_dry * (tk.T_air - Ts), rel=1e-6)
    assert d["Q_VL"] == pytest.approx(0.0, abs=1e-6)
    esperado = d["Q_in"] / tk.dH_LV * (1 - tk.rho_V_sat / tk.cryogen.rho_L)
    assert d["BOG"] == pytest.approx(esperado, rel=1e-3)
    assert d["B_L"] == pytest.approx(d["Q_in"] / tk.dH_LV, rel=1e-6)


def test_vapor_calentado_eta_w_0():
    """eta_w = 0 (NH3, caso 3, U_V x 5 para vapor bien calentado): T_j crece con z,
    T_BOG > Tv_avg > T_sat y G_min > 0."""
    c = CASO3
    U = c["q"] / (c["T_air"] - 239.83)
    tk = crear_horizontal("Ammonia", c["P"], c["V"], c["d"], c["L"], c["LF"], c["T_air"], eta_w=0.0,
                          U=U, U_V=5 * U, time_interval=3600.0)
    tk.evaporate(100 * 3600.0)
    Ts = tk.cryogen.T_sat
    for i in (1, 5, 20, 100):
        T = tk.sol.y[1:, i]
        assert np.all(np.diff(T) > 0), "T_j no monótona en t = %d h" % i
        assert tk.data["T_BOG"][i] > tk.data["Tv_avg"][i] > Ts
    assert tk.data["Tv_avg"][-1] - Ts > 5.0  # realmente calentado
    assert tk.data["G_min"] > 0
    # propiedades contra CoolProp también con el vapor muy calentado (rango amplio de T)
    r = verif_props(tk, tk.diagnostico(100 * 3600.0))
    assert max(r["rho"], r["h"]) < 1e-6 and r["k"] < 1e-4, r


def _todas_las_verificaciones(tk, t):
    d = tk.diagnostico(t)
    g = verif_geometria(tk, d)
    q = verif_Q(tk, d)
    p = verif_props(tk, d)
    return max(g.values()), max(q["Q_L"], q["Q_V"], q["Q_in"]), max(p["rho"], p["h"], p["k"], p["rho_L"], p["dH_LV"])


def test_deteccion_P_x3(monkeypatch):
    """Con integral_P x 3 los balances siguen cerrando, pero las verificaciones fallan."""
    orig = GeometriaHorizontal.integral_P
    monkeypatch.setattr(GeometriaHorizontal, "integral_P", lambda self, a, b: 3 * orig(self, a, b))
    tk = _crear("Ammonia", t_f=5 * 3600.0)
    assert tk.data["error_masa"] < 1e-3 and tk.data["error_energia"] < 1e-3  # no lo ven
    g, q, _ = _todas_las_verificaciones(tk, 5 * 3600.0)
    assert g > 1e-2, "la verificación geométrica no detectó P x 3 (%.2e)" % g
    assert q > 1e-4, "la verificación de Q_in no detectó P x 3 (%.2e)" % q


def test_deteccion_k_x50(monkeypatch):
    """Con k del vapor x 50 los balances cierran, pero la comparación con CoolProp falla."""
    orig = TablaVapor.k
    monkeypatch.setattr(TablaVapor, "k", lambda self, T: 50 * orig(self, T))
    tk = _crear("Ammonia", t_f=5 * 3600.0)
    assert tk.data["error_masa"] < 1e-3 and tk.data["error_energia"] < 1e-3
    _, _, p = _todas_las_verificaciones(tk, 5 * 3600.0)
    assert p > 1e-2, "la verificación de propiedades no detectó k x 50 (%.2e)" % p
