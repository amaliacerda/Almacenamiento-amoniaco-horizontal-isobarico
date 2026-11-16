"""Barrido de llenados, validación, sensibilidad y comparación de fluidos (Proyecto 12).

Uso (desde la raíz del repositorio):
    python3 scripts/barrido.py                 # corre todo (~2-3 min con 4 procesos)
    python3 scripts/barrido.py --usar-cache    # rehace solo las figuras con scripts/.cache_barrido.pkl

Solo usa la API pública de src/horizontal (HorizontalTank, diagnostico) y no modifica
el modelo ni los tests. Escribe figuras (PDF vectorial + PNG) y CSV en figs/.

Supuestos del caso 3 (PROVISIONALES, AGENTS.md): amoníaco con las dimensiones del
caso 2 (d_i = d_o = 1.676 m, L = 2.217 m, V = 4.89 m^3, tapas planas), flujo de calor
equivalente 0.5 W/m^2 a T_air = 293.15 K (U = 0.5/(293.15 - T_sat)), eta_w = 0.97,
P = 101325 Pa. T_air = 288.15 K (sur) y 303.15 K (norte) son valores ILUSTRATIVOS.
Vertical del mismo volumen: misma d_i, altura H = V/(pi R^2) (= L = 2.217 m con tapas
planas), mismo U; como variante, techo con Q_roof = U A_T (T_air - T_sat) (cota superior)
para igualar el área de envolvente total del horizontal.
"""
import argparse
import csv
import os
import pickle
import sys
import warnings
from multiprocessing import Pool

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(RAIZ, "figs")
CACHE = os.path.join(RAIZ, "scripts", ".cache_barrido.pkl")
RUTA_WANG = os.path.join(RAIZ, "data", "Data_wang_merida.xlsx")

# ---- Parámetros (AGENTS.md)
P_NH3 = 101325.0
V_TOT, D, L = 4.89, 1.676, 2.217
ETA_W = 0.97
T_AIR = 293.15
Q_EQ = 0.5                      # W/m^2 a T_air = 293.15 K (equivalente del caso 2)
T_AIR_SUR, T_AIR_NORTE = 288.15, 303.15
LFS = np.round(np.arange(0.1, 0.91, 0.1), 1)
T_F = 200 * 3600.0
H_BOR = 24                      # BOR de referencia: instante t = 24 h
TIEMPOS_PERFIL_H = [1, 10, 50, 200]
NOTA3 = ("Caso 3 PROVISIONAL: amoníaco con dimensiones del caso 2 (d = 1.676 m, L = 2.217 m, "
         "V = 4.89 m³), U ≡ 0.5 W/m² a 293.15 K, η_w = 0.97, P = 1 atm")


# ------------------------------------------------------------------ simulación
def simular(cfg):
    """Corre una simulación y devuelve un dict con series y resúmenes (picklable)."""
    from cryoevap.cryogens import Cryogen
    from horizontal import HorizontalTank

    cr = Cryogen(name=cfg["fluido"])
    cr.set_coolprops(cfg["P"])
    U = cfg.get("U")
    if U is None:
        U = cfg["q"] / (cfg["T_ref"] - cr.T_sat)
    vertical = cfg.get("orient") == "vertical"
    if vertical:
        tk = HorizontalTank(D, D, cfg["V"], LF=cfg["LF"], orientacion="vertical")
    else:
        tk = HorizontalTank(cfg["d"], cfg["d"], cfg["V"], LF=cfg["LF"], L=cfg["L"])
    tk.cryogen = cr
    q_roof = U * tk.geo.A_T * (cfg["T_air"] - cr.T_sat) if cfg.get("techo") else 0.0
    tk.set_HeatTransProps(U_L=U, U_V=U, T_air=cfg["T_air"],
                          Q_b_fixed=None if vertical else 0, Q_roof=q_roof, eta_w=cfg["eta_w"])
    tk.time_interval = cfg.get("dt", 3600.0)
    if cfg.get("n_celdas"):
        tk.n_celdas = cfg["n_celdas"]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        tk.evaporate(cfg["t_f"])
    d = tk.data
    out = dict(cfg=cfg, U=U, T_sat=cr.T_sat, avisos=[str(x.message)[:80] for x in w
                                                      if issubclass(x.category, RuntimeWarning)],
               t=d["Time"].copy(), BOG=d["BOG"].copy(), B_L=d["B_L"].copy(), LF=d["LF"].copy(),
               m_L=d["m_L"].copy(), m_V=d["m_V"].copy(), Tv_avg=d["Tv_avg"].copy(),
               T_BOG=d["T_BOG"].copy(), Q_in=d["Q_in"].copy(), Q_L=d["Q_L"].copy(),
               Q_VL=d["Q_VL"].copy(), Q_V=d["Q_V"].copy(), h=d["h"].copy(),
               error_masa=d["error_masa"], error_energia=d["error_energia"], G_min=d["G_min"])
    if cfg.get("perfiles"):
        perf = {}
        for th in cfg["perfiles"]:
            dg = tk.diagnostico(th * 3600.0)
            perf[th] = dict(z_c=dg["z_c"], T_j=dg["T_j"], h=dg["h"])
        out["perfiles"] = perf
    return out


def bor(r, t_h=H_BOR):
    """BOR instantáneo [%/día] = BOG / m_L a t_h horas."""
    i = int(np.argmin(np.abs(r["t"] / 3600 - t_h)))
    return 100 * r["BOG"][i] * 86400 / r["m_L"][i]


def cfg_base(fluido="Ammonia", LF=0.3, **kw):
    c = dict(fluido=fluido, P=P_NH3, V=V_TOT, d=D, L=L, LF=float(LF), T_air=T_AIR, T_ref=T_AIR,
             q=Q_EQ, eta_w=ETA_W, t_f=T_F, orient="horizontal", techo=False, dt=3600.0)
    c.update(kw)
    return c


def construir_casos():
    """Lista de (clave, cfg). Las claves identifican cada corrida."""
    from cryoevap.cryogens import Cryogen
    cr = Cryogen(name="Ammonia")
    cr.set_coolprops(P_NH3)
    U0 = Q_EQ / (T_AIR - cr.T_sat)  # U base del amoníaco, también usado para N2 (misma aislación)
    casos = []
    for fl in ("Ammonia", "Nitrogen"):
        for LF in LFS:
            for orient, techo in (("horizontal", False), ("vertical", False), ("vertical", True)):
                tag = "H" if orient == "horizontal" else ("V" if not techo else "Vt")
                perf = TIEMPOS_PERFIL_H if (fl == "Ammonia" and orient == "horizontal") else None
                casos.append((("barrido", fl, tag, float(LF)),
                              cfg_base(fl, LF, orient=orient, techo=techo, U=U0, perfiles=perf)))
    # Sensibilidad (amoníaco, horizontal). U fijo en T_air variable (misma aislación).
    for LF in LFS:
        for f in (0.5, 2.0):
            casos.append((("sens", "U", f, float(LF)), cfg_base(LF=LF, U=U0 * f)))
        for e in (0.0, 0.5):
            casos.append((("sens", "eta", e, float(LF)), cfg_base(LF=LF, U=U0, eta_w=e)))
        for Ta in (T_AIR_SUR, T_AIR_NORTE):
            casos.append((("sens", "Tair", Ta, float(LF)), cfg_base(LF=LF, U=U0, T_air=Ta)))
    # Convergencia de malla (amoníaco, LF = 0.3, horizontal)
    for n in (20, 40, 80, 160):
        casos.append((("malla", n), cfg_base(LF=0.3, U=U0, n_celdas=n)))
    # Validación caso 2 (hidrógeno)
    casos.append((("caso2",), dict(fluido="Hydrogen", P=101300.0, V=4.89, d=1.676, L=2.217, LF=0.30,
                                   T_air=293.15, T_ref=293.15, q=0.5, eta_w=0.97, t_f=T_F,
                                   orient="horizontal", techo=False, dt=600.0)))
    return casos, U0


def _worker(item):
    return item[0], simular(item[1])


def correr_todo():
    casos, U0 = construir_casos()
    with Pool(min(4, os.cpu_count() or 1)) as pool:
        res = dict(pool.imap_unordered(_worker, casos))
    res["_U0"] = U0
    with open(CACHE, "wb") as f:
        pickle.dump(res, f)
    return res


# ------------------------------------------------------------------ utilidades
def estilo():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": "#dddddd",
        "grid.linewidth": 0.5, "lines.linewidth": 1.8, "figure.dpi": 100,
        "pdf.fonttype": 42, "savefig.bbox": "tight", "axes.axisbelow": True,
    })


AZUL, NARANJA, AQUA, VIOLETA, ROJO, GRIS = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948", "#52514e"


def guardar(fig, nombre, cabecera, columnas, nota=None):
    """Guarda PDF + PNG + CSV (columnas: dict nombre -> arreglo; largos iguales o se rellena con vacío)."""
    if nota:
        fig.text(0.01, -0.02, nota, fontsize=7, color=GRIS, ha="left", va="top", wrap=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, nombre + "." + ext), dpi=200)
    n = max(len(v) for v in columnas.values())
    with open(os.path.join(FIGS, nombre + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# " + cabecera])
        w.writerow(list(columnas))
        for i in range(n):
            w.writerow(["%.8g" % v[i] if i < len(v) and v[i] == v[i] else "" for v in columnas.values()])
    import matplotlib.pyplot as plt
    plt.close(fig)
    print("  figura:", nombre)


# ------------------------------------------------------------------ figuras
def fig_geometria():
    import matplotlib.pyplot as plt
    from horizontal.geometria import GeometriaHorizontal, GeometriaVertical
    gh = GeometriaHorizontal(D, D, L)
    gv = GeometriaVertical(D, D, V_TOT)
    lf = np.linspace(0.005, 0.995, 199)
    hh = np.array([gh.h_desde_V(x * V_TOT) for x in lf])
    hv = lf * gv.H
    A_tot_h = gh.area_total()
    A_wet_h = np.array([gh.A_wet(h) for h in hh])
    # vertical: manto mojado + fondo (siempre mojado); total = manto + fondo + techo
    A_tot_v = np.pi * D * gv.H + 2 * gv.A_T
    A_wet_v = np.array([gv.A_wet(h) for h in hv]) + gv.A_T
    A_if_h = np.array([gh.A(h) for h in hh])
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.2))
    ax[0].plot(lf, A_wet_h / A_tot_h, color=AZUL, label="Horizontal")
    ax[0].plot(lf, A_wet_v / A_tot_v, color=NARANJA, ls="--", label="Vertical (manto + fondo)")
    ax[0].plot(lf, lf, color=GRIS, lw=0.8, ls=":", label="Referencia: fracción = LF")
    ax[0].set_xlabel("Llenado LF [-]")
    ax[0].set_ylabel("Fracción de pared mojada A_wet / A_total [-]")
    ax[0].set_title("(a) Pared mojada")
    ax[0].legend(loc="upper left")
    ax[1].plot(lf, A_if_h, color=AZUL, label="Horizontal, A(h)")
    ax[1].axhline(gv.A_T, color=NARANJA, ls="--", label="Vertical, π d²/4 (constante)")
    ax[1].set_xlabel("Llenado LF [-]")
    ax[1].set_ylabel("Área de la interfase [m²]")
    ax[1].set_title("(b) Interfase líquido-vapor")
    ax[1].set_ylim(0, 4.6)
    ax[1].legend(loc="lower center")
    fig.suptitle("Geometría: d = 1.676 m, V = 4.89 m³; vertical con H = V/(πR²) = %.3f m" % gv.H, fontsize=9.5)
    fig.tight_layout()
    guardar(fig, "fig01_geometria_pared_mojada",
            "LF; fraccion pared mojada horizontal (planas) y vertical (manto+fondo, total con techo); area interfase [m2]",
            dict(LF=lf, h_horiz_m=hh, frac_mojada_horiz=A_wet_h / A_tot_h, frac_mojada_vert=A_wet_v / A_tot_v,
                 A_interfase_horiz_m2=A_if_h, A_interfase_vert_m2=np.full_like(lf, gv.A_T)),
            nota="Resultado geométrico exacto (geometria.py), independiente de la física y de la malla.")


def leer_wang():
    import openpyxl
    wb = openpyxl.load_workbook(RUTA_WANG, data_only=True)

    def cols(hoja, *idx):
        filas = [r for r in wb[hoja].iter_rows(min_row=2, values_only=True) if r[idx[0]] is not None]
        return [np.array([r[i] for r in filas], dtype=float) for i in idx]
    return dict(T=cols("12a_T", 4, 5), LF=cols("12b_LF", 0, 1), dML=cols("12d_dMLdt", 2, 3),
                BOG=cols("12e_BOG", 2, 3))


def fig_validacion(res):
    import matplotlib.pyplot as plt
    r = res[("caso2",)]
    w = leer_wang()
    th = r["t"] / 3600
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 5.6))
    a = ax[0, 0]
    a.plot(th, r["BOG"] * 3600, color=AZUL, label="Modelo")
    a.plot(*w["BOG"], "o", color=NARANJA, ms=4, label="Wang y Mérida (digitalizado)")
    a.set_ylabel("BOG [kg/h]"); a.set_title("(a) BOG"); a.legend(loc="lower right")
    a = ax[0, 1]
    a.plot(th, -r["B_L"] * 3600, color=AZUL, label="Modelo")
    a.plot(w["dML"][0], -w["dML"][1], "o", color=NARANJA, ms=4, label="Wang y Mérida")
    a.set_ylabel("Pérdida de masa líquida −dm_L/dt [kg/h]"); a.set_title("(b) Evaporación del líquido")
    a.legend(loc="lower right")
    a = ax[1, 0]
    a.plot(th, 100 * r["LF"], color=AZUL, label="Modelo")
    a.plot(*w["LF"], "o", color=NARANJA, ms=4, label="Wang y Mérida")
    a.set_ylabel("Llenado LF [%]"); a.set_title("(c) Llenado"); a.legend(loc="lower left")
    a = ax[1, 1]
    a.plot(th, r["Tv_avg"], color=AZUL, label="Modelo, T media del vapor")
    a.plot(th, r["T_BOG"], color=AZUL, ls="--", lw=1.2, label="Modelo, T del vapor en el techo")
    a.plot(*w["T"], "o", color=NARANJA, ms=4, label="Wang y Mérida, T vapor")
    a.set_ylabel("Temperatura del vapor [K]"); a.set_title("(d) Temperatura del vapor")
    a.legend(loc="lower right", fontsize=7)
    for a in ax.ravel():
        a.set_xlabel("Tiempo [h]"); a.set_xlim(0, 200)
    ax[0, 0].set_ylim(0, 0.2); ax[0, 1].set_ylim(0, 0.2)   # misma escala (a) y (b)
    fig.suptitle("Validación con HIDRÓGENO (caso 2, horizontal, LF₀ = 0.30, 0.5 W/m², η_w = 0.97)", fontsize=9.5)
    fig.tight_layout()
    ib = np.interp(w["BOG"][0], th, r["BOG"] * 3600)
    guardar(fig, "fig02_validacion_caso2_h2",
            "caso 2 hidrogeno: modelo (cada 600 s) y puntos digitalizados Wang y Merida (hojas 12e, 12d, 12b, 12a)",
            dict(t_modelo_h=th, BOG_modelo_kgh=r["BOG"] * 3600, perdida_liquido_modelo_kgh=-r["B_L"] * 3600,
                 LF_modelo_pct=100 * r["LF"], Tv_avg_modelo_K=r["Tv_avg"], T_techo_modelo_K=r["T_BOG"],
                 t_BOG_wang_h=w["BOG"][0], BOG_wang_kgh=w["BOG"][1], BOG_modelo_en_t_wang_kgh=ib,
                 t_dML_wang_h=w["dML"][0], perdida_wang_kgh=-w["dML"][1],
                 t_LF_wang_h=w["LF"][0], LF_wang_pct=w["LF"][1],
                 t_Tv_wang_h=w["T"][0], Tv_wang_K=w["T"][1]),
            nota="La validación es con H2; la aplicación es con NH3. El modelo queda ~26 % bajo Wang y Mérida en BOG: "
                 "no se ajustó ningún parámetro.")
    dif = (ib - w["BOG"][1]) / w["BOG"][1]
    return dict(dif_media=float(np.mean(np.abs(dif[w["BOG"][0] >= 15]))), dif_t0=float(dif[0]),
                bog_modelo_fin=float(r["BOG"][-1] * 3600), bog_wang_fin=float(w["BOG"][1][-1]),
                perdida0_modelo=float(-r["B_L"][0] * 3600), perdida0_wang=float(-w["dML"][1][0]),
                perdida_fin_modelo=float(-r["B_L"][-1] * 3600), perdida_fin_wang=float(-w["dML"][1][-1]),
                LF_fin_modelo=float(100 * r["LF"][-1]), LF_fin_wang=float(w["LF"][1][-1]),
                Tv_fin_modelo=float(r["Tv_avg"][-1]), Tv_fin_wang=float(w["T"][1][-1]))


def serie_bor(res, fl, tag, t_h=H_BOR):
    return np.array([bor(res[("barrido", fl, tag, float(LF))], t_h) for LF in LFS])


def fig_barrido(res):
    import matplotlib.pyplot as plt
    bh, bv, bvt = (serie_bor(res, "Ammonia", t) for t in ("H", "V", "Vt"))
    b200 = {t: serie_bor(res, "Ammonia", t, 200) for t in ("H", "V", "Vt")}
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.4))
    a = ax[0]
    a.plot(LFS, bh, "o-", color=AZUL, label="Horizontal")
    a.plot(LFS, bv, "s--", color=NARANJA, label="Vertical, techo adiabático (convención CryoEvap)")
    a.plot(LFS, bvt, "^:", color=AQUA, label="Vertical, con calor por el techo (cota sup.)")
    a.set_xlabel("Llenado inicial LF [-]"); a.set_ylabel("BOR a t = 24 h [%/día]")
    a.set_title("(a) BOR de amoníaco"); a.legend(loc="upper right", fontsize=7)
    a = ax[1]
    a.plot(LFS, bh / bv, "o-", color=NARANJA, label="Horizontal / vertical (techo adiabático)")
    a.plot(LFS, bh / bvt, "^:", color=AQUA, label="Horizontal / vertical (con techo)")
    a.axhline(1, color=GRIS, lw=0.8)
    a.set_xlabel("Llenado inicial LF [-]"); a.set_ylabel("Razón de BOR [-]")
    a.set_title("(b) Razón horizontal / vertical"); a.legend(loc="best", fontsize=7)
    fig.suptitle("Barrido de llenados, mismo volumen (4.89 m³) y misma aislación (U = %.4f W/m²/K)" % res["_U0"], fontsize=9.5)
    fig.tight_layout()
    guardar(fig, "fig03_bor_vs_lf_nh3",
            "BOR [%/dia] = BOG/m_L a t=24 h y a t=200 h; NH3; H horizontal, V vertical techo adiabatico, Vt vertical con Q_roof",
            dict(LF=LFS, BOR24_horiz=bh, BOR24_vert=bv, BOR24_vert_techo=bvt,
                 BOR200_horiz=b200["H"], BOR200_vert=b200["V"], BOR200_vert_techo=b200["Vt"],
                 razon_H_V=bh / bv, razon_H_Vt=bh / bvt),
            nota=NOTA3 + ". Vertical: misma d, H = V/(πR²) = 2.217 m.")
    return dict(bh=bh, bv=bv, bvt=bvt, b200=b200)


def fig_bor_tiempo(res):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    cmap = mpl.colormaps["Blues"]
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.3), sharey=True)
    cols = {}
    for k, (tag, tit) in enumerate((("H", "(a) Horizontal"), ("V", "(b) Vertical (techo adiabático)"))):
        for i, LF in enumerate(LFS):
            r = res[("barrido", "Ammonia", tag, float(LF))]
            y = 100 * r["BOG"] * 86400 / r["m_L"]
            ax[k].plot(r["t"] / 3600, y, color=cmap(0.3 + 0.7 * i / 8), lw=1.4,
                       label="LF₀ = %.1f" % LF)
            cols["%s_LF%.1f" % (tag, LF)] = y
        ax[k].set_xlabel("Tiempo [h]"); ax[k].set_title(tit); ax[k].set_yscale("log")
    ax[0].set_ylabel("BOR instantáneo BOG/m_L [%/día]")
    ax[1].legend(ncol=1, fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle("Transitorio del BOR (amoníaco, mismas condiciones que la fig. 3)", fontsize=9.5)
    fig.tight_layout()
    cols = dict(t_h=res[("barrido", "Ammonia", "H", 0.1)]["t"] / 3600, **cols)
    guardar(fig, "fig04_bor_vs_tiempo_nh3", "BOR(t) [%/dia] NH3, columnas H_LFx.x horizontal y V_LFx.x vertical",
            cols, nota=NOTA3 + ".")


def fig_perfiles(res):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(8.2, 3.4), sharey=True, sharex=True)
    estilos = [(AZUL, "-"), (NARANJA, "--"), (AQUA, "-."), (VIOLETA, ":")]
    cols = {}
    Ts = None
    for a, LF in zip(ax, (0.3, 0.5, 0.8)):
        r = res[("barrido", "Ammonia", "H", LF)]
        Ts = r["T_sat"]
        for th, (c, ls) in zip(TIEMPOS_PERFIL_H, estilos):
            p = r["perfiles"][th]
            z = np.concatenate([[p["h"]], p["z_c"]])
            T = np.concatenate([[Ts], p["T_j"]])
            a.plot(z, T, color=c, ls=ls, label="t = %d h" % th)
            cols["z_LF%.1f_t%dh_m" % (LF, th)] = z
            cols["T_LF%.1f_t%dh_K" % (LF, th)] = T
        a.axvline(D, color=GRIS, lw=0.8)
        a.text(D - 0.03, 0.03, "techo →", transform=a.get_xaxis_transform(), ha="right", va="bottom", fontsize=7, color=GRIS)
        a.set_title("LF₀ = %.1f" % LF); a.set_xlabel("Altura z desde el fondo [m]")
        a.set_xlim(0, D)
    ax[0].set_ylabel("Temperatura del vapor [K]")
    ax[0].legend(loc="upper left", fontsize=7)
    fig.suptitle("Perfiles de temperatura del vapor (amoníaco, horizontal, T_sat = %.1f K, T_air = 293.15 K)" % Ts, fontsize=9.5)
    fig.tight_layout()
    guardar(fig, "fig05_perfiles_Tv_nh3", "T de cada celda de vapor (z_c, T_j) a varios tiempos, NH3 horizontal; primer punto = interfase (T_sat)",
            cols, nota=NOTA3 + ". El perfil parte de la interfase en T_sat. "
                       "Con 20 a 160 celdas la T del techo cambia < 0.02 K (efecto numérico despreciable).")


def fig_sensibilidad(res):
    import matplotlib.pyplot as plt
    U0 = res["_U0"]
    base = serie_bor(res, "Ammonia", "H")

    def serie(par, v):
        return np.array([bor(res[("sens", par, v, float(LF))]) for LF in LFS])
    grupos = [
        ("U", "(a) Coeficiente global U", [(0.5, "U × 0.5 (%.4f)" % (0.5 * U0), AQUA, "--"), (1.0, "U base (%.4f)" % U0, AZUL, "-"),
                                          (2.0, "U × 2 (%.4f)" % (2 * U0), NARANJA, ":")], "W/m²/K"),
        ("eta", "(b) Fracción η_w", [(0.0, "η_w = 0", AQUA, "--"), (0.5, "η_w = 0.5", NARANJA, ":"), (0.97, "η_w = 0.97 (base)", AZUL, "-")], ""),
        ("Tair", "(c) T_air (ilustrativa)", [(T_AIR_SUR, "288.15 K (sur)", AQUA, "--"), (T_AIR, "293.15 K (base)", AZUL, "-"),
                                             (T_AIR_NORTE, "303.15 K (norte)", NARANJA, ":")], ""),
    ]
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 3.4), sharey=True)
    cols = dict(LF=LFS)
    rel = {}
    for a, (par, tit, items, _) in zip(ax, grupos):
        for v, lab, c, ls in items:
            y = base if (par == "U" and v == 1.0) or (par == "eta" and v == 0.97) or (par == "Tair" and v == T_AIR) else serie(par, v)
            a.plot(LFS, y, color=c, ls=ls, marker="o", ms=3, label=lab)
            cols["BOR24_%s_%s" % (par, v)] = y
            rel[(par, v)] = y / base
        a.set_title(tit); a.set_xlabel("Llenado inicial LF [-]"); a.set_yscale("log")
        a.legend(fontsize=7, loc="upper right")
    ax[0].set_ylabel("BOR a t = 24 h [%/día]")
    fig.suptitle("Sensibilidad (amoníaco, horizontal)", fontsize=9.5)
    fig.tight_layout()
    guardar(fig, "fig06_sensibilidad_nh3", "BOR a 24 h [%/dia], NH3 horizontal; U en W/m2/K, eta_w, T_air en K (U fijo al variar T_air)",
            cols, nota=NOTA3 + ". T_air = 288.15 K (sur) y 303.15 K (norte) son valores ILUSTRATIVOS; U fijo al variar T_air.")
    return rel


def fig_fluidos(res):
    import matplotlib.pyplot as plt
    b = {(fl, t): serie_bor(res, fl, t) for fl in ("Ammonia", "Nitrogen") for t in ("H", "V", "Vt")}
    g = {(fl, t): np.array([res[("barrido", fl, t, float(LF))]["BOG"][H_BOR] * 3600 for LF in LFS])
         for fl in ("Ammonia", "Nitrogen") for t in ("H", "V", "Vt")}
    fig, ax = plt.subplots(1, 3, figsize=(8.6, 3.4))
    for fl, c, nom in (("Ammonia", AZUL, "NH₃"), ("Nitrogen", NARANJA, "N₂")):
        ax[0].plot(LFS, b[(fl, "H")], "o-", color=c, label=nom + " horizontal")
        ax[0].plot(LFS, b[(fl, "V")], "s--", color=c, lw=1.2, ms=3.5, label=nom + " vertical")
        ax[1].plot(LFS, g[(fl, "H")], "o-", color=c, label=nom + " horizontal")
        ax[1].plot(LFS, g[(fl, "V")], "s--", color=c, lw=1.2, ms=3.5, label=nom + " vertical")
        ax[2].plot(LFS, b[(fl, "H")] / b[(fl, "V")], "o-", color=c, label=nom + " / techo adiabático")
        ax[2].plot(LFS, b[(fl, "H")] / b[(fl, "Vt")], "^:", color=c, label=nom + " / con techo")
    ax[0].set_yscale("log"); ax[1].set_yscale("log")
    ax[0].set_ylabel("BOR a t = 24 h [%/día]"); ax[0].set_title("(a) BOR")
    ax[1].set_ylabel("BOG a t = 24 h [kg/h]"); ax[1].set_title("(b) BOG")
    ax[2].set_ylabel("Razón BOR horizontal / vertical [-]"); ax[2].set_title("(c) Razón H/V")
    ax[2].axhline(1, color=GRIS, lw=0.8)
    for a in ax:
        a.set_xlabel("Llenado inicial LF [-]")
    ax[0].legend(fontsize=6.5, loc="upper right"); ax[2].legend(fontsize=6.5, loc="best")
    fig.suptitle("Nitrógeno frente a amoníaco en el mismo estanque (misma U = %.4f W/m²/K, T_air = 293.15 K, η_w = 0.97)" % res["_U0"], fontsize=9)
    fig.tight_layout()
    cols = dict(LF=LFS)
    for (fl, t), y in b.items():
        cols["BOR24_%s_%s" % (fl, t)] = y
    for (fl, t), y in g.items():
        cols["BOG24_kgh_%s_%s" % (fl, t)] = y
    guardar(fig, "fig07_comparacion_n2_nh3", "BOR [%/dia] y BOG [kg/h] a 24 h; H horizontal, V vertical techo adiabatico, Vt vertical con techo",
            cols, nota=NOTA3 + ". Misma U absoluta (misma aislación): el N₂ recibe ~4 veces más calor por su mayor ΔT. "
                       "No hay validación propia con N₂ horizontal (solo regresión vertical, caso 1).")
    return b, g


def tablas_csv(res, val):
    """CSV de errores de conservación y de convergencia de malla."""
    with open(os.path.join(FIGS, "tabla_conservacion.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["fluido", "orientacion", "LF0", "error_masa", "error_energia", "G_min_kg_s", "avisos_G_neg"])
        for k, r in sorted((k, r) for k, r in res.items() if k[0] == "barrido"):
            w.writerow([k[1], k[2], k[3], "%.3e" % r["error_masa"], "%.3e" % r["error_energia"], "%.3e" % r["G_min"],
                        len(r["avisos"])])
    filas = []
    with open(os.path.join(FIGS, "tabla_malla.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n_celdas", "BOG_24h_kg_h", "BOG_200h_kg_h", "Tv_avg_200h_K", "T_techo_200h_K", "error_masa", "error_energia"])
        for n in (20, 40, 80, 160):
            r = res[("malla", n)]
            fila = [n, r["BOG"][24] * 3600, r["BOG"][-1] * 3600, r["Tv_avg"][-1], r["T_BOG"][-1], r["error_masa"], r["error_energia"]]
            filas.append(fila)
            w.writerow(["%.8g" % x for x in fila])
    return filas


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--usar-cache", action="store_true")
    args = ap.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    if args.usar_cache and os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            res = pickle.load(f)
    else:
        print("Corriendo simulaciones ...")
        res = correr_todo()
    estilo()
    print("Figuras:")
    fig_geometria()
    val = fig_validacion(res)
    sw = fig_barrido(res)
    fig_bor_tiempo(res)
    fig_perfiles(res)
    rel = fig_sensibilidad(res)
    b, g = fig_fluidos(res)
    malla = tablas_csv(res, val)
    # Salida de texto con los números clave
    np.set_printoptions(precision=4, linewidth=150)
    print("LF     ", LFS)
    print("BOR24 H ", sw["bh"]); print("BOR24 V ", sw["bv"]); print("BOR24 Vt", sw["bvt"])
    print("H/V    ", sw["bh"] / sw["bv"]); print("H/Vt   ", sw["bh"] / sw["bvt"])
    print("validacion:", val)
    for k, v in rel.items():
        print("sens", k, np.round(v, 3))
    for k, v in b.items():
        print("BOR24", k, v)
    print("malla:", malla)
    mx = max(r["error_masa"] for k, r in res.items() if k != "_U0")
    me = max(r["error_energia"] for k, r in res.items() if k != "_U0")
    gm = min(r["G_min"] for k, r in res.items() if k != "_U0")
    nav = sum(len(r["avisos"]) for k, r in res.items() if k != "_U0")
    print("max err masa %.2e, max err energia %.2e, G_min global %.3e, avisos %d" % (mx, me, gm, nav))


if __name__ == "__main__":
    sys.path.insert(0, os.path.join(RAIZ, "src"))
    main()
