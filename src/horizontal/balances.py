"""Balances globales de masa y energía de una corrida (docs/ecuaciones.md §5).

Masa:     E_m(t) = [m_L(t) + m_V(t)] - [m_L(0) + m_V(0)] + int_0^t B dt
Energía:  E_E(t) = H(t) - H(0) + int_0^t B h_V(T_techo) dt - int_0^t Q_in dt
          H = m_L h_L + sum_j rho_j h(T_j) V_j     [J]
          Q_in = Q_L + Q_b + U_V int P (T_air - T) dz + Q_roof  (calor del ambiente) [W]
Errores relativos (adimensionales):
          error_masa    = |E_m(t_f)| / int_0^t_f B dt        (masa evaporada/venteada)
          error_energia = |E_E(t_f)| / int_0^t_f Q_in dt      (calor total recibido)

IMPORTANTE: estos errores miden la consistencia del esquema numérico y del
integrador en el tiempo, NO son una validación física. En modo 'conservativo' el
esquema semidiscreto conserva masa y energía exactamente por construcción (BOG = G
en el techo, cp = dh/dT), así que los errores solo reflejan la tolerancia de
solve_ivp y la cuadratura (típicamente 1e-8). En modo 'uniforme' o con propiedades
'promedio' el error de energía sí mide la inconsistencia de esas aproximaciones.
La validez física se juzga contra CryoEvap (caso 1) y Wang y Mérida (caso 2).

Las integrales en el tiempo se hacen con cuadratura de Gauss-Legendre
(tank.cuadratura_gauss puntos) sobre cada paso interno del integrador, usando la
salida densa de solve_ivp; B y h_V(T_techo) se evalúan con el mismo modelo (no por
diferencias finitas). T_techo es la temperatura de la celda superior, que es la que
sale por el venteo con el esquema upwind.
"""
import numpy as np


def _integrandos(tank, y):
    """Devuelve (B [kg/s], B h_techo [W], Q_in [W], m_total [kg], H [J]) para el estado y."""
    _, d = tank._modelo(y)
    T = y[1:]
    h_techo = tank.tabla.h(T[-1])
    Q_in = d["Q_L"] + d["Q_b"] + d["Q_V"] + d["Q_roof"]
    m_L = tank.cryogen.rho_L * y[0]
    m_tot = m_L + d["m_V"]
    H = m_L * tank.cryogen.h_L + np.sum(d["rho"] * d["hT"] * d["V_j"])
    return d["BOG"], d["BOG"] * h_techo, Q_in, m_tot, H


def calcular_balances(tank):
    """Calcula los errores de conservación y los guarda en tank.data.

    Claves escritas:
        data['error_masa'], data['error_energia']        errores relativos finales [-]
        data['error_masa_t'], data['error_energia_t']    errores relativos en cada
                                                         tiempo de data['Time'] (nan en t = 0)
        data['m_BOG_acum'] [kg], data['Q_in_acum'] [J]   integrales acumuladas
    Devuelve (error_masa, error_energia).
    """
    sol = tank.sol_completa
    xg, wg = np.polynomial.legendre.leggauss(tank.cuadratura_gauss)
    t_pasos = sol.t
    t_dat = tank.sol.t
    # Nodos de corte: pasos internos + tiempos guardados (para acumular en cada uno)
    cortes = np.unique(np.concatenate([t_pasos, t_dat]))
    I = np.zeros((len(cortes), 3))  # acumulados de B, B h_techo, Q_in
    for i in range(len(cortes) - 1):
        a, b = cortes[i], cortes[i + 1]
        tg = 0.5 * (b - a) * xg + 0.5 * (a + b)
        Y = sol.sol(tg)
        s = np.zeros(3)
        for q in range(len(tg)):
            B, BH, Q, _, _ = _integrandos(tank, Y[:, q])
            s += wg[q] * np.array([B, BH, Q])
        I[i + 1] = I[i] + 0.5 * (b - a) * s
    idx = np.searchsorted(cortes, t_dat)
    I_dat = I[idx]

    _, _, _, m0, H0 = _integrandos(tank, tank.sol.y[:, 0])
    em, ee = [], []
    for k in range(len(t_dat)):
        _, _, _, m, H = _integrandos(tank, tank.sol.y[:, k])
        E_m = (m - m0) + I_dat[k, 0]
        E_E = (H - H0) + I_dat[k, 1] - I_dat[k, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            em.append(abs(E_m) / I_dat[k, 0] if I_dat[k, 0] > 0 else np.nan)
            ee.append(abs(E_E) / I_dat[k, 2] if I_dat[k, 2] > 0 else np.nan)
    tank.data["error_masa_t"] = np.array(em)
    tank.data["error_energia_t"] = np.array(ee)
    tank.data["error_masa"] = float(em[-1])
    tank.data["error_energia"] = float(ee[-1])
    tank.data["m_BOG_acum"] = I_dat[:, 0]
    tank.data["Q_in_acum"] = I_dat[:, 2]
    return tank.data["error_masa"], tank.data["error_energia"]
