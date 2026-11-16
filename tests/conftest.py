"""Utilidades y fixtures comunes de los tests del Proyecto 12.

Parámetros base de AGENTS.md (casos 1, 2 y 3). Los fixtures con scope de módulo
evitan repetir corridas dentro de un mismo archivo de tests.
"""
import os

import pytest
from cryoevap.cryogens import Cryogen
from horizontal import HorizontalTank

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_WANG = os.path.join(RAIZ, "data", "Data_wang_merida.xlsx")

# Caso 2 (Wang y Mérida): geometría y condiciones
CASO2 = dict(P=101300.0, V=4.89, d=1.676, L=2.217, LF=0.30, q=0.5, eta_w=0.97,
             T_air=293.15, t_f=200 * 3600.0)
# Caso 3 (mientras no haya datos definitivos): misma geometría y U equivalente del caso 2
CASO3 = dict(P=101325.0, V=4.89, d=1.676, L=2.217, LF=0.30, q=0.5, eta_w=0.97,
             T_air=293.15)


def crear_horizontal(fluido, P, V, d, L, LF, T_air, eta_w, q=None, U=None, U_V=None,
                     n_celdas=None, time_interval=600.0):
    """HorizontalTank horizontal (tapas planas, d_o = d_i) listo para evaporate.

    U_L = U_V = U; si no se da U, U = q / (T_air - T_sat) (flujo de calor q [W/m^2]).
    U_V explícito (p. ej. 0) reemplaza solo el coeficiente del vapor.
    """
    tk = HorizontalTank(d, d, V, LF=LF, L=L)
    cr = Cryogen(name=fluido)
    cr.set_coolprops(P)
    tk.cryogen = cr
    if U is None:
        U = q / (T_air - cr.T_sat)
    tk.set_HeatTransProps(U_L=U, U_V=U if U_V is None else U_V, T_air=T_air,
                          Q_b_fixed=0, Q_roof=0, eta_w=eta_w)
    if n_celdas is not None:
        tk.n_celdas = n_celdas
    tk.time_interval = time_interval
    return tk


@pytest.fixture(scope="module")
def tanque_h2_caso2():
    """Caso 2 (hidrógeno, 200 h)."""
    c = CASO2
    tk = crear_horizontal("Hydrogen", c["P"], c["V"], c["d"], c["L"], c["LF"], c["T_air"],
                          c["eta_w"], q=c["q"])
    tk.evaporate(c["t_f"])
    return tk


@pytest.fixture(scope="module")
def tanque_nh3_caso3():
    """Caso 3 de prueba (amoníaco, 200 h)."""
    c = CASO3
    tk = crear_horizontal("Ammonia", c["P"], c["V"], c["d"], c["L"], c["LF"], c["T_air"],
                          c["eta_w"], q=c["q"])
    tk.evaporate(200 * 3600.0)
    return tk


@pytest.fixture(scope="module")
def tanque_n2_horizontal():
    """Nitrógeno, geometría y U del caso 2/3, 100 h."""
    c = CASO3
    tk = crear_horizontal("Nitrogen", 101325.0, c["V"], c["d"], c["L"], c["LF"], c["T_air"],
                          c["eta_w"], q=c["q"])
    tk.evaporate(100 * 3600.0)
    return tk
