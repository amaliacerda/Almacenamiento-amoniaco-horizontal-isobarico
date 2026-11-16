"""test_sin_calor_al_vapor: con U_V = 0 el vapor queda a T_sat y BOG = B_L.

Sin calor de pared al vapor y con T(h) = T_sat el vapor no tiene fuente: T = T_sat en
todas las celdas y Q_VL = 0. El BOG (= B_L - dm_V/dt) difiere de B_L solo porque el
volumen de vapor crece al bajar el nivel: dm_V/dt = rho_V_sat dV_V/dt = rho_V_sat B_L/rho_L,
es decir BOG = B_L (1 - rho_V_sat/rho_L) (0.13 % en NH3, 0.57 % en N2). Se comprueba
BOG ~ B_L con 1 % (criterio de la tarea) y la fórmula exacta con 1e-3; B_L inicial
contra el valor analítico U_L A_wet (T_air - T_sat) / dH_LV. Fluidos: amoníaco (caso 3) y nitrógeno.
"""
import numpy as np
import pytest

from conftest import CASO3, crear_horizontal


@pytest.fixture(scope="module", params=[("Ammonia", 101325.0), ("Nitrogen", 101325.0)],
                ids=["amoniaco", "nitrogeno"])
def tk(request):
    fluido, P = request.param
    c = CASO3
    t = crear_horizontal(fluido, P, c["V"], c["d"], c["L"], c["LF"], c["T_air"], c["eta_w"],
                         q=c["q"], U_V=0.0)
    t.evaporate(50 * 3600.0)
    return t


def test_vapor_a_Tsat(tk):
    T_sat = tk.cryogen.T_sat
    assert np.max(np.abs(tk.sol.y[1:] - T_sat)) < 1e-4  # K
    assert np.max(np.abs(tk.data["Tv_avg"] - T_sat)) < 1e-4
    assert np.max(np.abs(tk.data["T_BOG"] - T_sat)) < 1e-4
    assert np.max(np.abs(tk.data["Q_VL"])) < 1e-6  # W


def test_BOG_igual_a_B_L(tk):
    d = tk.data
    assert np.max(np.abs(d["BOG"] - d["B_L"]) / d["B_L"]) < 0.01
    # BOG = B_L (1 - rho_V_sat / rho_L): el vapor solo ocupa el volumen que libera el líquido
    esperado = d["B_L"] * (1 - tk.rho_V_sat / tk.cryogen.rho_L)
    assert np.max(np.abs(d["BOG"] - esperado) / esperado) < 1e-3
    # valor esperado concreto de B_L inicial: U_L A_wet(h0) (T_air - T_sat) / dH_LV
    A_wet0 = tk.geo.A_wet(tk.geo.h_desde_V(tk.LF * tk.V))
    B_L0 = tk.U_L * A_wet0 * (tk.T_air - tk.cryogen.T_sat) / tk.dH_LV
    assert d["B_L"][0] == pytest.approx(B_L0, rel=1e-4)
    assert d["BOG"][0] == pytest.approx(B_L0, rel=0.01)
