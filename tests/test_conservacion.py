"""test_conservacion_masa y test_conservacion_energia (criterio AGENTS.md: error < 1 %).

Cada test hace dos comprobaciones: (a) el error que informa el modelo
(data['error_masa'], data['error_energia']) y (b) una verificación independiente
armada solo con las salidas guardadas (data y sol), con integración trapezoidal.
Casos: amoníaco (caso 3), nitrógeno (misma geometría) y el caso 2 (hidrógeno, que
reutiliza el fixture de la validación).
"""
import numpy as np
import pytest

TOL = 0.01


@pytest.fixture(params=["tanque_nh3_caso3", "tanque_n2_horizontal", "tanque_h2_caso2"],
                ids=["amoniaco_caso3", "nitrogeno", "hidrogeno_caso2"])
def tanque(request):
    return request.getfixturevalue(request.param)


def test_conservacion_masa(tanque):
    d = tanque.data
    t = d["Time"]
    m0 = d["m_L"][0] + d["m_V"][0]
    mf = d["m_L"][-1] + d["m_V"][-1]
    bog_int = np.trapezoid(d["BOG"], t)
    assert bog_int > 0
    err_indep = abs((m0 - mf) - bog_int) / bog_int
    assert err_indep < TOL, "error masa independiente = %.3e" % err_indep
    assert d["error_masa"] < TOL, "error_masa del modelo = %.3e" % d["error_masa"]


def test_conservacion_energia(tanque):
    """Delta H del estanque + int B h_V(T_techo) dt = int Q_in dt."""
    d = tanque.data
    t = d["Time"]
    sol = tanque.sol
    # Entalpía total H(t) = m_L h_L + sum rho h V_j en cada tiempo guardado
    H = []
    h_techo = []
    for i in range(len(t)):
        y = sol.y[:, i]
        _, dd = tanque._modelo(y)
        H.append(tanque.cryogen.rho_L * y[0] * tanque.cryogen.h_L + np.sum(dd["rho"] * dd["hT"] * dd["V_j"]))
        h_techo.append(float(tanque.tabla.h(y[-1])))
    H = np.array(H)
    sale = np.trapezoid(d["BOG"] * np.array(h_techo), t)
    entra = np.trapezoid(d["Q_in"], t)
    assert entra > 0
    err_indep = abs((H[-1] - H[0]) + sale - entra) / entra
    assert err_indep < TOL, "error energía independiente = %.3e" % err_indep
    assert d["error_energia"] < TOL, "error_energia del modelo = %.3e" % d["error_energia"]
