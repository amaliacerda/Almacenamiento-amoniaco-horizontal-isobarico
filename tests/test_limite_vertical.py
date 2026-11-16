"""test_limite_vertical: el modelo nuevo en modo vertical reproduce CryoEvap original.

Caso 1 de AGENTS.md (notebook LN2_small): nitrógeno, P = 1e5 Pa, d_i = 0.201 m,
d_o = 0.204 m, V = 6.75e-3 m^3, LF = 0.278, U_L = U_V = 0.026, T_air = 298.15 K,
eta_w = 0.963, dz = 0.005 m, t = 2 tau, registro cada 60 s. Ambos modelos usan la
misma malla (caras linspace(0, 1, n_z)). La corrida de CryoEvap (RK45) tarda ~30 s
y se hace una sola vez (fixture de módulo).

Tolerancia: 1 % (AGENTS.md) en BOG y Tv_avg. El BOG en t = 0 de CryoEvap sale de una
reconstrucción por diferencias finitas con salto inicial artificial (3.82e-6 frente a
4.08e-6 kg/s del valor analítico); por eso la serie temporal se compara desde el
primer registro (t = 60 s), y el valor final se compara por separado.
"""
import numpy as np
import pytest
from cryoevap.cryogens import Cryogen
from cryoevap.storage_tanks import Tank

from horizontal import HorizontalTank

TOL = 0.01
PAR = dict(d_i=0.201, d_o=0.204, V=6.75e-3, LF=0.278, U=0.026, T_air=298.15, eta_w=0.963, P=100000.0)


def _preparar(tk):
    c = Cryogen(name="nitrogen")
    c.set_coolprops(PAR["P"])
    tk.cryogen = c
    tk.set_HeatTransProps(U_L=PAR["U"], U_V=PAR["U"], T_air=PAR["T_air"], Q_b_fixed=0, Q_roof=0,
                          eta_w=PAR["eta_w"]) if isinstance(tk, HorizontalTank) else \
        tk.set_HeatTransProps(PAR["U"], PAR["U"], PAR["T_air"], Q_b_fixed=0, Q_roof=0, eta_w=PAR["eta_w"])
    tau_min = np.floor(tk.tau / 60) + 1
    n_z = 1 + int(np.round(tk.l_V / 0.005))
    tk.z_grid = np.linspace(0, 1, n_z)
    tk.U_roof = 0
    tk.time_interval = 60
    return 60 * tau_min * 2


@pytest.fixture(scope="module")
def corridas():
    nuevo = HorizontalTank(PAR["d_i"], PAR["d_o"], PAR["V"], LF=PAR["LF"], orientacion="vertical")
    t_f = _preparar(nuevo)
    nuevo.evaporate(t_f)
    orig = Tank(PAR["d_i"], PAR["d_o"], PAR["V"], PAR["LF"])
    t_f2 = _preparar(orig)
    assert t_f == t_f2
    orig.evaporate(t_f2)
    return nuevo, orig


def test_limite_vertical_BOG(corridas):
    nuevo, orig = corridas
    a, b = np.asarray(nuevo.data["BOG"]), np.asarray(orig.data["BOG"])
    assert a.shape == b.shape
    assert a[-1] == pytest.approx(b[-1], rel=TOL), "BOG final"
    dif = np.abs(a[1:] - b[1:]) / np.abs(b[1:])
    assert dif.max() < TOL, "max dif relativa BOG (t >= 60 s) = %.3e" % dif.max()


def test_limite_vertical_Tv_avg(corridas):
    nuevo, orig = corridas
    a, b = np.asarray(nuevo.data["Tv_avg"]), np.asarray(orig.data["Tv_avg"])
    assert a.shape == b.shape
    assert a[-1] == pytest.approx(b[-1], rel=TOL), "Tv_avg final"
    dif = np.abs(a - b) / b
    assert dif.max() < TOL, "max dif relativa Tv_avg = %.3e" % dif.max()


def test_limite_vertical_otras_magnitudes(corridas):
    """Nivel y calor líquido-vapor también coinciden (misma física)."""
    nuevo, orig = corridas
    assert nuevo.data["V_L"][-1] == pytest.approx(orig.data["V_L"][-1], rel=TOL)
    assert nuevo.data["B_L"][-1] == pytest.approx(orig.data["B_L"][-1], rel=TOL)
