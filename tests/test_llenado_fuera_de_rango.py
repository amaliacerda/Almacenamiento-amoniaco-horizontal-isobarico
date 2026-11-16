"""test_llenado_fuera_de_rango: LF fuera del intervalo abierto (0, 1) falla con mensaje claro."""
import pytest
from cryoevap.cryogens import Cryogen

from horizontal import HorizontalTank

D, V, L = 1.676, 4.89, 2.217


@pytest.mark.parametrize("LF", [0.0, 1.0, -0.1, 1.2])
def test_LF_invalido_en_constructor(LF):
    with pytest.raises(ValueError, match="LF.*fuera del intervalo abierto"):
        HorizontalTank(D, D, V, LF=LF, L=L)


@pytest.mark.parametrize("LF", [0.0, 1.0, -0.1, 1.2])
def test_LF_invalido_al_simular(LF):
    """Si LF se asigna mal tras crear el tanque, la simulación no corre."""
    tk = HorizontalTank(D, D, V, LF=0.3, L=L)
    cr = Cryogen(name="Ammonia")
    cr.set_coolprops(101325.0)
    tk.cryogen = cr
    tk.set_HeatTransProps(U_L=0.01, U_V=0.01, T_air=293.15, Q_b_fixed=0, Q_roof=0, eta_w=0.97)
    with pytest.raises(ValueError, match="LF.*fuera del intervalo abierto"):
        tk.LF = LF
        tk.evaporate(3600.0)


def test_LF_vertical_invalido():
    with pytest.raises(ValueError, match="LF"):
        HorizontalTank(0.201, 0.204, 6.75e-3, LF=1.5, orientacion="vertical")


def test_LF_valido_no_falla():
    tk = HorizontalTank(D, D, V, LF=0.999, L=L)
    assert tk.LF == pytest.approx(0.999)
