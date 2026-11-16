"""test_geometria: propiedades exactas de GeometriaHorizontal / GeometriaVertical."""
import numpy as np
import pytest
from scipy.integrate import quad

from horizontal import GeometriaHorizontal, GeometriaVertical

TOL = 1e-8
R = 0.838
L = 2.217
H_NIVELES = [0.05, 0.3, 0.838, 1.2, 1.6, 1.67]


def _geos():
    return [
        ("planas_do=di", GeometriaHorizontal(2 * R, 2 * R, L, "planas")),
        ("planas_do>di", GeometriaHorizontal(2 * R, 2 * R + 0.02, L, "planas")),
        ("hemisfericas", GeometriaHorizontal(2 * R, 2 * R + 0.02, L, "hemisfericas")),
    ]


@pytest.mark.parametrize("nombre,geo", _geos())
def test_volumen_total(nombre, geo):
    """V_L(2R) = pi R^2 L (planas) o + 4/3 pi R^3 (hemisféricas)."""
    esperado = np.pi * R ** 2 * L
    if geo.tapas == "hemisfericas":
        esperado += 4 / 3 * np.pi * R ** 3
    assert geo.V_L(2 * R) == pytest.approx(esperado, rel=1e-12)
    assert geo.volumen_total() == pytest.approx(esperado, rel=1e-12)


@pytest.mark.parametrize("nombre,geo", _geos())
def test_area_mojada_total(nombre, geo):
    """A_wet(2R) = área total de pared (manto con R_o + tapas)."""
    Ro = geo.R_o
    esperado = 2 * np.pi * Ro * L + (2 * np.pi * R ** 2 if geo.tapas == "planas" else 4 * np.pi * R ** 2)
    assert geo.A_wet(2 * R) == pytest.approx(esperado, rel=1e-12)
    assert geo.area_total() == pytest.approx(esperado, rel=1e-12)
    assert geo.A_dry(2 * R) == pytest.approx(0.0, abs=1e-12)
    assert geo.A_wet(0.0) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("nombre,geo", _geos())
@pytest.mark.parametrize("h", H_NIVELES)
def test_V_y_A_wet_vs_integral_numerica(nombre, geo, h):
    """V_L(h) = int A dz y A_wet(h) = int P dz (cuadratura numérica), rel 1e-8.

    P es singular (1/sqrt) en z = 0: se integra en el ángulo, z = R (1 - cos t),
    con lo que dz = R sin t dt y el integrando es suave.
    """
    t_h = np.arccos((R - h) / R)

    def en_angulo(f):
        return quad(lambda t: f(R * (1 - np.cos(t))) * R * np.sin(t), 0, t_h,
                    epsabs=0, epsrel=1e-12, limit=200)[0]

    V_num = en_angulo(lambda z: float(geo.A(z)))
    A_num = en_angulo(lambda z: float(geo.P(z)))
    assert geo.V_L(h) == pytest.approx(V_num, rel=TOL)
    assert geo.A_wet(h) == pytest.approx(A_num, rel=TOL)
    # integrales por celda [z0, z1]
    z0, z1 = 0.4 * h, h
    t0 = np.arccos((R - z0) / R)
    celda_A = quad(lambda t: float(geo.A(R * (1 - np.cos(t)))) * R * np.sin(t), t0, t_h,
                   epsabs=0, epsrel=1e-12)[0]
    celda_P = quad(lambda t: float(geo.P(R * (1 - np.cos(t)))) * R * np.sin(t), t0, t_h,
                   epsabs=0, epsrel=1e-12)[0]
    assert geo.integral_A(z0, z1) == pytest.approx(celda_A, rel=TOL)
    assert geo.integral_P(z0, z1) == pytest.approx(celda_P, rel=TOL)


@pytest.mark.parametrize("nombre,geo", _geos())
@pytest.mark.parametrize("h", [1e-3, 0.05, 0.4, 0.838, 1.3, 1.675])
def test_h_desde_V_ida_y_vuelta(nombre, geo, h):
    """h_desde_V(V_L(h)) = h."""
    assert geo.h_desde_V(float(geo.V_L(h))) == pytest.approx(h, rel=1e-9)


def test_h_desde_V_fuera_de_rango():
    geo = GeometriaHorizontal(2 * R, 2 * R, L)
    with pytest.raises(ValueError):
        geo.h_desde_V(-1e-3)
    with pytest.raises(ValueError):
        geo.h_desde_V(geo.volumen_total() * 1.01)


def test_media_altura_mitad_volumen():
    """Por simetría, V_L(R) = V_total / 2 y A_wet(R) = A_total / 2."""
    geo = GeometriaHorizontal(2 * R, 2 * R, L)
    assert geo.V_L(R) == pytest.approx(geo.volumen_total() / 2, rel=1e-12)
    assert geo.A_wet(R) == pytest.approx(geo.area_total() / 2, rel=1e-12)


def test_geometria_vertical():
    """Vertical: V_L = pi d_i^2/4 h, A_wet = pi d_o h, integrales lineales."""
    d_i, d_o, V = 0.201, 0.204, 6.75e-3
    g = GeometriaVertical(d_i, d_o, V)
    assert g.H == pytest.approx(V / (np.pi * d_i ** 2 / 4), rel=1e-12)
    assert g.volumen_total() == pytest.approx(V, rel=1e-12)
    assert g.A_wet(0.05) == pytest.approx(np.pi * d_o * 0.05, rel=1e-12)
    assert g.integral_P(0.02, 0.07) == pytest.approx(np.pi * d_o * 0.05, rel=1e-12)
    assert g.h_desde_V(g.V_L(0.1)) == pytest.approx(0.1, rel=1e-12)
