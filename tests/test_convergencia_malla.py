"""test_convergencia_malla: el BOG cambia < 0.5 % entre las dos mallas más finas.

Mallas de 20, 40, 80 y 160 celdas de vapor. Se comparan las dos más finas (80 y 160)
en dos instantes: transitorio (t = 2 h) y final (t = 40 h), para amoníaco (caso 3,
fluido de aplicación) y nitrógeno. Además se exige que el error disminuya al refinar
(|BOG_80 - BOG_160| <= |BOG_40 - BOG_160|).
"""
import numpy as np
import pytest

from conftest import CASO3, crear_horizontal

TOL = 0.005
MALLAS = [20, 40, 80, 160]
T_F = 40 * 3600.0


@pytest.fixture(scope="module", params=["Ammonia", "Nitrogen"], ids=["amoniaco", "nitrogeno"])
def bogs(request):
    c = CASO3
    out = {}
    for n in MALLAS:
        tk = crear_horizontal(request.param, 101325.0, c["V"], c["d"], c["L"], c["LF"], c["T_air"],
                              c["eta_w"], q=c["q"], n_celdas=n, time_interval=3600.0)
        tk.evaporate(T_F)
        out[n] = tk.data["BOG"]
    return out


@pytest.mark.parametrize("idx,nombre", [(2, "t = 2 h"), (-1, "t = 40 h")])
def test_convergencia_malla(bogs, idx, nombre):
    fino, previo, medio = bogs[160][idx], bogs[80][idx], bogs[40][idx]
    dif = abs(previo - fino) / fino
    assert dif < TOL, "BOG(80) vs BOG(160) en %s: %.3e" % (nombre, dif)
    assert abs(previo - fino) <= abs(medio - fino), "no monótona al refinar"
