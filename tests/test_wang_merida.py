"""test_wang_merida: caso 2 (hidrógeno) contra Wang y Mérida 2024, figura 12.

Datos digitalizados en data/Data_wang_merida.xlsx, hoja 12e_BOG, columnas
t_horiz_h y BOG_horiz_kg_h (kg/h). Se compara el BOG del modelo (interpolado en los
tiempos digitalizados) con criterio AGENTS.md: diferencia relativa <= 10 %.
Se evalúa (a) el régimen estacionario, t >= 15 h, y (b) toda la curva (promedio de
la diferencia relativa en los 15 puntos). Sin tolerancias relajadas.
"""
import numpy as np
import openpyxl
import pytest

from conftest import RUTA_WANG

TOL = 0.10


@pytest.fixture(scope="module")
def wang():
    ws = openpyxl.load_workbook(RUTA_WANG, data_only=True)["12e_BOG"]
    filas = [r for r in ws.iter_rows(min_row=2, values_only=True) if r[2] is not None]
    return np.array([r[2] for r in filas]), np.array([r[3] for r in filas])


def _bog_kg_h(tk, t_h):
    return np.interp(t_h, tk.data["Time"] / 3600, tk.data["BOG"] * 3600)


def test_wang_merida_estacionario(tanque_h2_caso2, wang):
    t, bog = wang
    m = t >= 15
    modelo = _bog_kg_h(tanque_h2_caso2, t[m])
    dif = np.abs(modelo - bog[m]) / bog[m]
    assert dif.mean() <= TOL, "dif media BOG (t >= 15 h) = %.1f %% (modelo %.4f, Wang %.4f kg/h)" % (
        100 * dif.mean(), modelo.mean(), bog[m].mean())
    assert dif.max() <= TOL, "dif máxima BOG (t >= 15 h) = %.1f %%" % (100 * dif.max())


def test_wang_merida_curva_completa(tanque_h2_caso2, wang):
    t, bog = wang
    modelo = _bog_kg_h(tanque_h2_caso2, t)
    dif = np.abs(modelo - bog) / bog
    assert dif.mean() <= TOL, "dif media BOG (curva completa) = %.1f %%" % (100 * dif.mean())
