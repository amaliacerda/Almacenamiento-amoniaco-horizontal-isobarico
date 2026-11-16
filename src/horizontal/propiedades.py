"""Propiedades del vapor sobrecalentado a presión fija, tabuladas desde CoolProp.

El modelo es isobárico, así que rho, h, cp y k del vapor dependen solo de T.
Se tabulan una vez con CoolProp (AbstractState HEOS con fase gas impuesta, que
evalúa sin ambigüedad en T = T_sat, sin el desplazamiento de 1e-3 K que usa
CryoEvap) y se interpolan con splines cúbicos. cp se toma como la derivada del
spline de h, de modo que cp = dh/dT exactamente: esto es lo que hace que el
balance de energía semidiscreto cierre (docs §5).
"""
import numpy as np
import CoolProp.CoolProp as CP
from scipy.interpolate import CubicSpline


class TablaVapor:
    """Tabla de propiedades del vapor a presión P [Pa] entre T_min y T_max [K].

    Métodos (vectorizados en T [K]):
        rho(T)   densidad [kg/m^3]
        drho(T)  d rho / dT [kg/m^3/K]
        h(T)     entalpía específica [J/kg] (misma referencia que CoolProp)
        cp(T)    dh/dT [J/kg/K]
        k(T)     conductividad térmica [W/m/K]
    """

    def __init__(self, fluido, P, T_sat, T_max, n=500):
        self.fluido = fluido
        self.P = float(P)
        self.T_sat = float(T_sat)
        AS = CP.AbstractState("HEOS", fluido)
        AS.specify_phase(CP.iphase_gas)
        # Malla en T más densa cerca de T_sat (donde las propiedades cambian más)
        s = np.linspace(0.0, 1.0, n)
        T = self.T_sat + (float(T_max) - self.T_sat) * s ** 2
        # Algunos puntos levemente subenfriados (vapor metaestable) para que el
        # spline no extrapole si T baja un poco de T_sat por error numérico.
        sub = self.T_sat - np.array([0.2, 0.1, 0.05])
        Tlist, rho, h, k = [], [], [], []
        for Ti in np.concatenate([np.sort(sub), T]):
            try:
                AS.update(CP.PT_INPUTS, self.P, Ti)
                r, hh, kk = AS.rhomass(), AS.hmass(), AS.conductivity()
            except ValueError:
                continue
            if not (np.isfinite(r) and np.isfinite(hh) and np.isfinite(kk)):
                continue
            Tlist.append(Ti); rho.append(r); h.append(hh); k.append(kk)
        Tlist = np.array(Tlist)
        self.T_tabla = Tlist
        self._rho = CubicSpline(Tlist, rho)
        self._drho = self._rho.derivative()
        self._h = CubicSpline(Tlist, h)
        self._cp = self._h.derivative()
        self._k = CubicSpline(Tlist, k)

    def rho(self, T):
        return self._rho(T)

    def drho(self, T):
        return self._drho(T)

    def h(self, T):
        return self._h(T)

    def cp(self, T):
        return self._cp(T)

    def k(self, T):
        return self._k(T)
