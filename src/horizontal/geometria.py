"""Geometría de estanques cilíndricos (horizontal y vertical) con interfaz común.

Referencia: docs/ecuaciones.md, sección 1. Unidades SI (m, m^2, m^3).
z se mide desde el fondo (z = 0) hasta el techo (z = H). h es el nivel del líquido.

Interfaz común (ambas clases, todas vectorizadas en z o h salvo las integrales
escalares y h_desde_V):
    H                   altura interior total [m] (2R en horizontal, V/A_T en vertical)
    A(z)                área de la rebanada horizontal a la altura z [m^2]
    P(z)                área de pared por unidad de altura, dA_pared/dz [m]
    V_L(h)              volumen bajo el nivel h [m^3]
    h_desde_V(V)        nivel h para un volumen de líquido V [m] (inversa de V_L)
    A_wet(h)            área de pared mojada por el líquido [m^2]
    A_dry(h)            área de pared en contacto con el vapor [m^2]
    area_total()        área de pared total [m^2] = A_wet(H)
    volumen_total()     volumen interior total [m^3] = V_L(H)
    integral_A(z0, z1)  volumen exacto entre z0 y z1: int_z0^z1 A(z) dz [m^3]
    integral_P(z0, z1)  pared exacta entre z0 y z1: int_z0^z1 P(z) dz [m^2]

Convención de pared (igual que CryoEvap): las áreas de pared que reciben calor del
ambiente usan el diámetro exterior d_o. En vertical P = pi d_o (docs §6). En
horizontal el manto usa R_o = d_o/2 a igual ángulo: dA_manto/dz = 2 L R_o / w(z);
las tapas planas usan el radio interior (con d_o = d_i se recupera docs §1 exacto).
"""
import numpy as np
from scipy.optimize import brentq


class Geometria:
    """Clase base: define h_desde_V y A_dry a partir de los métodos de cada geometría."""

    H = None
    orientacion = None

    def h_desde_V(self, V):
        """Nivel h [m] tal que V_L(h) = V [m^3]. Resuelve con brentq en [0, H].

        Lanza ValueError si V está fuera de [0, volumen_total()].
        """
        Vt = self.volumen_total()
        if V < 0 or V > Vt:
            raise ValueError("Volumen de líquido %.6g m^3 fuera de [0, %.6g] m^3" % (V, Vt))
        if V == 0:
            return 0.0
        if V == Vt:
            return float(self.H)
        return brentq(lambda h: self.V_L(h) - V, 0.0, self.H, xtol=1e-14 * self.H, rtol=1e-15)

    def A_dry(self, h):
        """Área de pared seca (en contacto con el vapor) [m^2] = area_total - A_wet(h)."""
        return self.area_total() - self.A_wet(h)

    def area_total(self):
        """Área de pared total [m^2]."""
        return float(self.A_wet(self.H))

    def volumen_total(self):
        """Volumen interior total [m^3]."""
        return float(self.V_L(self.H))


class GeometriaVertical(Geometria):
    """Cilindro vertical de CryoEvap (rama main).

    A(z) = pi d_i^2 / 4 constante, P(z) = pi d_o constante, H = V / A_T.
    La pared considerada es solo el manto lateral (el fondo y el techo se tratan
    aparte con Q_b y Q_roof, como en CryoEvap): A_wet(h) = pi d_o h.

    Parámetros: d_i, d_o [m] diámetros interior y exterior; V [m^3] volumen total.
    """

    orientacion = "vertical"

    def __init__(self, d_i, d_o, V):
        self.d_i = float(d_i)
        self.d_o = float(d_o)
        self.A_T = np.pi * self.d_i ** 2 / 4  # [m^2] sección transversal
        self.H = float(V) / self.A_T          # [m] altura interior

    def A(self, z):
        """A(z) = pi d_i^2/4 [m^2]."""
        return self.A_T * np.ones_like(np.asarray(z, dtype=float))

    def P(self, z):
        """P(z) = pi d_o [m] (manto por unidad de altura)."""
        return np.pi * self.d_o * np.ones_like(np.asarray(z, dtype=float))

    def V_L(self, h):
        """V_L(h) = A_T h [m^3]."""
        return self.A_T * np.asarray(h, dtype=float)

    def h_desde_V(self, V):
        """h = V / A_T [m] (inversa exacta, sin brentq)."""
        Vt = self.volumen_total()
        if V < 0 or V > Vt * (1 + 1e-12):
            raise ValueError("Volumen de líquido %.6g m^3 fuera de [0, %.6g] m^3" % (V, Vt))
        return float(V) / self.A_T

    def A_wet(self, h):
        """A_wet(h) = pi d_o h [m^2] (solo manto)."""
        return np.pi * self.d_o * np.asarray(h, dtype=float)

    def integral_A(self, z0, z1):
        """int_z0^z1 A dz = A_T (z1 - z0) [m^3]."""
        return self.A_T * (np.asarray(z1, dtype=float) - np.asarray(z0, dtype=float))

    def integral_P(self, z0, z1):
        """int_z0^z1 P dz = pi d_o (z1 - z0) [m^2]."""
        return np.pi * self.d_o * (np.asarray(z1, dtype=float) - np.asarray(z0, dtype=float))


class GeometriaHorizontal(Geometria):
    """Cilindro horizontal de radio R = d_i/2 y largo L (parte cilíndrica).

    Ecuaciones (docs/ecuaciones.md §1), con w(z) = sqrt(z (2R - z)),
    phi(z) = acos((R - z)/R) y S(z) = R^2 phi(z) - (R - z) w(z) (área del segmento):
        A(z)     = 2 L w(z)                       [+ pi w^2 con tapas hemisféricas]
        P(z)     = 2 L R_o / w(z) + 4 w(z)        [tapas planas]
                 = 2 L R_o / w(z) + 2 pi R        [tapas hemisféricas]
        V_L(h)   = L S(h)                         [+ pi h^2 (3R - h)/3 hemisféricas]
        A_wet(h) = 2 L R_o phi(h) + 2 S(h)        [planas]
                 = 2 L R_o phi(h) + 2 pi R h      [hemisféricas]
        int A dz = L [S(z1) - S(z0)]              [+ pi (R z^2 - z^3/3) evaluado]
        int P dz = 2 L R_o [phi(z1) - phi(z0)] + 2 [S(z1) - S(z0)]   [planas]
    P es singular en z = 0 y z = 2R (devuelve inf allí), pero su integral es finita.

    Parámetros: d_i, d_o [m]; L [m] largo de la parte cilíndrica;
    tapas: 'planas' (por defecto) o 'hemisfericas'.
    """

    orientacion = "horizontal"

    def __init__(self, d_i, d_o, L, tapas="planas"):
        if tapas not in ("planas", "hemisfericas"):
            raise ValueError("tapas debe ser 'planas' o 'hemisfericas', no %r" % (tapas,))
        self.d_i = float(d_i)
        self.d_o = float(d_o)
        self.R = self.d_i / 2           # [m] radio interior
        self.R_o = self.d_o / 2         # [m] radio exterior (manto)
        self.L = float(L)               # [m] largo
        self.tapas = tapas
        self.H = 2 * self.R             # [m] altura interior

    # Funciones auxiliares (z se recorta a [0, 2R] solo para evitar NaN por redondeo)
    def _zc(self, z):
        return np.clip(np.asarray(z, dtype=float), 0.0, self.H)

    def w(self, z):
        """Media cuerda w(z) = sqrt(z (2R - z)) [m]."""
        z = self._zc(z)
        return np.sqrt(np.clip(z * (self.H - z), 0.0, None))

    def phi(self, z):
        """Ángulo phi(z) = acos((R - z)/R) [rad], en [0, pi]."""
        z = self._zc(z)
        return np.arccos(np.clip((self.R - z) / self.R, -1.0, 1.0))

    def S(self, z):
        """Área del segmento circular bajo z: S = R^2 phi - (R - z) w [m^2]."""
        z = self._zc(z)
        return self.R ** 2 * self.phi(z) - (self.R - z) * self.w(z)

    # Interfaz común
    def A(self, z):
        """A(z) = 2 L w(z) (+ pi w^2 si hemisféricas) [m^2]."""
        w = self.w(z)
        a = 2 * self.L * w
        if self.tapas == "hemisfericas":
            a = a + np.pi * w ** 2
        return a

    def P(self, z):
        """P(z) = 2 L R_o / w(z) + 4 w(z) (planas) o + 2 pi R (hemisféricas) [m]."""
        w = self.w(z)
        with np.errstate(divide="ignore"):
            manto = 2 * self.L * self.R_o / w
        if self.tapas == "planas":
            return manto + 4 * w
        return manto + 2 * np.pi * self.R * np.ones_like(w)

    def V_L(self, h):
        """V_L(h) = L [R^2 theta - (R - h) sqrt(2Rh - h^2)] (+ casquete esférico) [m^3]."""
        v = self.L * self.S(h)
        if self.tapas == "hemisfericas":
            hc = self._zc(h)
            v = v + np.pi * hc ** 2 * (3 * self.R - hc) / 3
        return v

    def A_wet(self, h):
        """A_wet(h) = 2 L R_o theta + 2 S(h) (planas) o + 2 pi R h (hemisféricas) [m^2]."""
        a = 2 * self.L * self.R_o * self.phi(h)
        if self.tapas == "planas":
            return a + 2 * self.S(h)
        return a + 2 * np.pi * self.R * self._zc(h)

    def _prim_A(self, z):
        """Primitiva de A: L S(z) (+ pi (R z^2 - z^3/3)) [m^3]."""
        p = self.L * self.S(z)
        if self.tapas == "hemisfericas":
            zc = self._zc(z)
            p = p + np.pi * (self.R * zc ** 2 - zc ** 3 / 3)
        return p

    def _prim_P(self, z):
        """Primitiva de P: 2 L R_o phi(z) + 2 S(z) (planas) o + 2 pi R z [m^2]."""
        p = 2 * self.L * self.R_o * self.phi(z)
        if self.tapas == "planas":
            return p + 2 * self.S(z)
        return p + 2 * np.pi * self.R * self._zc(z)

    def integral_A(self, z0, z1):
        """int_z0^z1 A(z) dz exacta [m^3] (volumen de la celda)."""
        return self._prim_A(z1) - self._prim_A(z0)

    def integral_P(self, z0, z1):
        """int_z0^z1 P(z) dz exacta [m^2] (pared de la celda)."""
        return self._prim_P(z1) - self._prim_P(z0)
