"""Estanque cilíndrico (horizontal o vertical) con evaporación isobárica.

Modelo 1D de docs/ecuaciones.md, escrito sobre CryoEvap (rama main, 859e9ba,
GPL-3.0, https://github.com/felipehuerta17/CryoEvap). Se hereda de
cryoevap.storage_tanks.Tank (API: evaporate, data, plot_*, tau, b_l_dot) y se
reutiliza cryoevap.cryogens.Cryogen (propiedades de saturación). Se reescriben
la geometría y la ecuación del vapor. Unidades SI.

Esquema numérico (vapor, z entre h y H):
  - Volúmenes finitos sobre una grilla móvil en la coordenada adimensional
    zeta = (z - h)/(H - h). Caras zeta_f fijas; caras físicas z_f = h + zeta_f (H - h).
  - Volumen V_j = int A dz y pared P_j = int P dz exactos por celda (geometria.py).
  - Conducción k_f A(z_f) (T_j+1 - T_j)/dz en caras interiores; en la interfase
    (T = T_sat) gradiente de segundo orden con las dos primeras celdas; en el
    techo flujo nulo (A(2R) = 0 en horizontal; techo aislado en vertical).
  - Advección escrita con el flujo de masa G_f (kg/s) que cruza cada cara móvil,
    upwind (G_f >= 0). Forma en T obtenida restando h * continuidad a la energía:
        rho_j cp_j V_j dT_j/dt = G_j-1/2 [h(T_j-1) - h(T_j)] + F_j+1/2 - F_j-1/2
                                  + (1 - eta_w) U_V P_j (T_air - T_j)
  - Integración con solve_ivp (BDF por defecto, o Radau) con jacobiano denso
    por diferencias finitas: el sistema está acoplado globalmente (ver evaporate).
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import OptimizeResult

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

from .geometria import GeometriaHorizontal, GeometriaVertical
from .propiedades import TablaVapor
from . import balances


class HorizontalTank(Tank):
    """Estanque cilíndrico horizontal (o vertical, para regresión contra CryoEvap).

    Constructor:
        HorizontalTank(d_i, d_o, V, LF=0.97, L=None, orientacion='horizontal',
                       tapas='planas')
        d_i, d_o [m]  diámetros interior y exterior
        V [m^3]       volumen total
        LF [-]        llenado inicial, en el intervalo abierto (0, 1); si no, ValueError
        L [m]         largo de la parte cilíndrica (solo horizontal). Si es None,
                      se calcula de V y d_i. Si se da, debe ser coherente con V (1e-3).
        orientacion   'horizontal' o 'vertical' (modo vertical = geometría de CryoEvap)
        tapas         'planas' (por defecto) o 'hemisfericas' (solo horizontal)

    Diagnóstico: tank.diagnostico(t) devuelve las cantidades por celda (z_f, V_j, P_j,
    T_j, G) y los calores Q_L, Q_wi, Q_VL, Q_in por separado (ver su docstring).
    data['G_min'] es el mínimo de G_f observado en la corrida (warning si < 0).

    Opciones numéricas (atributos, modificables antes de evaporate):
        n_celdas (40)            número de celdas de vapor. Alternativa compatible con
                                 CryoEvap: asignar tank.z_grid = np.linspace(0, 1, n)
                                 define las n caras (n - 1 celdas).
        flujo_masa               'conservativo' (por defecto): G_f desde la continuidad
                                 local; 'uniforme': mdot = B_L uniforme (docs §3).
        propiedades              'locales' (por defecto) o 'promedio' (ponderadas por A).
        metodo ('BDF'), rtol (1e-8), atol_T (1e-6 K), atol_V (1e-12 * V)
        time_interval (3600 s)   paso de los datos guardados en data
        cuadratura_gauss (3)     puntos de Gauss por paso del integrador para los balances

    Resultados: tank.sol (t, y en los tiempos guardados; y[0] = V_L, y[1:] = T de las
    celdas), tank.z_grid (centros de celda en zeta), tank.data (ver _reconstruct).
    """

    def __init__(self, d_i, d_o, V, LF=0.97, L=None, orientacion="horizontal", tapas="planas"):
        super().__init__(d_i, d_o, V, LF=0.5)  # LF se valida abajo con el setter
        self._zeta_caras = None  # Tank.__init__ fija z_grid de 100 nodos; se usa n_celdas
        if orientacion == "vertical":
            self.geo = GeometriaVertical(d_i, d_o, V)
        elif orientacion == "horizontal":
            R = d_i / 2
            if tapas == "planas":
                L_V = V / (np.pi * R ** 2)
            else:
                L_V = (V - 4 / 3 * np.pi * R ** 3) / (np.pi * R ** 2)
            if L is None:
                L = L_V
            elif abs(L - L_V) > 1e-3 * L_V:
                raise ValueError("L = %.6g m no es coherente con V = %.6g m^3 y d_i = %.6g m "
                                 "(L esperado = %.6g m)" % (L, V, d_i, L_V))
            self.geo = GeometriaHorizontal(d_i, d_o, L, tapas=tapas)
        else:
            raise ValueError("orientacion debe ser 'horizontal' o 'vertical', no %r" % (orientacion,))
        self.orientacion = orientacion
        self.V = float(self.geo.volumen_total())
        self.l = float(self.geo.H)  # [m] altura interior (en vertical igual a Tank.l)
        self.LF = LF

        # Opciones numéricas
        self.n_celdas = 40
        self._zeta_caras = None
        self.flujo_masa = "conservativo"
        self.propiedades = "locales"
        self.metodo = "BDF"
        self.rtol = 1e-8
        self.atol_T = 1e-6
        self.atol_V = None
        self.cuadratura_gauss = 3
        self.tabla = None
        self.Q_b_fixed = None
        self.Q_roof = 0.0
        self.eta_w = 0.0

    # ------------------------------------------------------------------ API
    @property
    def LF(self):
        """Llenado (fracción de volumen líquido) [-]. Debe estar en (0, 1)."""
        return self._LF

    @LF.setter
    def LF(self, valor):
        valor = float(valor)
        if not (0.0 < valor < 1.0):
            raise ValueError("LF = %r fuera del intervalo abierto (0, 1): el llenado debe ser "
                             "mayor que 0 y menor que 1." % (valor,))
        self._LF = valor

    @property
    def z_grid(self):
        """Centros de celda en zeta = (z - h)/(H - h) [-] (largo = número de celdas).

        Asignar z_grid = caras (array creciente de 0 a 1, como en CryoEvap) fija la
        malla con len(z_grid) - 1 celdas.
        """
        zf = self.zeta_caras
        return 0.5 * (zf[1:] + zf[:-1])

    @z_grid.setter
    def z_grid(self, caras):
        caras = np.asarray(caras, dtype=float)
        if caras.ndim != 1 or len(caras) < 4 or caras[0] != 0 or caras[-1] != 1 \
                or np.any(np.diff(caras) <= 0):
            raise ValueError("z_grid debe ser un arreglo creciente de caras de 0 a 1 con >= 4 puntos")
        self._zeta_caras = caras

    @property
    def zeta_caras(self):
        """Caras de las celdas en zeta [-], largo n_celdas + 1."""
        if getattr(self, "_zeta_caras", None) is not None:
            return self._zeta_caras
        return np.linspace(0.0, 1.0, getattr(self, "n_celdas", 99) + 1)

    def set_HeatTransProps(self, *, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0, eta_w=0):
        """Propiedades de transferencia de calor. SOLO argumentos con nombre.

        U_L, U_V [W/m^2/K]  coeficientes globales de pared mojada y seca
        T_air [K]           temperatura ambiente
        Q_b_fixed [W]       calor por el fondo. Vertical: si es None, U_L A_T (T_air - T_sat)
                            como CryoEvap. Horizontal: si es None, 0 (no hay fondo aparte).
        Q_roof [W]          calor fijo por el techo, se suma a la celda superior
        eta_w [-]           fracción del calor de pared seca que va directo a la interfase
        """
        super().set_HeatTransProps(U_L, U_V, T_air, Q_b_fixed=Q_b_fixed, Q_roof=Q_roof, eta_w=eta_w)

    @property
    def Q_b(self):
        """Calor por el fondo [W] (ver set_HeatTransProps)."""
        if self.Q_b_fixed is not None:
            return self.Q_b_fixed
        if self.orientacion == "vertical":
            return self.U_L * self.geo.A_T * (self.T_air - self.cryogen.T_sat)
        return 0.0

    # Propiedades heredadas reescritas con la geometría (usadas por tau, v_z, b_l_dot)
    @property
    def h_L(self):
        """Nivel de líquido para el llenado actual LF [m]."""
        return self.geo.h_desde_V(self.LF * self.V)

    @property
    def l_V(self):
        """Altura del espacio de vapor H - h [m]."""
        return self.geo.H - self.h_L

    @property
    def A_T(self):
        """Área de la interfase A(h) [m^2] (en vertical, pi d_i^2/4)."""
        if not hasattr(self, "geo"):
            return np.pi * self.d_i ** 2 / 4
        return float(self.geo.A(self.h_L))

    @A_T.setter
    def A_T(self, valor):
        pass  # Tank.__init__ lo asigna; aquí se calcula desde la geometría

    @property
    def A_L(self):
        """Pared mojada A_wet(h) [m^2]."""
        return float(self.geo.A_wet(self.h_L))

    @property
    def A_V(self):
        """Pared seca A_dry(h) [m^2]."""
        return float(self.geo.A_dry(self.h_L))

    # ------------------------------------------------------------ evaporate
    def evaporate(self, t_f):
        """Simula la evaporación isobárica durante t_f [s].

        Guarda tank.sol (tiempos 0, time_interval, ..., <= t_f), tank.sol_completa
        (salida de solve_ivp con dense_output) y llena tank.data, incluidos
        data['error_masa'] y data['error_energia'] (balances.py).
        """
        self.LF = self.LF  # revalida
        cr = self.cryogen
        if not isinstance(cr, Cryogen) or not cr.rho_L:
            raise ValueError("Asigna tank.cryogen = Cryogen(name=...) con set_coolprops(P) antes de evaporar")
        for nombre in ("U_L", "U_V", "T_air"):
            if not hasattr(self, nombre):
                raise ValueError("Llama set_HeatTransProps(U_L=..., U_V=..., T_air=...) antes de evaporar")
        if self.flujo_masa not in ("conservativo", "uniforme"):
            raise ValueError("flujo_masa debe ser 'conservativo' o 'uniforme'")
        if self.propiedades not in ("locales", "promedio"):
            raise ValueError("propiedades debe ser 'locales' o 'promedio'")

        T_max = max(self.T_air, cr.T_sat) + 20.0
        if self.tabla is None or self.tabla.fluido != cr.name or self.tabla.P != cr.P \
                or self.tabla.T_tabla[-1] < T_max - 1:
            self.tabla = TablaVapor(cr.name, cr.P, cr.T_sat, T_max)
        tb = self.tabla
        self.h_V_sat = float(tb.h(cr.T_sat))       # entalpía del vapor saturado [J/kg]
        self.dH_LV = self.h_V_sat - cr.h_L          # calor latente [J/kg]
        self.rho_V_sat = float(tb.rho(cr.T_sat))
        self.k_V_sat = float(tb.k(cr.T_sat))

        self._zf = self.zeta_caras.copy()
        self._zc = 0.5 * (self._zf[1:] + self._zf[:-1])
        n = len(self._zc)

        VL_0 = self.V * self.LF
        IC = np.concatenate([[VL_0], np.full(n, cr.T_sat)])

        # Patrón del jacobiano. El sistema está acoplado globalmente:
        #  - dV_L/dt depende de todas las T (Q_wi = eta_w U_V sum P_j (T_air - T_j));
        #  - cada dT_j/dt depende de B_L (=> de todas las T vía Q_wi) y de dh/dt
        #    (término de grilla móvil), y en modo 'conservativo' además de todas las
        #    T_i con i < j por la recursión de G.
        # Por lo tanto el patrón correcto es denso (cualquier patrón más ralo
        # contamina el jacobiano por diferencias finitas agrupadas). Se usa
        # jac_sparsity=None (denso): n + 1 evaluaciones por jacobiano, que se
        # recalcula pocas veces con BDF/Radau.
        S = None
        self.jac_sparsity = S

        atol_V = self.atol_V if self.atol_V is not None else 1e-12 * self.V
        atol = np.concatenate([[atol_V], np.full(n, self.atol_T)])

        sol = solve_ivp(self.sys_isobaric, (0, t_f), IC, method=self.metodo,
                        rtol=self.rtol, atol=atol, jac_sparsity=S, dense_output=True)
        if not sol.success:
            raise RuntimeError("solve_ivp falló: " + sol.message)
        self.sol_completa = sol

        t_eval = np.arange(0, t_f + 1, self.time_interval)
        t_eval = t_eval[t_eval <= t_f]
        if t_eval[-1] < t_f:
            t_eval = np.append(t_eval, t_f)  # t_f siempre incluido (balances en t_f)
        self.sol = OptimizeResult(t=t_eval, y=sol.sol(t_eval), nfev=sol.nfev, njev=sol.njev,
                                  nlu=sol.nlu, message=sol.message, success=sol.success)
        self.sol.y[:, 0] = IC
        self._G_min = np.inf
        self._reconstruct()
        balances.calcular_balances(self)
        # Diagnóstico upwind (M2): el esquema supone G_f >= 0 en todas las caras
        self.data["G_min"] = float(self._G_min)
        if self._G_min < 0:
            warnings.warn("HorizontalTank: el flujo de masa G_f a través de una cara fue negativo "
                          "(G_min = %.3e kg/s). El esquema upwind supone flujo ascendente "
                          "(G >= 0); los resultados con G < 0 no son confiables." % self._G_min,
                          RuntimeWarning, stacklevel=2)

    # ---------------------------------------------------- ecuaciones del modelo
    def sys_isobaric(self, t, y):
        """Lado derecho del sistema: [dV_L/dt, dT_1/dt, ..., dT_N/dt] (m^3/s, K/s)."""
        return self._modelo(y)[0]

    def _modelo(self, y):
        """Evalúa el modelo para el estado y = [V_L, T_1..T_N].

        Devuelve (dy, d) con d un diccionario de diagnóstico (W, kg/s, m, ...).

        Líquido (docs §2):
            Q_L  = U_L A_wet(h) (T_air - T_sat)                                [W]
            Q_wi = eta_w U_V sum_j P_j (T_air - T_j)                           [W]
            Q_VL = k_sat A(h) dT/dz|_h (2.º orden, T(h) = T_sat)               [W]
            dV_L/dt = -(Q_L + Q_b + Q_VL + Q_wi) / (rho_L dH_LV)               [m^3/s]
            B_L = -rho_L dV_L/dt                                               [kg/s]
        Grilla móvil: dh/dt = (dV_L/dt)/A(h); velocidad de la cara f:
            zdot_f = (1 - zeta_f) dh/dt; dV_j/dt = A_j+1/2 zdot_j+1/2 - A_j-1/2 zdot_j-1/2.
        Flujo de masa relativo a la cara móvil, G_f [kg/s]:
            'conservativo': G_1/2 = B_L; G_j+1/2 = G_j-1/2 - d(rho_j V_j)/dt
                            (continuidad por celda; G en el techo = BOG exacto)
            'uniforme':     G_f = mdot - rho_f A_f zdot_f, con mdot uniforme en z en el
                            marco fijo más el término de grilla móvil (1 - zeta) dh/dt dT/dz
                            de docs §3. Se usa mdot = B_L (1 - rho_V_sat/rho_L), que es el
                            flujo de vapor en el marco fijo justo sobre la interfase (salto
                            de masa: el vapor además rellena el volumen que deja el líquido,
                            rho_V_sat A zdot_0 = -B_L rho_V_sat/rho_L). Así G_0 = B_L
                            exactamente. Con mdot = B_L literal (docs §3) se tendría
                            G_0 = B_L (1 + rho_V/rho_L), que viola el salto en la interfase
                            en una fracción rho_V/rho_L (~1,8 % en H2, ~0,1 % en NH3).
        Vapor por celda (forma conservativa en T, upwind):
            rho_j cp_j V_j dT_j/dt = G_j-1/2 [h(T_j-1) - h(T_j)] + F_j+1/2 - F_j-1/2
                                     + (1 - eta_w) U_V P_j (T_air - T_j) [+ Q_roof en j = N]
            F_f = k_f A(z_f) (T_j+1 - T_j)/dz (W, positivo hacia abajo); F_1/2 = Q_VL;
            F_N+1/2 = 0 (área nula en el techo horizontal; techo aislado en vertical).
        """
        cr, tb, geo = self.cryogen, self.tabla, self.geo
        V_L = y[0]
        T = y[1:]
        n = len(T)
        h = geo.h_desde_V(V_L)
        H = geo.H
        Lv = H - h
        zf = h + self._zf * Lv
        zf[-1] = H
        zc = h + self._zc * Lv
        A_f = geo.A(zf)
        V_j = geo.integral_A(zf[:-1], zf[1:])
        P_j = geo.integral_P(zf[:-1], zf[1:])

        # Propiedades
        rho = tb.rho(T)
        hT = tb.h(T)
        V_V = np.sum(V_j)
        if self.propiedades == "locales":
            cp = tb.cp(T)
            k = tb.k(T)
            k_f = 0.5 * (k[1:] + k[:-1])
            k_0 = self.k_V_sat
            dh_up = np.empty(n)
            dh_up[0] = self.h_V_sat - hT[0]
            dh_up[1:] = hT[:-1] - hT[1:]
            rho_cap = rho
        else:
            cp_m = np.sum(tb.cp(T) * V_j) / V_V
            k_m = np.sum(tb.k(T) * V_j) / V_V
            rho_m = np.sum(rho * V_j) / V_V
            cp = np.full(n, cp_m)
            k_f = np.full(n - 1, k_m)
            k_0 = k_m
            Tu = np.concatenate([[cr.T_sat], T[:-1]])
            dh_up = cp_m * (Tu - T)
            rho_cap = np.full(n, rho_m)

        # Conducción
        dzc = np.diff(zc)
        F = np.zeros(n + 1)  # F[f]: calor que baja por la cara f [W]
        F[1:-1] = k_f * A_f[1:-1] * (T[1:] - T[:-1]) / dzc
        d1 = zc[0] - h
        d2 = zc[1] - h
        dTdz_i = (T[0] - cr.T_sat) * d2 / (d1 * (d2 - d1)) - (T[1] - cr.T_sat) * d1 / (d2 * (d2 - d1))
        Q_VL = k_0 * A_f[0] * dTdz_i
        F[0] = Q_VL
        F[-1] = 0.0

        # Fuentes de pared
        Q_wall_j = self.U_V * P_j * (self.T_air - T)
        S = (1 - self.eta_w) * Q_wall_j
        S[-1] += self.Q_roof
        Q_wi = self.eta_w * np.sum(Q_wall_j)
        Q_L = self.U_L * geo.A_wet(h) * (self.T_air - cr.T_sat)
        Q_b = self.Q_b

        # Líquido
        dV_L = -(Q_L + Q_b + Q_VL + Q_wi) / (cr.rho_L * self.dH_LV)
        B_L = -cr.rho_L * dV_L

        # Grilla móvil
        dh = dV_L / A_f[0]
        zdot = (1 - self._zf) * dh
        AZ = A_f * zdot
        AZ[-1] = 0.0
        dV_j = AZ[1:] - AZ[:-1]

        C = F[1:] - F[:-1] + S          # entra por arriba - sale por abajo + fuente [W]
        M = rho_cap * cp * V_j          # capacidad térmica de cada celda [J/K]
        drho = tb.drho(T)
        G = np.empty(n + 1)
        if self.flujo_masa == "conservativo":
            G[0] = B_L
            dT = np.empty(n)
            for j in range(n):
                dT[j] = (G[j] * dh_up[j] + C[j]) / M[j]
                G[j + 1] = G[j] - (drho[j] * V_j[j] * dT[j] + rho[j] * dV_j[j])
        else:
            rho_f = np.concatenate([[self.rho_V_sat], rho[:-1], [rho[-1]]])
            mdot = B_L * (1 - self.rho_V_sat / cr.rho_L)  # marco fijo; da G_0 = B_L
            G[:] = mdot - rho_f * AZ
            dT = (G[:-1] * dh_up + C) / M
        dm_V = np.sum(drho * V_j * dT + rho * dV_j)  # dm_V/dt exacto [kg/s]
        BOG = B_L - dm_V

        if hasattr(self, "_G_min"):
            self._G_min = min(self._G_min, float(np.min(G)))
        dy = np.empty(n + 1)
        dy[0] = dV_L
        dy[1:] = dT
        d = dict(h=h, zf=zf, zc=zc, A_f=A_f, V_j=V_j, P_j=P_j, rho=rho, hT=hT, Q_VL=Q_VL, Q_wi=Q_wi, Q_L=Q_L, Q_b=Q_b,
                 Q_V=np.sum(Q_wall_j), Q_roof=self.Q_roof, B_L=B_L, BOG=BOG, dm_V=dm_V, G=G,
                 dV_L=dV_L, V_V=V_V, m_V=np.sum(rho * V_j), drho=drho, dT=dT, dV_j=dV_j,
                 k_V_avg=np.sum(tb.k(T) * V_j) / V_V, cp_V_avg=np.sum(tb.cp(T) * V_j) / V_V)
        return dy, d

    # ------------------------------------------------------------ resultados
    def _reconstruct(self):
        """Llena self.data en los tiempos de self.sol (claves de CryoEvap y nuevas).

        Claves: Time [s], V_L [m^3], h [m], LF [-], Tv_avg [K] (promedio ponderado por A),
        rho_V_avg [kg/m^3] (= m_V/V_V), m_V, m_L [kg], Q_VL, Q_L, Q_V (pared seca total),
        Q_Vw (= Q_wi), Q_b, Q_roof, Q_tot (Q_L + Q_b + Q_wi + Q_VL, como CryoEvap),
        Q_in (calor total que entra del ambiente) [W], B_L, BOG [kg/s], T_BOG [K]
        (celda del techo), drho_V_avg, dV_L (derivadas analíticas, no diferencias finitas).
        """
        claves = ["h", "Q_VL", "Q_wi", "Q_L", "Q_b", "Q_V", "B_L", "BOG", "dV_L", "V_V", "m_V",
                  "k_V_avg", "cp_V_avg"]
        acc = {c: [] for c in claves}
        Tv_avg, drho_avg, T_BOG = [], [], []
        for i in range(len(self.sol.t)):
            y = self.sol.y[:, i]
            _, d = self._modelo(y)
            for c in claves:
                acc[c].append(d[c])
            T = y[1:]
            Tv_avg.append(np.sum(T * d["V_j"]) / d["V_V"])
            T_BOG.append(T[-1])
            # rho_avg = m_V / V_V -> d rho_avg/dt = (dm_V - rho_avg dV_V)/V_V, dV_V = -dV_L
            drho_avg.append((d["dm_V"] + d["m_V"] / d["V_V"] * d["dV_L"]) / d["V_V"])
        a = {c: np.array(v) for c, v in acc.items()}
        self.data["Time"] = self.sol.t
        self.data["V_L"] = self.sol.y[0]
        self.data["h"] = a["h"]
        self.data["LF"] = self.sol.y[0] / self.V
        self.data["Tv_avg"] = np.array(Tv_avg)
        self.data["rho_V_avg"] = a["m_V"] / a["V_V"]
        self.data["m_V"] = a["m_V"]
        self.data["m_L"] = self.cryogen.rho_L * self.sol.y[0]
        self.data["Q_VL"] = a["Q_VL"]
        self.data["Q_L"] = a["Q_L"]
        self.data["Q_V"] = a["Q_V"]
        self.data["Q_Vw"] = a["Q_wi"]
        self.data["Q_b"] = a["Q_b"]
        self.data["Q_tot"] = a["Q_L"] + a["Q_b"] + a["Q_wi"] + a["Q_VL"]
        self.data["Q_in"] = a["Q_L"] + a["Q_b"] + a["Q_V"] + self.Q_roof
        self.data["B_L"] = a["B_L"]
        self.data["BOG"] = a["BOG"]
        self.data["T_BOG"] = np.array(T_BOG)
        self.data["drho_V_avg"] = np.array(drho_avg)
        self.data["dV_L"] = a["dV_L"]
        # Actualiza el Cryogen con el último estado (como CryoEvap)
        self.cryogen.T_V = self.sol.y[1:, -1]
        self.cryogen.Tv_avg = Tv_avg[-1]
        self.cryogen.rho_V_avg = self.data["rho_V_avg"][-1]
        self.cryogen.k_V_avg = a["k_V_avg"][-1]
        self.cryogen.cp_V_avg = a["cp_V_avg"][-1]

    def diagnostico(self, t):
        """Cantidades por celda y globales del estado en el tiempo t [s] (0 <= t <= t_f).

        Si t es uno de los tiempos guardados (tank.sol.t) se usa ese estado; si no,
        la salida densa del integrador. No cambia el estado del estanque.

        Devuelve un diccionario con (unidades SI):
            t [s], V_L [m^3], h [m] nivel del líquido
            z_f   (N+1,) caras físicas [m], de h a H
            z_c   (N,)   centros de celda [m]
            A_f   (N+1,) área de rebanada en las caras [m^2] (A_f[0] = A(h) interfase)
            V_j   (N,)   volumen exacto de cada celda, int A dz [m^3]
            P_j   (N,)   área de pared exacta de cada celda, int P dz [m^2]
            T_j   (N,)   temperatura de cada celda [K]
            rho_j (N,)   densidad local [kg/m^3];  h_j (N,) entalpía local [J/kg]
            G     (N+1,) flujo de masa relativo a cada cara móvil [kg/s] (G[0] = B_L;
                         G[-1] = BOG en modo 'conservativo')
            dT_dt (N,) [K/s], dV_L_dt [m^3/s]
            Q_L   = U_L A_wet(h) (T_air - T_sat)                [W]
            Q_wi  = eta_w U_V sum_j P_j (T_air - T_j)           [W]
            Q_VL  = k A(h) dT/dz|_h                             [W]
            Q_V   = U_V sum_j P_j (T_air - T_j)  (pared seca)   [W]
            Q_Vv  = (1 - eta_w) Q_V  (parte que entra al vapor) [W]
            Q_b, Q_roof [W]
            Q_in  = Q_L + Q_b + Q_V + Q_roof (calor del ambiente) [W]
            B_L, BOG, dm_V_dt [kg/s], m_V, m_L [kg]
        """
        if self.sol is None:
            raise ValueError("No hay solución: ejecuta tank.evaporate(t_f) primero")
        t = float(t)
        k = np.flatnonzero(np.isclose(self.sol.t, t, rtol=0, atol=1e-9 * max(1.0, self.sol.t[-1])))
        if k.size:
            y = self.sol.y[:, k[0]]
        else:
            if t < 0 or t > self.sol_completa.t[-1]:
                raise ValueError("t = %g s fuera de [0, %g] s" % (t, self.sol_completa.t[-1]))
            y = self.sol_completa.sol(t)
        dy, d = self._modelo(y)
        Q_in = d["Q_L"] + d["Q_b"] + d["Q_V"] + d["Q_roof"]
        return dict(t=t, V_L=float(y[0]), h=d["h"], z_f=d["zf"], z_c=d["zc"], A_f=d["A_f"],
                    V_j=d["V_j"], P_j=d["P_j"], T_j=np.array(y[1:]), rho_j=d["rho"], h_j=d["hT"],
                    G=d["G"], dT_dt=dy[1:], dV_L_dt=dy[0], Q_L=d["Q_L"], Q_wi=d["Q_wi"],
                    Q_VL=d["Q_VL"], Q_V=d["Q_V"], Q_Vv=(1 - self.eta_w) * d["Q_V"], Q_b=d["Q_b"],
                    Q_roof=d["Q_roof"], Q_in=Q_in, B_L=d["B_L"], BOG=d["BOG"], dm_V_dt=d["dm_V"],
                    m_V=d["m_V"], m_L=self.cryogen.rho_L * float(y[0]))

    def evap_rate(self):
        """Tasa de evaporación B_L [kg/s] reconstruida (igual que data['B_L'])."""
        if self.sol is None:
            print("No solution available, perform tank.evaporate() to construct solution")
            return
        return self.data["B_L"]
