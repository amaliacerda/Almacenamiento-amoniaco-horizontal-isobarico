# Ecuaciones del modelo: estanque cilíndrico horizontal isobárico

Referencia de ecuaciones para todos los agentes. Derivación paso a paso
explicada en el documento "Ecuaciones del estanque horizontal paso a paso".
Unidades SI. Cambiar este archivo requiere aprobación de Amalia (vía orquestador).

## Notación
- R: radio interior; D = 2R; L: largo de la parte cilíndrica; d_o: diámetro exterior.
- z: altura desde el fondo (z = 0 fondo, z = 2R techo). h: nivel del líquido.
- w(z) = sqrt(z (2R - z)): media cuerda. theta(h) = acos((R - h)/R).
- Tapas planas salvo que AGENTS.md diga otra cosa.

## 1. Geometría
- Área de una rebanada horizontal: A(z) = 2 L w(z). Área de interfase: A(h).
- Volumen de líquido: V_L(h) = L [R^2 theta - (R - h) sqrt(2Rh - h^2)].
  h se obtiene de V_L numéricamente (brentq).
- Pared por unidad de altura: P(z) = 2 L R / w(z) + 4 w(z)
  (manto + dos tapas planas). P es singular en z = 0 y z = 2R, pero su integral es finita.
- Área mojada: A_wet(h) = 2 L R theta + 2 [R^2 theta - (R - h) sqrt(2Rh - h^2)].
- Área seca: A_dry = 2 pi R L + 2 pi R^2 - A_wet(h).
- Integrales por celda [z0, z1] (usar estas, no evaluar A o P en un punto):
  - volumen: integral de A(z) dz
  - pared: integral de P(z) dz = 2 L R [phi(z1) - phi(z0)] + 2 [S(z1) - S(z0)],
    con phi(z) = acos((R - z)/R) y S(z) = R^2 phi(z) - (R - z) w(z).
- Con tapas hemisféricas: sumar a A(z) el término pi w(z)^2, al volumen el de
  la esfera cortada, y a la pared el área de la esfera cortada (2 pi R por unidad de altura).

## 2. Líquido (saturado a T_sat(P), bien mezclado)
  dV_L/dt = -(Q_L + Q_VL + Q_wi) / (rho_L dH_LV)
  B_L     = (Q_L + Q_VL + Q_wi) / dH_LV
  Q_L     = U_L A_wet(h) (T_air - T_sat)
  Q_VL    = k A(h) dT/dz |_(z = h)
  Q_wi    = eta_w U_V integral_h^2R P(z) (T_air - T) dz
  Q_b     = 0 (no hay "fondo" separado en horizontal)

## 3. Vapor (z entre h y 2R)
Forma conservativa (la que se implementa, en volúmenes finitos):
  rho cp A dT/dt = d/dz(k A dT/dz) - mdot cp dT/dz + (1 - eta_w) U_V P(z) (T_air - T)
- mdot: flujo de masa ascendente; caso base mdot = B_L(t), uniforme en z.
- Bordes: T(h) = T_sat; en z = 2R el área de la cara es 0, así que no hay flujo.
- Grilla móvil: con zeta = (z - h)/(2R - h) aparece el término (1 - zeta) dh/dt dT/dz,
  con dh/dt = (dV_L/dt) / A(h).
- Propiedades del vapor: locales o promedios ponderados por A(z) (decisión en AGENTS.md).

Forma dividida por rho cp A (solo para explicar; NO implementar así, diverge en z = 2R):
  dT/dt = alpha [d2T/dz2 + (A'/A) dT/dz] - (mdot / (rho A)) dT/dz
          + (1 - eta_w) U_V (P/A) (T_air - T) / (rho cp)
  A'/A = (R - z) / (z (2R - z)),   P/A = R / (z (2R - z)) + 2 / L

## 4. BOG
  B = B_L - dm_V/dt,  m_V = integral_h^2R rho A dz = rho_avg V_V
  B = B_L + rho_avg dV_L/dt - V_V d(rho_avg)/dt

## 5. Balances globales (verificación, error < 1 %)
- Masa: [m_L(0) + m_V(0)] - [m_L(t) + m_V(t)] = integral_0^t B dt.
- Energía (P y volumen constantes, dH = dU):
  H(t) - H(0) + integral_0^t B h_V(T_techo) dt = integral_0^t Q_total dt,
  con H = m_L h_L + integral rho h(T) A dz.

## 6. Límite vertical (regresión)
Con A(z) = pi d_i^2 / 4 constante y P(z) = pi d_o, las ecuaciones se reducen al
modelo de CryoEvap (rama main): S_wall = 4 U_V d_o / d_i^2 (T_air - T)(1 - eta_w).
