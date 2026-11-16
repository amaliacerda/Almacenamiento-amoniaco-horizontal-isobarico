# AGENTS.md: reglas y parámetros comunes del Proyecto 12

Proyecto 12 (IIQ3843, 2026): modelo 1D de evaporación isobárica de amoníaco en
estanques cilíndricos horizontales. Este archivo es la fuente única de verdad
para todos los agentes. Si un agente necesita un parámetro que no está aquí,
lo pide al orquestador; nunca lo inventa.

## Repositorios
- Repositorio de trabajo (el único donde se escribe):
  https://github.com/amaliacerda/Almacenamiento-amoniaco-horizontal-isobarico
- Referencia, SOLO LECTURA: https://github.com/felipehuerta17/CryoEvap
  - Rama `main` (commit 859e9ba): modelo vertical. Es la BASE del modelo nuevo.
  - Rama `HorizontalTank`: no se usa como base. Solo se puede leer el archivo
    `data/Data_wang_merida.xlsx` para validar.
  - Nunca hacer push, PR, issues ni comentarios en el repositorio de CryoEvap.
- Licencia: CryoEvap es GPL-3.0. Si se copia o adapta código suyo, el
  repositorio de trabajo debe ser GPL-3.0 y citar el original.

## Decisión de diseño
Se parte desde cero sobre el modelo vertical de CryoEvap (`main`): se instala
CryoEvap sin modificarlo y se escribe una clase propia `HorizontalTank` que
reutiliza `Cryogen` (propiedades con CoolProp) y la estructura de `Tank`
(`evaporate`, `data`, `_reconstruct`), cambiando toda la geometría y la
ecuación del vapor. CryoEvap vertical sin cambios sirve como caso de regresión.

## Estructura del repositorio de trabajo
```
src/horizontal/      modelador      geometria.py, tanque_horizontal.py, balances.py
tests/               tester-fisico  test_*.py
scripts/ figs/ reports/  graficador
data/                Data_wang_merida.xlsx (copiado de la rama HorizontalTank)
docs/ecuaciones.md   referencia de ecuaciones (todos la leen, nadie la cambia sin aprobación)
AGENTS.md, CLAUDE.md
.claude/agents/      modelador.md, tester_fisico.md, graficador.md, revisor_numerico.md
```

## Entorno
conda Cryosoft: python 3, numpy, scipy, matplotlib, CoolProp, openpyxl, pytest.
CryoEvap instalado con `pip install git+https://github.com/felipehuerta17/CryoEvap.git@859e9ba`.

## Convenciones
- Unidades SI en todo el código (Pa, K, m, kg, s, W, J).
- z se mide desde el fondo del estanque; z = 2R es el techo. h = nivel del líquido.
- Calores positivos cuando entran al estanque o al líquido.
- Nombre CoolProp del amoníaco: "Ammonia"; del nitrógeno: "Nitrogen".

## Casos base

### Caso 1: regresión vertical (CryoEvap original, notebook LN2_small)
Nitrógeno, P = 100000 Pa, d_i = 0.201 m, d_o = 0.204 m, V = 6.75e-3 m^3,
LF = 0.278, U_L = U_V = 0.026 W/m^2/K, T_air = 298.15 K, Q_b = 0, Q_roof = 0,
eta_w = 0.963, dz = 0.005 m, t = 2 * tau (como en el notebook).

### Caso 2: validación horizontal (Wang y Mérida 2024, figura 12)
Hidrógeno (solo para validar), P = 101300 Pa, V = 4.89 m^3, d_i = d_o = 1.676 m,
L = 2.217 m, tapas planas, LF = 0.30, flujo de calor 0.5 W/m^2
(U = 0.5 / (T_air - T_sat)), eta_w = 0.97, T_air = 293.15 K (supuesto),
t = 200 h. Datos: data/Data_wang_merida.xlsx (hojas 12a_T, 12b_LF, 12c_mL,
12d_dMLdt, 12e_BOG). Son curvas digitalizadas.

### Caso 3: aplicación (amoníaco, estanque horizontal)
Amoníaco, P = 101325 Pa (refrigerado, isobárico con venteo).
Dimensiones, U y T_air: POR DEFINIR con Amalia y el profesor.
Mientras tanto, para pruebas: mismas d_i, L y U equivalente del caso 2.

## Fluidos
- Hidrógeno: solo caso 2 (validación).
- Nitrógeno: comparación y regresión.
- Amoníaco: fluido de aplicación.

## Criterios de aceptación
- Conservación de masa y de energía: error relativo < 1 % (exigencia del curso).
- Caso 1: el código nuevo en modo vertical reproduce CryoEvap original
  (diferencia en BOG y Tv_avg < 1 %).
- Caso 2: diferencia en BOG frente a Wang y Mérida <= 10 %.
- Convergencia de malla: el BOG cambia < 0.5 % entre las dos mallas más finas.

## Decisiones pendientes
- Confirmar con el profesor el alcance y la aplicación (operación real del
  amoníaco horizontal y validez del supuesto isobárico).
- Dimensiones reales del caso 3.
- Tipo de tapas del caso 3 (planas por defecto).

## Registro de decisiones
(El orquestador anota aquí cada decisión de modelado: fecha, decisión, motivo.)

- 2026-11-16. Flujo de masa del vapor conservativo por defecto (`flujo_masa='conservativo'`):
  G se obtiene celda a celda de la continuidad local, G = B_L en la interfase y G = BOG en el techo.
  La versión uniforme mdot = B_L de docs/ecuaciones.md queda como `flujo_masa='uniforme'`.
  Motivo: con mdot uniforme el balance de energía falla en 1,47 % (H2, eta_w = 0); conservativo, 7e-9.
- 2026-11-16. Propiedades del vapor locales por celda por defecto (`propiedades='locales'`);
  la ponderación por A(z) queda como `propiedades='promedio'` (error de energía 0,55 % en la misma prueba).
  Propiedades tabuladas con CoolProp a P fija y splines (`TablaVapor`), cp = dh/dT.
- 2026-11-16. Advección escrita en entalpías, upwind: G [h(T_j-1) - h(T_j)].
- 2026-11-16. Gradiente en la interfase de segundo orden (T_sat, T_1, T_2), mismo para Q_VL y la celda 1.
- 2026-11-16. Manto con R exterior (2 L R_o / w), como pi d_o en CryoEvap; tapas con R interior.
  Con d_o = d_i (casos 2 y 3) coincide con docs/ecuaciones.md.
- 2026-11-16. Jacobiano denso por diferencias finitas (jac_sparsity=None): el acoplamiento es global (B_L y Q_wi dependen
  de todas las T), así que el patrón tridiagonal + columna de V_L era incompleto. Resultados idénticos, tiempo ~1 s.
- 2026-11-16. Modo 'uniforme': mdot = B_L (1 - rho_V_sat/rho_L) para cumplir el salto de masa en la interfase
  (G_0 = B_L). Se aparta en rho_V/rho_L (~1,8 % H2, ~0,1 % NH3) del mdot = B_L de docs §3.
- 2026-11-16. Diagnóstico: `HorizontalTank.diagnostico(t)` entrega cantidades por celda y calores separados; data['G_min']
  y RuntimeWarning si G < 0 (upwind); t_f se incluye siempre en t_eval.
- 2026-11-16. Paquete instalable con pyproject.toml (layout src, `pip install -e .`). Licencia GPL-3.0 (LICENSE) por adaptar CryoEvap.
- 2026-11-16. Caso 2 abierto: con tapas planas y 0,5 W/m^2 el BOG queda ~26 % bajo Wang y Mérida (Q_in 8,0 W frente a ~10,6 W
  implícitos en su régimen). El revisor descarta un error de implementación y descarta el área como causa principal: la
  tasa inicial de pérdida de líquido del modelo (0,0634 kg/h) coincide con la de Wang en t = 0 (0,0634 kg/h, hoja 12d), y
  en Wang esa tasa sube ~35 % en 4 h. Hipótesis: en Wang la entrada de calor crece en el tiempo o el flujo se define de
  otra forma. No se ajustó ningún parámetro; falta revisar en el paper cómo se define la entrada de calor.
- 2026-11-16. Los errores de masa y energía cierran casi por construcción (miden consistencia del esquema y del
  integrador, no validan la física). Se agregan verificaciones independientes (Q_in con fórmulas cerradas, suma de
  P_j = A_dry, límite eta_w = 1, inyección de errores).
- 2026-11-16. Modo vertical frente a CryoEvap: pasa el caso 1 (< 1 %) porque eta_w = 0,963 hace pesar poco al vapor.
  Con eta_w = 0 el BOG difiere ~20 % por diferencias de modelo de CryoEvap (v_z con rho_sat, término de grilla, Q_wi con
  Tv_avg, BOG por diferencias finitas). El modelo nuevo es más consistente. docs/ecuaciones.md §3 y §6 deben actualizarse
  (forma conservativa con G(z)) con aprobación de Amalia.
