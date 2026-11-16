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
