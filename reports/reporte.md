# Reporte del Proyecto 12: evaporación isobárica en estanques cilíndricos horizontales

Autor de las figuras y del reporte: subagente graficador. Fecha: 2026-11-16.
Todo se generó con `python3 scripts/barrido.py` (≈ 40 s con 4 procesos; 113 simulaciones). No se modificó `src/` ni `tests/`.

**Dos advertencias que condicionan todo el reporte**

1. **La validación es con hidrógeno; la aplicación es con amoníaco.** El caso 2 (H2, Wang y Mérida 2024) es la única comparación horizontal contra una referencia externa. No existe ninguna validación del modelo con amoníaco.
2. **El caso 3 es PROVISIONAL.** Aún no tiene dimensiones reales: se usan las del caso 2 (d_i = d_o = 1.676 m, L = 2.217 m, V = 4.89 m³, tapas planas), U equivalente a 0.5 W/m² a T_air = 293.15 K (U = 0.0094 W/m²/K), η_w = 0.97 y P = 101325 Pa. Cada figura del caso 3 lo indica al pie. Las temperaturas ambiente 288.15 K (sur) y 303.15 K (norte) son valores **ilustrativos**, no datos climáticos.

## 1. Modelo y resumen de cambios

Archivos creados (nada más fue tocado):

| Archivo | Contenido |
|---|---|
| `scripts/barrido.py` | Barrido de llenados, validación, sensibilidad, comparación de fluidos; escribe `figs/`. Función `simular(cfg)` sobre `HorizontalTank` y `diagnostico(t)`. `--usar-cache` rehace solo figuras (usa `scripts/.cache_barrido.pkl`, regenerable). |
| `figs/fig01..fig07_*.pdf/.png/.csv` | 7 figuras (PDF vectorial, PNG de 200 dpi, CSV con los datos) |
| `figs/tabla_conservacion.csv`, `figs/tabla_malla.csv` | Errores de conservación por corrida y convergencia de malla |
| `reports/reporte.md` | Este reporte |

Definiciones usadas: **BOR [%/día] = 100 · BOG · 86400 / m_L**, evaluado a t = 24 h (fuera del transitorio inicial; ver fig. 4) y también a 200 h en los CSV. Simulaciones de 200 h, malla de 40 celdas, salida cada 1 h (600 s en el caso 2).
**Vertical del mismo volumen:** misma d_i, altura H = V/(πR²) = 2.217 m (coincide con L por usar tapas planas), mismo U. Se corren dos variantes: (i) techo adiabático, convención de CryoEvap (Q_roof = 0; el fondo sí recibe calor U·A_T·ΔT); (ii) "con techo", Q_roof = U·A_T·(T_air − T_sat), una **cota superior** del calor por el techo, que iguala el área de envolvente total del horizontal.

## 2. Implementación y conservación

### Errores de conservación (todas las corridas del barrido)

| Fluido / geometría | Corridas | máx. error masa | máx. error energía | mín. G_f [kg/s] |
|---|---|---|---|---|
| NH3 horizontal | 9 (LF 0.1–0.9) | 1.4e-6 | 1.9e-7 | 5.7e-6 |
| NH3 vertical (techo adiabático) | 9 | 5.3e-7 | 3.7e-8 | 4.9e-6 |
| NH3 vertical con techo | 9 | 8.8e-8 | 5.3e-8 | 4.2e-6 |
| N2 horizontal | 9 | 1.4e-7 | 8.4e-8 | 1.6e-4 |
| N2 vertical (adiabático / con techo) | 18 | 5.4e-8 | 2.7e-8 | 1.4e-4 |
| Caso 2 (H2) | 1 | 1.3e-8 | 1.1e-8 | 1.7e-5 |
| Sensibilidad NH3 (54 corridas) | 54 | 1.7e-5 | 2.4e-6 | todas > 0 |

Máximo global: 1.7e-5 (masa), 2.4e-6 (energía), frente al criterio < 1e-2. No hubo avisos de flujo negativo (G_min > 0 en las 113 corridas), de modo que el esquema upwind es válido.

**Qué significan estos errores.** Miden la consistencia del esquema y del integrador, no validan la física: en modo conservativo el semidiscreto cierra masa y energía casi por construcción (docstring de `balances.py`, registro de decisiones). Que sean ~1e-8 no dice nada sobre si el modelo reproduce la realidad. La evidencia física son la regresión contra CryoEvap (caso 1) y la comparación con Wang y Mérida (caso 2).

### Convergencia de malla (NH3, LF = 0.3, horizontal, 200 h; `figs/tabla_malla.csv`)

| n celdas | BOG 24 h [kg/h] | BOG 200 h [kg/h] | T_v media 200 h [K] | T techo 200 h [K] |
|---|---|---|---|---|
| 20 | 0.0209344 | 0.0208711 | 240.4865 | 241.0645 |
| 40 | 0.0209339 | 0.0208705 | 240.4869 | 241.0753 |
| 80 | 0.0209337 | 0.0208701 | 240.4872 | 241.0806 |
| 160 | 0.0209337 | 0.0208700 | 240.4873 | 241.0833 |

Entre 80 y 160 celdas el BOG cambia 4e-6 relativo (criterio < 5e-3): la malla de 40 celdas usada en el barrido está sobrada para el BOG. La temperatura del techo, que es lo más sensible a la resolución, varía < 0.02 K entre 20 y 160 celdas: **los artefactos numéricos cerca del techo son despreciables** frente a las diferencias físicas discutidas abajo.

## 3. Resultado de los tests (`python3 -m pytest -q`)

81 pasan, 2 fallan (83 en total; otro agente agregó tests de verificación independiente, que no se tocaron).

| Grupo (archivo) | Resultado |
|---|---|
| Conservación de masa y energía (`test_conservacion.py`) | pasan |
| Convergencia de malla (`test_convergencia_malla.py`) | pasan |
| Regresión vertical frente a CryoEvap, caso 1 (`test_limite_vertical.py`) | pasan (ver nota b) |
| Geometría, llenado fuera de rango, sin calor al vapor y otros | pasan |
| `test_wang_merida.py::test_wang_merida_estacionario` | **FALLA**: dif. media de BOG (t ≥ 15 h) = 25.5 % (modelo 0.0627, Wang 0.0841 kg/h), criterio ≤ 10 % |
| `test_wang_merida.py::test_wang_merida_curva_completa` | **FALLA**: dif. media de BOG, curva completa = 37.8 % |

Los tests de conservación y de regresión vertical pasan, por lo que corrí el barrido (regla de graficador.md). Las dos fallas de Wang y Mérida no se ocultan: se discuten en la sección 4.

Notas de honestidad (hallazgos del revisor):
- (a) Los errores de masa y energía no validan la física (sección 2).
- (b) El caso 1 pasa contra CryoEvap (< 1 %) porque η_w = 0.963 hace pesar poco al vapor. Con η_w = 0 el BOG difiere ~20 % por diferencias de modelo de CryoEvap (v_z con ρ_sat, término de grilla, Q_wi con Tv_avg, BOG por diferencias finitas); el modelo nuevo es más consistente, pero esa coincidencia no es una validación fuerte del vapor.
- (c) Caso 2: ver abajo.

## 4. Validación y sensibilidad

### Fig. 1: geometría (`figs/fig01_geometria_pared_mojada.pdf`)
La fracción de pared mojada en horizontal sigue una curva en S frente a LF: queda sobre la línea A_wet/A_total = LF a LF < 0.5 y bajo ella a LF > 0.5; en vertical (manto + fondo) es casi lineal y siempre incluye el fondo. A igual LF la diferencia es pequeña (< 0.07 en fracción). Lo que sí difiere mucho es la **interfase**: horizontal 1.0–3.7 m² según LF (máx. a LF = 0.5), vertical constante 2.2 m². Es un resultado geométrico exacto, sin malla ni física.

### Fig. 2: validación con hidrógeno, caso 2 (`figs/fig02_validacion_caso2_h2.pdf`)
**Es el resultado negativo del proyecto y se muestra completo.** El BOG del modelo se estabiliza en 0.0623 kg/h frente a 0.084 kg/h de Wang y Mérida (−26 %, no cumple el 10 %); en el transitorio inicial Wang muestra un BOG de 0.18 kg/h que el modelo no reproduce (el modelo parte con vapor saturado). Hallazgo del revisor, confirmado con la hoja 12d: la **tasa inicial de pérdida de líquido del modelo coincide con la de Wang en t = 0** (0.0633 frente a 0.0634 kg/h), pero en Wang sube ~35 % en unas 4 h (hasta 0.0856 kg/h) mientras el modelo queda en 0.0635 kg/h. Por eso el llenado de Wang cae a 25.2 % a las 200 h y el del modelo a 26.3 %. La temperatura del vapor de Wang (≈ 20.93 K al final) queda entre la T media del modelo (20.63 K) y la T del techo (20.95 K).
Causa propuesta (hipótesis, no verificada): en Wang la entrada de calor crece con el tiempo (p. ej. un coeficiente o un flujo que depende de la temperatura de pared, o un área efectiva mayor) o el flujo de 0.5 W/m² se define de otra forma. **No se ajustó ningún parámetro.** Falta revisar en el paper cómo se define la entrada de calor. Hasta entonces, el modelo horizontal no está validado cuantitativamente.

### Fig. 3: BOR frente a LF, horizontal frente a vertical, amoníaco (`figs/fig03_bor_vs_lf_nh3.pdf`), figura clave
El BOR baja de 0.151 %/día (LF = 0.1) a 0.017 %/día (LF = 0.9) en horizontal, casi como 1/LF, porque el calor que entra es casi independiente de LF y se reparte en más masa. Frente al vertical con techo adiabático, el horizontal tiene **~15 % más BOR en todo el rango** (razón 1.148–1.158). Esa diferencia se explica casi por completo por el área: con techo adiabático el vertical recibe 6.88 W (LF = 0.3) frente a 7.97 W del horizontal, es decir, el horizontal "tiene" el área del techo. Al darle calor por el techo al vertical (cota superior), la razón cae a **1.00–1.07** (máx. ~1.065 a LF = 0.5–0.6): la geometría horizontal por sí sola aporta como mucho ~6 % de BOR adicional. La conclusión depende fuertemente de qué se suponga para el techo vertical; la verdad está entre ambas curvas.

### Fig. 4: transitorio del BOR (`figs/fig04_bor_vs_tiempo_nh3.pdf`)
El BOR instantáneo cambia poco en 200 h: cae ~3–5 % en las primeras ~10 h (el vapor se estratifica y absorbe calor sensible) y luego es casi constante; a LF bajo, sube levemente al final porque disminuye m_L. Por eso el BOR a 24 h es representativo. Con η_w = 0.97 casi todo el calor va directo a la interfase y el vapor se calienta apenas ~1 K, así que el transitorio térmico del vapor pesa poco en el BOR.

### Fig. 5: perfiles de temperatura del vapor (`figs/fig05_perfiles_Tv_nh3.pdf`)
El vapor se estratifica: T_sat = 239.8 K en la interfase y hasta 241.1 K en el techo (LF = 0.3, 200 h); a LF = 0.8 solo 240.2 K. El perfil llega a régimen casi estacionario a las 50 h (curvas de 50 y 200 h superpuestas) para LF = 0.3–0.8. A LF = 0.8 el espacio de vapor es una cuña delgada y el gradiente es menor. Rango de temperaturas muy pequeño (≈ 1 K): el vapor de amoníaco queda casi saturado en este caso, lo cual es consecuencia directa de η_w = 0.97 (calibrado con hidrógeno, **no validado para amoníaco**). Efecto numérico: la celda del techo es estrecha, pero la T del techo cambia < 0.02 K con la malla.

### Fig. 6: sensibilidad, amoníaco horizontal (`figs/fig06_sensibilidad_nh3.pdf`)
- **U:** el BOR es casi proporcional (×0.50 y ×1.97–2.00 para U×0.5 y U×2), es decir, domina el aporte de calor.
- **T_air:** 288.15 K (sur, ilustrativo) da −9.4 % y 303.15 K (norte, ilustrativo) +18.7 % de BOR, independiente de LF; coincide con ΔT = T_air − T_sat (48.4/53.4 y 63.4/53.4). U se mantuvo fijo.
- **η_w:** con η_w = 0 el BOR baja 28 % a LF = 0.1 y solo 2.5 % a LF = 0.9 (η_w = 0.5: −14 % y −1 %). Es el parámetro más incierto (no tiene respaldo para NH3) y su efecto es mayor a LF bajo, donde hay más pared seca.

### Fig. 7: nitrógeno frente a amoníaco en el mismo estanque (`figs/fig07_comparacion_n2_nh3.pdf`)
Con la misma U absoluta (misma aislación) el N2 recibe ~4 veces más calor (ΔT = 216 K frente a 53 K) y tiene un calor latente volumétrico ~5.8 veces menor, de modo que su BOR es ~23 veces el del amoníaco (p. ej. 1.18 frente a 0.050 %/día a LF = 0.3; BOG 0.58 frente a 0.021 kg/h). La razón horizontal/vertical es casi la misma en ambos fluidos (~1.15 con techo adiabático): es una propiedad de la geometría y del área de envolvente, no del fluido. Con el techo vertical calentado (cota superior) el N2 apenas cambia respecto del techo adiabático (razón H/V ≈ 1.15 en casi todo el rango, 1.13 a LF = 0.9), mientras que en NH3 sí cambia (1.00–1.07). No investigué la causa; es un punto a revisar (hipótesis: el calor inyectado en la celda del techo llega a la interfase solo por conducción, lo que depende de k y de los gradientes del fluido). La comparación con N2 horizontal no tiene validación propia (solo regresión vertical, caso 1).

## 5. Limitaciones y próximos pasos

Limitaciones:
- Caso 2 abierto: BOG −26 % bajo Wang y Mérida; sin causa confirmada (hipótesis en sección 4).
- Caso 3 provisional (dimensiones, U, T_air, tapas); las conclusiones cuantitativas de BOR no son aplicables a un estanque real hasta definirlo con Amalia y el profesor.
- η_w = 0.97 viene del caso de hidrógeno; no hay respaldo para amoníaco y el BOR a LF bajo depende de él.
- La comparación horizontal/vertical depende de cómo se trate el techo vertical (razón 1.15 frente a 1.00–1.07); falta fijar una convención justa y aprobada.
- Los errores de conservación no validan la física. El caso 1 pasa por el peso pequeño del vapor con η_w = 0.963 (con η_w = 0 difiere ~20 % por diferencias de modelo de CryoEvap).
- Modelo 1D con vapor sin convección natural ni efectos de la fase líquida (líquido bien mezclado a T_sat), supuesto isobárico con venteo; la validez para la operación real del amoníaco horizontal sigue pendiente (AGENTS.md).
- Las T_air norte/sur son ilustrativas.

Próximos pasos:
1. Revisar en Wang y Mérida (2024) cómo se define y cómo evoluciona la entrada de calor, y repetir el caso 2 antes de cualquier ajuste.
2. Obtener las dimensiones reales, U, T_air y tipo de tapas del caso 3; rehacer figs. 3, 5, 6 y 7 (basta cambiar constantes en `scripts/barrido.py`).
3. Decidir el tratamiento del calor por el techo en el vertical de comparación.
4. Estudiar η_w para amoníaco (literatura o ensayo) y, si es posible, buscar un dato con amoníaco para validar.
5. Actualizar docs/ecuaciones.md §3 y §6 (forma conservativa con G(z)), con aprobación de Amalia.
