# CLAUDE.md

Claude Code carga este archivo al abrir el repositorio. La sesión principal
actúa como el orquestador del proyecto.

@AGENTS.md
@docs/ecuaciones.md

## Rol de la sesión principal
Sigue las instrucciones de `orquestador.md` (copiado en `.claude/orquestador.md`).
Los subagentes están en `.claude/agents/`: modelador, tester-fisico, graficador
y revisor-numerico. Delega en ellos con la herramienta de agentes; no hagas su trabajo.

## Instalación de los agentes en el repositorio
1. Copiar `AGENTS.md`, `CLAUDE.md` y `docs/ecuaciones.md` a la raíz del repositorio
   amaliacerda/Almacenamiento-amoniaco-horizontal-isobarico.
2. Copiar `modelador.md`, `tester_fisico.md`, `graficador.md` y
   `revisor_numerico.md` a `.claude/agents/`.
3. Copiar `orquestador.md` a `.claude/orquestador.md` (no va en `.claude/agents/`
   porque es la sesión principal, no un subagente).
