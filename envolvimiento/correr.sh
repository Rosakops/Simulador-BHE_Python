#!/bin/sh
# =============================================================================
#  Lanzador del proyecto BHE / SENACYT
# =============================================================================
#  USO
#    sh correr.sh            pruebas + veredictos + desglose + figuras
#    sh correr.sh pruebas    solo el marcador de las 4 suites
#    sh correr.sh tabla      solo veredictos y desglose
#    sh correr.sh figuras    solo regenerar los PNG
#    sh correr.sh dendrimero ventana geométrica del PAMAM + sus 6 figuras
#    sh correr.sh teoricos   dendrímeros + polímeros de las fichas (SINTÉTICOS)
#    sh correr.sh polimero  suelo del polímero macizo y la micela (G.1b)
#    sh correr.sh web       genera la web estática en ../web  (no calcula: muestra)
#    sh correr.sh notas      la prosa: lectura de resultados y qué cambió
#    sh correr.sh farmaco    tamano molecular del fingolimod (tarea C9)
#    sh correr.sh colab      regenerar el Colab y verificar equivalencia
#    sh correr.sh cribado    ../cribado_bhe_tres_liposomas.py  (tiene el bug P5)
#    sh correr.sh acople     ../acople/acople_fingolimod_liposoma.py
#
#  DESACTIVADOS el 2026-08-24 (decisión de Jhovan): 'dendrimero', 'teoricos' y
#  'polimero' de la lista de arriba ya NO son subcomandos activos — el
#  simulador es exclusivo para liposomas. Se dejan comentados como registro
#  histórico de cuando el simulador incluía esas clases; ver más abajo.
#
#  POR QUÉ ESTE ARCHIVO EXISTE
#  Las tareas de Zed descartan el campo `args` y parten los comandos con
#  comillas. Metiendo todo aquí, Zed solo pasa una palabra y no hay nada que
#  romper. Ver la cabecera de ../.zed/tasks.json.
#
#  MPLBACKEND=Agg es obligatorio: envolvimiento_colab.py se ejecuta entero al
#  importarse y acaba en plt.show(), que en una máquina con escritorio abriría
#  una ventana y dejaría el script colgado.
#  Para verlas de forma interactiva:  MPLBACKEND= sh correr.sh figuras
# =============================================================================

set -e
cd "$(dirname "$0")"
export MPLBACKEND=Agg

PY=python3
titulo() { printf '\n===============================================================================\n %s\n===============================================================================\n' "$1"; }

# -----------------------------------------------------------------------------
#  DEPENDENCIAS. Comprobadas aquí a propósito: al cambiar de máquina (laptop ->
#  PC) faltaba numpy y lo que salía era un traceback de importación en mitad de
#  un título, que no dice qué hacer. Esto lo dice.
# -----------------------------------------------------------------------------
faltan=""
for _m in numpy matplotlib; do
    $PY -c "import $_m" 2>/dev/null || faltan="$faltan $_m"
done
if [ -n "$faltan" ]; then
    echo
    echo "  FALTAN DEPENDENCIAS DE PYTHON:$faltan"
    echo
    echo "  Arch / EndeavourOS / Manjaro:"
    echo "    sudo pacman -S python-numpy python-matplotlib python-plotly"
    echo "  Debian / Ubuntu / Mint:"
    echo "    sudo apt install python3-numpy python3-matplotlib python3-plotly"
    echo "  Fedora:"
    echo "    sudo dnf install python3-numpy python3-matplotlib python3-plotly"
    echo
    echo "  (plotly solo hace falta para 'sh correr.sh acople')"
    echo
    exit 3
fi

pruebas() { $PY pruebas.py; }

tabla() {
    titulo "VEREDICTOS"
    $PY rutas.py --callado
    $PY rutas.py --detalle --callado
}

figuras() {
    titulo "FIGURAS"
    $PY rutas.py --figuras --callado
    $PY envolvimiento_script.py --figuras >/dev/null
    echo "    envolvimiento_radio_critico.png"
    echo "    envolvimiento_barrera.png"
    echo "    envolvimiento_G_de_D.png"
    # Las de dendrímero, polímero y teóricos NO las hacía 'figuras': cada una
    # salía solo con su propio comando. Como 'web' se limita a COPIAR los PNG
    # que encuentre, la web acababa mezclando figuras generadas en momentos
    # distintos y con versiones distintas del código, y se contradecían entre
    # sí. Bug detectado por Jhovan el 2026-08-12. Ahora 'figuras' las genera
    # TODAS y 'web' llama a 'figuras' antes de construir.
    # DESACTIVADO 2026-08-24 (decisión de Jhovan): el simulador es exclusivo
    # para liposomas; 'rutas.py --dendrimero/--polimero/--teoricos' ya no
    # existen. Registro histórico, no se borra.
    # $PY rutas.py --dendrimero >/dev/null
    # $PY rutas.py --polimero   >/dev/null
    # $PY rutas.py --teoricos   >/dev/null
    # echo "    + dendrímero, polímero y teóricos"
    echo
    echo "  en $(pwd)"
}

# DESACTIVADAS el 2026-08-24 (decisión de Jhovan): el simulador es exclusivo
# para liposomas. 'rutas.py --dendrimero/--teoricos/--polimero' ya no
# existen en rutas.py. Registro histórico, no se borran.
# dendrimero() {
#     titulo "DENDRÍMERO PAMAM · tarea G.1a"
#     $PY rutas.py --dendrimero
# }
#
# teoricos() {
#     titulo "TRANSPORTADORES TEÓRICOS DEL SIMULADOR DE ACOPLE · datos sintéticos"
#     $PY rutas.py --teoricos
# }
#
# polimero() {
#     titulo "POLÍMERO MACIZO Y MICELA · tarea G.1b"
#     $PY rutas.py --polimero
# }

web() {
    # REGENERA LAS FIGURAS ANTES DE COPIARLAS. construir_web.py ejecuta el
    # simulador para los NÚMEROS, pero las figuras son PNG en disco y solo los
    # copia: si estaban viejos, la web mostraba figuras que se contradecían con
    # sus propias tablas. Con esto no puede volver a pasar.
    figuras
    titulo "WEB ESTÁTICA · se genera ejecutando el simulador"
    $PY construir_web.py
}

notas() { $PY rutas.py --notas; }

farmaco() {
    # rdkit NO se comprueba en el bloque general de dependencias (arriba):
    # solo lo usa tamano_farmaco.py (línea 84), 'farmaco' es un subcomando
    # puntual, no parte de 'pruebas' ni de 'web' (la corrida normal del
    # simulador nunca lo importa). Bloquear todo el script por una
    # dependencia que el 90% de las corridas no toca sería peor que dar el
    # aviso aquí, igual que ya se hace con plotly para 'acople'.
    if ! $PY -c "import rdkit" 2>/dev/null; then
        echo
        echo "  FALTA rdkit (lo necesita tamano_farmaco.py)"
        echo
        echo "  Arch / EndeavourOS / Manjaro:"
        echo "    sudo pacman -S rdkit"
        echo "  Debian / Ubuntu / Mint:"
        echo "    sudo apt install python3-rdkit"
        echo "  Fedora:"
        echo "    sudo dnf install python3-rdkit"
        echo "  (NO uses 'pip install rdkit': en Arch el Python del sistema"
        echo "   es externally-managed (PEP 668) y pip falla)"
        echo
        exit 3
    fi
    $PY tamano_farmaco.py
}

colab() {
    $PY construir_colab.py
    $PY verifica_equivalencia.py
}

cribado() {
    titulo "CRIBADO DE ADHESIÓN · ojo, bug P5 (ver PENDIENTES.md)"
    $PY ../cribado_bhe_tres_liposomas.py
}

acople() {
    titulo "ACOPLE FÁRMACO-LIPOSOMA"
    if ! $PY -c "import plotly" 2>/dev/null; then
        echo "  falta 'plotly':  pip install plotly"
        return 0
    fi
    # Acaba abriendo una animación 3D en el navegador. Sin navegador esa última
    # línea revienta aunque el cálculo haya ido bien, y set -e tumbaría todo.
    if $PY ../acople/acople_fingolimod_liposoma.py; then
        return 0
    fi
    echo
    echo "  (si el informe salió completo, falló solo la animación 3D:"
    echo "   necesita navegador. Los números son válidos)"
    return 0
}

case "${1:-todo}" in
    pruebas) pruebas ;;
    tabla)   tabla ;;
    figuras) figuras ;;
    # dendrimero) dendrimero ;;   # DESACTIVADO 2026-08-24, ver nota arriba
    # teoricos) teoricos ;;       # DESACTIVADO 2026-08-24, ver nota arriba
    # polimero) polimero ;;       # DESACTIVADO 2026-08-24, ver nota arriba
    web)     web ;;
    notas)   notas ;;
    farmaco) farmaco ;;
    colab)   colab ;;
    cribado) cribado ;;
    acople)  acople ;;
    todo)
        pruebas
        tabla
        figuras
        echo
        echo "  sh correr.sh notas       lectura de los resultados y qué cambió"
        # echo "  sh correr.sh dendrimero  ventana geométrica del PAMAM (G.1a)"  # DESACTIVADO 2026-08-24
        echo "  sh correr.sh web         genera la web estática en ../web"
        echo "  sh correr.sh cribado  ·  sh correr.sh acople   (módulos de fuera)"
        echo
        ;;
    *)
        echo "No conozco la opción '$1'."
        echo "Usa: todo | pruebas | tabla | figuras | web | notas | farmaco | colab | cribado | acople"
        exit 2
        ;;
esac
