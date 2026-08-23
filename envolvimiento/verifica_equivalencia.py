#!/usr/bin/env python3
# =============================================================================
#  VERIFICACIÓN CRUZADA: la versión de Colab y la de script deben dar
#  EXACTAMENTE los mismos números.
#
#  Es la salvaguarda contra el error clásico de mantener dos copias: que se
#  arregle un bug en una y no en la otra. Si este archivo falla, la versión de
#  Colab está desactualizada -> vuelve a ejecutar construir_colab.py
#
#  Uso:  python3 verifica_equivalencia.py
# =============================================================================

import importlib.util
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
#  OBLIGATORIO Y ANTES DE NADA: forzar el backend no interactivo.
#
#  envolvimiento_colab.py es un script de Colab: se ejecuta ENTERO al importarlo
#  y termina en plt.show(). Este archivo lo importa para comparar los números.
#  En una máquina con escritorio (Wayland/Qt), ese plt.show() abre una ventana y
#  DEJA EL PROCESO BLOQUEADO esperando a que se cierre a mano; la verificación
#  parece colgada y nunca imprime su resultado.
#
#  Con MPLBACKEND=Agg no hay ventana y plt.show() no hace nada. Se pone por
#  variable de entorno porque tiene que estar puesto ANTES de que colab importe
#  pyplot, y aquí todavía no lo hemos importado.
#
#  Se fuerza SIN excepción (no setdefault): esta verificación es numérica y no
#  tiene ningún motivo para dibujar. El plt.show() de la plantilla de Colab es
#  correcto allí, donde la figura sale en línea; el problema es solo importarlo
#  desde un escritorio.
# ---------------------------------------------------------------------------
os.environ["MPLBACKEND"] = "Agg"

import numpy as np

AQUI = Path(__file__).resolve().parent


def cargar(nombre, ruta):
    spec = importlib.util.spec_from_file_location(nombre, ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[nombre] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    ruta_colab = AQUI / "envolvimiento_colab.py"
    if not ruta_colab.exists():
        print("FALTA envolvimiento_colab.py: ejecuta primero construir_colab.py")
        return 1

    core = cargar("_core", AQUI / "envolvimiento_core.py")

    # La versión de Colab se ejecuta al importarla: imprime su informe y ADEMÁS
    # intenta dibujar sus figuras con plt.show(). Aquí solo queremos los
    # números, así que se silencian las tres salidas que produce:
    #
    #   · stdout  -> su informe completo
    #   · stderr  -> el UserWarning "FigureCanvasAgg is non-interactive, and
    #                thus cannot be shown", que es ESPERADO y no indica ningún
    #                problema: es matplotlib avisando de que con backend Agg no
    #                hay ventana que mostrar. Es justo lo que queremos.
    #   · warnings -> por si acaso escapa alguno más
    #
    # Se silencia por higiene: ese aviso aparecía en medio de la verificación y
    # daba la impresión de que algo había fallado cuando todo iba bien.
    import contextlib
    import io
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            colab = cargar("_colab", ruta_colab)

    casos = [(R, z, p, h, k, s)
             for R in (8.0, 20.0, 45.0, 75.0)
             for z in (-30.0, -5.0, 0.0, 5.0, 30.0)
             for p in (0.0, 5.0)
             for h in (3.0e-21, 6.5e-21)
             for k in (15.0, 50.0)
             for s in (0.003, 0.3)]

    fallos = 0
    for c in casos:
        R, z, p, h, k, s = c
        a = core.clasificar(R, z, p, h, k, s)
        b = colab.clasificar(R, z, p, h, k, s)
        for clave, va in a.items():
            vb = b[clave]
            if isinstance(va, bool):
                igual = va == vb
            elif isinstance(va, (int, float)):
                igual = (np.isinf(va) and np.isinf(vb)) or np.isclose(
                    va, vb, rtol=0, atol=0, equal_nan=True)
            else:
                igual = va == vb
            if not igual:
                fallos += 1
                if fallos <= 10:
                    print(f"  DIFIERE {clave} en {c}: core={va!r} colab={vb!r}")

    print("=" * 70)
    print(f" Casos comparados: {len(casos)}  |  campos por caso: {len(a)}")
    if fallos == 0:
        print(" RESULTADO: idénticos, bit a bit. Las dos versiones son la misma física.")
    else:
        print(f" RESULTADO: {fallos} diferencias. Vuelve a ejecutar construir_colab.py")
    print("=" * 70)
    return 0 if fallos == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
