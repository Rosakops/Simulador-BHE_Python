#!/usr/bin/env python3
# =============================================================================
#  CORREDOR DE PRUEBAS: una línea por suite
# =============================================================================
#  Corre las cuatro suites del proyecto y da UNA línea por cada una.
#  Si alguna falla, y SOLO entonces, despliega su salida completa: cuando todo
#  va bien el detalle no aporta nada, y cuando algo se rompe lo quieres entero.
#
#  Uso:  python3 pruebas.py          compacto
#        python3 pruebas.py --todo   despliega el detalle aunque pase
# =============================================================================

import contextlib
import importlib.util
import io
import os
import re
import sys
from pathlib import Path

os.environ["MPLBACKEND"] = "Agg"   # antes de que nadie importe pyplot

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))


def _cargar(nombre, ruta):
    spec = importlib.util.spec_from_file_location(nombre, ruta)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[nombre] = mod
    spec.loader.exec_module(mod)
    return mod


def _marcador(salida):
    """Marcador corto de cada suite: 'N/M', o el nº de casos si no hay conteo."""
    for linea in salida.splitlines():
        if "RESULTADO" in linea:
            m = re.search(r"(\d+/\d+)", linea)
            if m:
                return m.group(1)
            break
    # la verificación de equivalencia no cuenta pruebas, cuenta casos
    m = re.search(r"Casos comparados:\s*(\d+)", salida)
    if m:
        return f"{m.group(1)} casos"
    return "—"


def correr(detalle=False):
    import rutas
    import glicocalix
    import envolvimiento_core as E

    def equivalencia():
        m = _cargar("_ve", AQUI / "verifica_equivalencia.py")
        return m.main() == 0

    suites = [
        ("rutas.py",            lambda: rutas.validar_contra_experimentos()),
        ("glicocalix.py",       lambda: glicocalix.test_glicocalix()),
        ("envolvimiento_core",  lambda: E.test_limites()),
        ("equivalencia colab",  equivalencia),
    ]

    filas, fallos = [], []
    for nombre, fn in suites:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            ok = bool(fn())
        salida = buf.getvalue()
        filas.append((nombre, ok, _marcador(salida)))
        if not ok:
            fallos.append((nombre, salida))

    print()
    print("  PRUEBAS")
    for nombre, ok, marca in filas:
        punto = "." * max(2, 26 - len(nombre))
        print(f"    {nombre} {punto} {marca:>16s}   {'OK' if ok else 'FALLA'}")

    if detalle and not fallos:
        for nombre, fn in suites:
            print(f"\n  --- detalle de {nombre} ---")
            fn()

    for nombre, salida in fallos:
        print(f"\n  ================ DETALLE DEL FALLO · {nombre} ================")
        print(salida)

    print()
    return not fallos


if __name__ == "__main__":
    sys.exit(0 if correr("--todo" in sys.argv[1:]) else 1)
