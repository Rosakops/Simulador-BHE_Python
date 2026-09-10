#!/usr/bin/env python3
# =============================================================================
#  grafica_dispersion_dataset_50.py · última día, 2026-09-09
# =============================================================================
#  Figura de dispersión (Diámetro vs Potencial zeta) de los 50 liposomas
#  sintéticos de dataset_50_liposomas.py, un panel por ruta (A-D). Sigue la
#  misma regla de oro que construir_dataset_50_web.py: este script NO
#  recalcula nada nuevo, solo ejecuta generar()/evaluar() de los módulos ya
#  validados y dibuja lo que obtiene.
#
#  Por qué un panel por ruta y no un solo scatter: cada liposoma tiene un
#  veredicto DISTINTO por ruta (A/B/C/D), así que no hay una sola categoría
#  de color que sirva para las 4 a la vez. La banda de la compuerta de
#  caveola (60-80 nm) solo se dibuja en el panel A porque g_caveola es la
#  única compuerta de las 4 rutas que depende de ese rango (ver RUTAS en
#  rutas.py) — ponerla en B/C/D sería engañoso, esas rutas no la usan.
#
#  Nota de lectura: en este dataset A, C y D dan el mismo veredicto en los
#  50 liposomas porque las tres comparten la compuerta "Difusión en espacio
#  extracelular (transportador)" (g_difusion_ecs(d, "transportador")), y es
#  la única que llega a FALLA en estos 50 diseños. No es un error del
#  script ni casualidad del muestreo: es la misma función evaluada sobre el
#  mismo diseño en las tres rutas. Ver rutas.py, definición de RUTAS.
#
#  SALE: web/img/dataset_50/dispersion_diametro_zeta_4rutas.png
# =============================================================================

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import rutas as R
from dataset_50_liposomas import generar

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
IMG_DIR = RAIZ / "web" / "img" / "dataset_50"

# Mismo código de color que usan las otras figuras del proyecto
# (_color_ruta en construir_dataset_50_web.py), con un tercer color por si
# algún diseño llegara a NO EXCLUIDA (no ocurre en este dataset, pero el
# código no debe asumirlo).
COLOR_VEREDICTO = {
    "EXCLUIDA": "#e02424",
    "NO EVALUABLE": "#2255cc",
    "NO EXCLUIDA": "#2e7d32",
}

XLIM = (20, 720)
YLIM = (-31, 9)


def _generar_figura(disenos, ruta_pathfile):
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    diametro = [d.diametro_nm for d in disenos]
    zeta = [d.zeta_mV for d in disenos]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
    fig.patch.set_facecolor("white")

    handles_por_etiqueta = {}
    for ax, (nombre_ruta, compuertas) in zip(axes.flat, R.RUTAS.items()):
        letra = nombre_ruta[0]
        veredictos_ruta = [R.evaluar_ruta(d, compuertas)[0] for d in disenos]

        puntos = {}
        for x, y, v in zip(diametro, zeta, veredictos_ruta):
            puntos.setdefault(v, ([], []))
            puntos[v][0].append(x)
            puntos[v][1].append(y)

        if letra == "A":
            ax.axvspan(60, 80, color="#f0a500", alpha=0.12, zorder=0)
            ax.text(70, YLIM[1] - 0.8, "caveola\n60–80 nm", ha="center",
                    va="top", fontsize=6.5, color="#a67c00")

        ax.axhline(0, color="#888888", lw=0.8, zorder=1)

        for v, (xs, ys) in puntos.items():
            color = COLOR_VEREDICTO.get(v, "#333333")
            h = ax.scatter(xs, ys, color=color, s=32, label=v, zorder=3)
            handles_por_etiqueta[v] = h

        n_excl = sum(1 for v in veredictos_ruta if v == "EXCLUIDA")
        n_noeval = sum(1 for v in veredictos_ruta if v == "NO EVALUABLE")
        n_noexcl = sum(1 for v in veredictos_ruta if v == "NO EXCLUIDA")
        partes_conteo = [f"excl. {n_excl}", f"no eval. {n_noeval}"]
        if n_noexcl:
            partes_conteo.append(f"no excl. {n_noexcl}")
        ax.set_title(f"{nombre_ruta}  ({' · '.join(partes_conteo)})", fontsize=10.5)

        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        ax.grid(True, alpha=0.3)

    fig.text(0.5, 0.035, "Diámetro (nm)", ha="center", fontsize=12, color="#c0392b")
    fig.text(0.02, 0.5, "Potencial zeta (mV)", va="center", rotation="vertical",
              fontsize=12, color="#c0392b")

    fig.suptitle("Diámetro vs Potencial zeta por ruta de entrada a la BHE\n"
                 f"({len(disenos)} liposomas sintéticos)",
                 fontsize=14, fontweight="bold", y=0.995)

    # Orden fijo para que la leyenda no cambie de orden entre corridas.
    orden = ["EXCLUIDA", "NO EVALUABLE", "NO EXCLUIDA"]
    etiquetas = [e for e in orden if e in handles_por_etiqueta]
    handles = [handles_por_etiqueta[e] for e in etiquetas]
    fig.legend(handles=handles, labels=etiquetas, loc="lower center", ncol=len(etiquetas),
               bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=10)

    fig.tight_layout(rect=[0.035, 0.06, 1, 0.96])
    fig.savefig(ruta_pathfile, dpi=160, bbox_inches="tight")
    plt.close(fig)


def generar_grafica(disenos=None):
    """Genera la figura y devuelve la ruta del PNG. Reutilizable desde
    construir_dataset_50_web.py u otro script."""
    disenos = disenos if disenos is not None else generar()
    ruta = IMG_DIR / "dispersion_diametro_zeta_4rutas.png"
    _generar_figura(disenos, ruta)
    return ruta


if __name__ == "__main__":
    ruta = generar_grafica()
    print(f"  {ruta}")
