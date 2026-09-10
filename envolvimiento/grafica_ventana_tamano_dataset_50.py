#!/usr/bin/env python3
# =============================================================================
#  grafica_ventana_tamano_dataset_50.py · última día, 2026-09-09
# =============================================================================
#  Ventanas de tamaño (compuertas que dependen del diámetro) con los 50
#  liposomas sintéticos superpuestos como puntos, no como líneas verticales.
#  Reutiliza las bandas y umbrales YA calculados por figuras() en rutas.py
#  (_ventana_fabricable, _umbrales_modelo, D_SIGILOSO_*) para que los rangos
#  nunca se desincronicen del motor real si esos umbrales cambian.
#
#  Por qué puntos y no las 50 líneas de figuras(): con 50 diseños las líneas
#  verticales con etiqueta rotada se amontonan y dejan de leerse (ver
#  discusión con Jhovan, 2026-09-09). Un punto por diseño y por compuerta,
#  con jitter vertical dentro de su fila y coloreado PASA/FALLA/DESCONOCIDA
#  para ESA compuerta específica, es la misma información sin el amontonado.
#
#  SALE: web/img/dataset_50/ventana_tamano_dataset_50.png
# =============================================================================

from pathlib import Path
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import rutas as R
from dataset_50_liposomas import generar

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
IMG_DIR = RAIZ / "web" / "img" / "dataset_50"

VERDE, ROJO, GRIS = "#2e7d32", "#c62828", "#9e9e9e"
COLOR_ESTADO = {"PASA": VERDE, "FALLA": ROJO, "DESCONOCIDA": GRIS}


def generar_grafica(disenos=None):
    disenos = disenos if disenos is not None else generar()
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    etq_fab, lo_fab, hi_fab, fab_indef = R._ventana_fabricable(disenos)
    u = R._umbrales_modelo()

    # (nombre visible, función que evalúa ESA compuerta para un diseño,
    #  límite inferior, límite superior de la banda permitida)
    filas = [
        (etq_fab, lambda d: R.g_transportador_fabricable(d), lo_fab, hi_fab),
        ("Tamiz del glicocálix", lambda d: R.g_glicocalix_tamiz(d), 1.0, u["glicocalix"]),
        ("Envolvimiento de membrana", lambda d: R.g_envolvimiento(d), u["envolvimiento"], 2000.0),
        ("Compuerta de caveola", lambda d: R.g_caveola(d), 1.0, 80.0),
        ("Difusión extracelular\n(transportador)", lambda d: R.g_difusion_ecs(d, "transportador"), 1.0, R.D_SIGILOSO_PASA_nm),
        ("Captación por macrófago", lambda d: R.g_captacion_fagocitica(d), 550.0, 2000.0),
    ]

    fig, ax = plt.subplots(figsize=(12.5, 6.4))
    rng = random.Random(42)  # mismo jitter reproducible en cada corrida

    for i, (nombre, fn, lo, hi) in enumerate(filas):
        y = len(filas) - 1 - i
        ax.barh(y, hi - lo, left=lo, height=0.62, color=VERDE, alpha=0.18,
                edgecolor="black", linewidth=0.5, zorder=0)

        for d in disenos:
            r = fn(d)
            color = COLOR_ESTADO.get(r.estado, "#333333")
            jitter = rng.uniform(-0.26, 0.26)
            ax.scatter(d.diametro_nm, y + jitter, color=color, s=16,
                       alpha=0.75, zorder=3, edgecolors="none")

    ax.set_yticks([len(filas) - 1 - i for i in range(len(filas))])
    ax.set_yticklabels([nom for nom, *_ in filas], fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(1, 2000)
    ax.set_ylim(-0.7, len(filas) - 0.3)
    ax.set_xlabel("Diámetro (nm)")
    ax.set_title("Ventanas de tamaño · 50 liposomas sintéticos\n"
                 "banda = rango permitido por la compuerta · punto = un liposoma, "
                 "coloreado por su resultado en ESA compuerta", fontsize=11)
    ax.grid(axis="x", alpha=0.3, which="both")

    ax.legend(handles=[
        plt.matplotlib.patches.Patch(facecolor=VERDE, alpha=0.75, label="PASA"),
        plt.matplotlib.patches.Patch(facecolor=ROJO, alpha=0.75, label="FALLA"),
        plt.matplotlib.patches.Patch(facecolor=GRIS, alpha=0.75, label="DESCONOCIDA"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=9, frameon=False)

    fig.tight_layout(rect=[0, 0.08, 1, 1])
    ruta = IMG_DIR / "ventana_tamano_dataset_50.png"
    fig.savefig(ruta, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return ruta


if __name__ == "__main__":
    ruta = generar_grafica()
    print(f"  {ruta}")
