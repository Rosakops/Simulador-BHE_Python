#!/usr/bin/env python3
# =============================================================================
#  grafica_matriz_dataset_50.py · última día, 2026-09-09
# =============================================================================
#  Matriz de veredictos (diseño x ruta) de los 50 liposomas sintéticos, como
#  heatmap compacto y TRANSPUESTO (4 filas = rutas, 50 columnas = liposomas)
#  en vez de la tabla de 50 filas x 4 columnas con texto en cada celda que
#  usa figuras() en rutas.py -- esa se vuelve una tira muy larga con 50
#  diseños. El color solo (sin texto por celda) es suficiente con este
#  número de columnas; el detalle por diseño ya está en dataset_50.html.
#
#  SALE: web/img/dataset_50/matriz_veredictos_dataset_50.png
# =============================================================================

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import rutas as R
from dataset_50_liposomas import generar

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
IMG_DIR = RAIZ / "web" / "img" / "dataset_50"

VERDE, ROJO, GRIS = "#2e7d32", "#c62828", "#9e9e9e"
COD = {"EXCLUIDA": 0, "NO EVALUABLE": 1, "NO EXCLUIDA": 2}
COLORES = [ROJO, GRIS, VERDE]  # mismo orden que COD


def generar_grafica(disenos=None):
    disenos = disenos if disenos is not None else generar()
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    nombres_rutas = list(R.RUTAS.keys())
    matriz = np.zeros((len(nombres_rutas), len(disenos)), dtype=int)
    for j, d in enumerate(disenos):
        veredictos = R.evaluar(d)
        for i, nombre_ruta in enumerate(nombres_rutas):
            v, _ = veredictos[nombre_ruta]
            matriz[i, j] = COD[v]

    cmap = matplotlib.colors.ListedColormap(COLORES)
    fig, ax = plt.subplots(figsize=(14, 3.6))
    ax.imshow(matriz, cmap=cmap, vmin=0, vmax=2, aspect="auto")

    ax.set_yticks(range(len(nombres_rutas)))
    ax.set_yticklabels([n.split(" (")[0] for n in nombres_rutas], fontsize=10)
    ax.set_xticks(range(len(disenos)))
    ax.set_xticklabels([str(i + 1) for i in range(len(disenos))], fontsize=6, rotation=90)
    ax.set_xlabel("Liposoma #", fontsize=10)
    ax.set_title("Matriz de veredictos · 50 liposomas sintéticos × 4 rutas de entrada a la BHE",
                 fontsize=12, fontweight="bold")

    # líneas finas entre celdas
    ax.set_xticks(np.arange(-0.5, len(disenos), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(nombres_rutas), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    ax.legend(handles=[
        plt.matplotlib.patches.Patch(facecolor=ROJO, label="EXCLUIDA"),
        plt.matplotlib.patches.Patch(facecolor=GRIS, label="NO EVALUABLE"),
        plt.matplotlib.patches.Patch(facecolor=VERDE, label="NO EXCLUIDA"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=3, fontsize=9, frameon=False)

    fig.tight_layout(rect=[0, 0.12, 1, 1])
    ruta = IMG_DIR / "matriz_veredictos_dataset_50.png"
    fig.savefig(ruta, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return ruta


if __name__ == "__main__":
    ruta = generar_grafica()
    print(f"  {ruta}")
