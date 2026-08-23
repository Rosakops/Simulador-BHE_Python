#!/usr/bin/env python3
# =============================================================================
#  DATASET DE 50 LIPOSOMAS SINTÉTICOS: último día, 2026-08-18
# =============================================================================
#  Genera 50 diseños de liposoma variando aleatoriamente los TRES parámetros
#  físico-químicos que el simulador usa realmente para esta clase (diametro_nm,
#  zeta_mV, peg_nm: los únicos campos de `Diseno` que dependen del diseño del
#  liposoma; farmaco_diametro_nm es propiedad del fármaco y clase="liposoma" es
#  fija). Los corre por las 4 rutas (A/B/C/D) con `evaluar()`, exactamente el
#  mismo motor que usa rutas.py para el CATALOGO.
#
#  RANGOS: los de los diseños reales ya validados en CATALOGO_REAL + teóricos
#  (Mao 2014, Gong 2022, Chow 2025, Muselman 2026 + los 3 teóricos):
#    diametro_nm  30  - 700  nm
#    zeta_mV     -28  - +7   mV
#    peg_nm        0  - 5    nm
#  Decisión de Jhovan (2026-08-18): no ampliar con fuentes nuevas por tiempo.
#
#  CARGA ÚTIL (G.2): queda FUERA del bucle. No es un campo de `Diseno`, es una
#  compuerta aparte (g_carga_util) que hoy devuelve DESCONOCIDA para todo el
#  CATALOGO por decisión del 2026-08-13, vigente. Los 50 diseños heredan esa
#  misma DESCONOCIDA en C y D (las únicas rutas que la usan), sin inventar ni
#  fijar un valor nuevo bajo presión de tiempo.
#
#  NOMBRES: "Liposoma 1" .. "Liposoma 50" (decisión de Jhovan, 2026-08-18).
#  sintetico=True en los 50 → nunca citables como predicción real (mismo
#  criterio que CATALOGO_TEORICO).
#
#  SEMILLA FIJA (42) para que el dataset sea reproducible.
# =============================================================================

import csv
import random

from rutas import Diseno, evaluar, RUTAS

SEMILLA = 42
N = 50

D_MIN, D_MAX = 30.0, 700.0
Z_MIN, Z_MAX = -28.0, 7.0
P_MIN, P_MAX = 0.0, 5.0


def generar(n=N, semilla=SEMILLA):
    rng = random.Random(semilla)
    disenos = []
    for i in range(1, n + 1):
        d = Diseno(
            nombre=f"Liposoma {i}",
            diametro_nm=round(rng.uniform(D_MIN, D_MAX), 2),
            zeta_mV=round(rng.uniform(Z_MIN, Z_MAX), 2),
            peg_nm=round(rng.uniform(P_MIN, P_MAX), 2),
            clase="liposoma",
            sintetico=True,
            nota="Dataset 50 aleatorio 2026-08-18, semilla 42, rango de CATALOGO real",
        )
        disenos.append(d)
    return disenos


def correr(disenos):
    filas = []
    for d in disenos:
        veredictos = evaluar(d)
        fila = {
            "nombre": d.nombre,
            "diametro_nm": d.diametro_nm,
            "zeta_mV": d.zeta_mV,
            "peg_nm": d.peg_nm,
        }
        for ruta_nombre in RUTAS:
            v, _ = veredictos[ruta_nombre]
            fila[ruta_nombre] = v
        filas.append(fila)
    return filas


def escribir_csv(filas, ruta="dataset_50_liposomas.csv"):
    campos = ["nombre", "diametro_nm", "zeta_mV", "peg_nm"] + list(RUTAS.keys())
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)
    return ruta


def resumen(filas):
    print(f"\n  {len(filas)} liposomas sintéticos, semilla {SEMILLA}\n")
    print(f"  {'nombre':<14}{'Ø nm':>8}{'ζ mV':>8}{'PEG nm':>8}   " +
          "  ".join(f"{r:<28}" for r in RUTAS))
    for fila in filas:
        print(f"  {fila['nombre']:<14}{fila['diametro_nm']:>8.1f}"
              f"{fila['zeta_mV']:>8.2f}{fila['peg_nm']:>8.2f}   " +
              "  ".join(f"{fila[r]:<28}" for r in RUTAS))

    print("\n  Conteo por ruta:")
    for r in RUTAS:
        cont = {}
        for fila in filas:
            cont[fila[r]] = cont.get(fila[r], 0) + 1
        resumen_str = ", ".join(f"{k}: {v}" for k, v in sorted(cont.items()))
        print(f"    {r:<32} {resumen_str}")


if __name__ == "__main__":
    disenos = generar()
    filas = correr(disenos)
    ruta = escribir_csv(filas)
    resumen(filas)
    print(f"\n  CSV escrito en {ruta}")
