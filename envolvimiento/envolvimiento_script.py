#!/usr/bin/env python3
# =============================================================================
#  ENVOLVIMIENTO DE MEMBRANA: versión SCRIPT (Python local)
#  Proyecto BHE / SENACYT. Tareas 2.4, 2.5 y 2.6.
#
#  Uso:
#      python3 envolvimiento_script.py                 # informe completo
#      python3 envolvimiento_script.py --solo-tests    # solo validación
#      python3 envolvimiento_script.py --figuras       # además guarda PNG
#
#  La física vive en envolvimiento_core.py. Este archivo solo presenta.
#  La versión de Colab se GENERA desde estos dos con construir_colab.py, para
#  que no puedan divergir.
# =============================================================================

import argparse
import numpy as np

import envolvimiento_core as E


# ---- fichas de los tres diseños del estudio --------------------------------
# AVISO: los zeta son valores de DISEÑO asignados en un rango, NO medidas.
# Esto es un cribado del espacio de diseño, no una predicción sobre
# formulaciones concretas.
LIPOSOMAS = [
    dict(nombre="Convencional", R_nm=20.0, zeta_mV=+5.0, peg_nm=0.0),
    dict(nombre="Furtivo/PEG",  R_nm=15.5, zeta_mV=+2.0, peg_nm=5.0),
    dict(nombre="Catiónico",    R_nm=17.5, zeta_mV=+6.7, peg_nm=0.0),
]

HAMAKER_REF = 4.5e-21   # valor central del barrido, para las tablas de ejemplo


def cabecera():
    print("=" * 79)
    print(" ENVOLVIMIENTO DE MEMBRANA: ¿puede el endotelio envolver el nanotransportador?")
    print("=" * 79)
    print(" Fuente: Deserno M., Phys. Rev. E 69, 031903 (2004).")
    print(f" Endotelio: zeta = {E.ZETA_BHE_mV} mV (Santa-Maria 2019, Fig. 4A). T = 37 °C.")
    print(f" Barridos de IGNORANCIA  kappa {E.KAPPA_kT} kT | sigma {E.SIGMA_mNm} mN/m"
          f" | Hamaker {[f'{h:.1e}' for h in E.HAMAKER_J]} J")
    print()
    print(" LO QUE ESTE MÓDULO AFIRMA:")
    print("   Un diseño EXCLUIDO no puede iniciar transcitosis por envolvimiento de")
    print("   membrana. Esa afirmación es fuerte.")
    print("   Un diseño NO EXCLUIDO es un CANDIDATO, no una predicción de éxito.")
    print("   Esa afirmación es débil, y así hay que escribirla.")
    print("=" * 79)


def tabla_compuertas(kappa_kT=25.0, sigma_mNm=0.03, hamaker_J=HAMAKER_REF):
    print(f"\n### COMPUERTAS · kappa={kappa_kT:.0f} kT, sigma={sigma_mNm} mN/m, "
          f"A_H={hamaker_J:.1e} J\n")
    print(f"{'diseño':14s} {'w(uN/m)':>9} {'R_min':>7} {'pozo':>9} {'barr.ent':>9}"
          f" {'barr.env':>9}  G1 G2 G3 G4  veredicto")
    print("-" * 79)
    for lip in LIPOSOMAS:
        r = E.clasificar(lip["R_nm"], lip["zeta_mV"], lip["peg_nm"],
                         hamaker_J, kappa_kT, sigma_mNm)
        marca = lambda b: " ✓" if b else " ✗"
        veredicto = "candidato" if r["NO_EXCLUIDO"] else "EXCLUIDO"
        print(f"{lip['nombre']:14s} {r['w_uNm']:9.0f} {r['R_min_nm']:6.1f}n "
              f"{r['pozo_dlvo_kT']:8.1f}k {r['barrera_entrada_kT']:8.2f}k "
              f"{r['barrera_envolv_kT']:8.1f}k "
              f"{marca(r['G1_radio_critico'])}{marca(r['G2_w_sobre_sigma'])}"
              f"{marca(r['G3_envolv_completo'])}{marca(r['G4_caveola'])}  {veredicto}")
    print("-" * 79)
    print(" G1 radio crítico R > sqrt(2k/w)   G2 w/sigma > 1.37   "
          "G3 envolvimiento completo   G4 caveola <= 80 nm")
    print(" barr.ent = barrera DLVO de entrada (cruzarla para llegar al contacto)")
    print(" barr.env = barrera de la transición de envolvimiento (compuerta CINÉTICA)")


def tabla_sensibilidad():
    print("\n### SENSIBILIDAD: por qué sigma es el parámetro que manda\n")
    lip = LIPOSOMAS[0]
    print(f" Diseño de referencia: {lip['nombre']} "
          f"(R={lip['R_nm']} nm, zeta={lip['zeta_mV']:+.1f} mV, PEG={lip['peg_nm']} nm)\n")
    print(f"{'sigma(mN/m)':>12} {'kappa(kT)':>10} {'sigma_tilde':>12} "
          f"{'R_min(nm)':>10} {'barrera env.(kT)':>18}")
    print("-" * 68)
    for s in E.SIGMA_mNm:
        for k in E.KAPPA_kT:
            r = E.clasificar(lip["R_nm"], lip["zeta_mV"], lip["peg_nm"],
                             HAMAKER_REF, k, s)
            print(f"{s:12.3f} {k:10.0f} {r['sigma_tilde']:12.4f} "
                  f"{r['R_min_nm']:10.1f} {r['barrera_envolv_kT']:18.1f}")
    print("-" * 68)
    print(" El barrido completo de kappa mueve R_min un factor ~1.8.")
    print(" El barrido completo de sigma mueve la barrera un factor ~50.")
    print(" => el análisis de sensibilidad de la Fase 4 debe centrarse en SIGMA.")


def lectura():
    print("\n" + "=" * 79)
    print(" LECTURA: qué discrimina y qué no")
    print("=" * 79)
    print(" 1. Las compuertas TERMODINÁMICAS no discriminan en este rango. La adhesión")
    print("    de van der Waals en el contacto es del orden de miles de uN/m, así que")
    print("    R_min sale ~7 nm y todo diseño por encima de eso la supera. Reportar")
    print("    esto como resultado: la barrera al envolvimiento NO es termodinámica.")
    print()
    print(" 2. Lo que sí discrimina es CINÉTICO:")
    print("      · la barrera DLVO de entrada, que solo aparece con PEG;")
    print("      · la barrera de la transición de envolvimiento, gobernada por sigma.")
    print()
    print(" 3. sigma no está medida para endotelio cerebral. Es decir: el parámetro")
    print("    que decide el resultado es justamente el que nadie ha medido. Eso hay")
    print("    que decirlo en el resumen, no esconderlo en las limitaciones.")
    print("=" * 79)
    print(" LÍMITES QUE HAY QUE DECLARAR EN EL ARTÍCULO:")
    print(" L1. Derjaguin supone superficies RÍGIDAS; Deserno modela la DEFORMACIÓN de")
    print("     la membrana. Encadenarlos es un híbrido, no una derivación limpia.")
    print(" L2. Deserno supone coloide RÍGIDO. Un liposoma es deformable, lo que abarata")
    print("     el envolvimiento. El modelo es CONSERVADOR en ese sentido.")
    print(" L3. El glicocálix (0.2-5 um, Walter 2021) no está en el modelo. Es de 6 a 160")
    print("     veces más grueso que estos liposomas y plausiblemente es el paso")
    print("     limitante real.")
    print(" L4. La barrera de envolvimiento es una interpolación de escalados de Deserno")
    print("     (Figs. 3 y 9), no una fórmula exacta. Usar como ORDEN DE MAGNITUD.")
    print(" L5. Los zeta son valores de DISEÑO, no medidas.")
    print(" L6. kappa y sigma entran como BARRIDO DE IGNORANCIA: van en limitaciones.")
    print("     Los ejes de diseño (zeta, radio, PEG) van en resultados. No mezclarlos.")
    print("=" * 79)


def figuras(prefijo="envolvimiento"):
    """Guarda las figuras con matplotlib. Solo con --figuras."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Fig 1: G(D) de los tres diseños
    fig, ax = plt.subplots(figsize=(7, 4.5))
    D = np.linspace(E.D0_nm, 6.0, 4000)
    for lip in LIPOSOMAS:
        G = E.energia_libre_J(D, lip["R_nm"], lip["zeta_mV"],
                              lip["peg_nm"], HAMAKER_REF) / E.KT_J
        ax.plot(D, G, label=lip["nombre"])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("separación D (nm)"); ax.set_ylabel("G(D)  [kT]")
    ax.set_title("Energía libre de interacción DLVO")
    ax.set_ylim(-30, 15); ax.legend(); fig.tight_layout()
    fig.savefig(f"{prefijo}_G_de_D.png", dpi=160); plt.close(fig)

    # Fig 2: radio crítico frente a w, para los tres kappa
    fig, ax = plt.subplots(figsize=(7, 4.5))
    w_grid = np.logspace(-5, -1, 400)
    for k in E.KAPPA_kT:
        ax.plot(w_grid * 1e6, [E.radio_critico_nm(w, k) for w in w_grid],
                label=f"kappa = {k:.0f} kT")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("w  (uN/m)"); ax.set_ylabel("R_min  (nm)")
    ax.set_title("Radio crítico de envolvimiento  R_min = sqrt(2k/w)")
    ax.legend(); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(f"{prefijo}_radio_critico.png", dpi=160); plt.close(fig)

    # Fig 3: barrera de envolvimiento frente a sigma_tilde
    fig, ax = plt.subplots(figsize=(7, 4.5))
    st = np.logspace(-3, 3, 600)
    for k in E.KAPPA_kT:
        ax.plot(st, E.barrera_envolvimiento_kT(st, k), label=f"kappa = {k:.0f} kT")
    ax.axhline(3.0, color="gray", ls="--", lw=0.8)
    ax.text(1e-3, 3.4, "~3 kT: cruzable térmicamente", fontsize=8, color="gray")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("tensión reducida  sigma_tilde"); ax.set_ylabel("barrera  (kT)")
    ax.set_title("Barrera de la transición de envolvimiento (orden de magnitud)")
    ax.legend(); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(f"{prefijo}_barrera.png", dpi=160); plt.close(fig)

    print(f"\n Figuras guardadas: {prefijo}_G_de_D.png, "
          f"{prefijo}_radio_critico.png, {prefijo}_barrera.png")


def main():
    ap = argparse.ArgumentParser(description="Módulo de envolvimiento · BHE/SENACYT")
    ap.add_argument("--solo-tests", action="store_true", help="solo la validación")
    ap.add_argument("--figuras", action="store_true", help="guarda las figuras PNG")
    args = ap.parse_args()

    if args.solo_tests:
        ok = E.test_limites()
        raise SystemExit(0 if ok else 1)

    cabecera()
    print()
    E.test_limites()
    tabla_compuertas()
    tabla_sensibilidad()
    lectura()
    if args.figuras:
        figuras()


if __name__ == "__main__":
    main()
