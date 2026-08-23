#!/usr/bin/env python3
# =============================================================================
#  GENERADOR de la versión para Google Colab
#
#  Por qué existe este archivo: mantener a mano dos copias de la misma física
#  (una para script y otra para Colab) es la forma más rápida de que diverjan y
#  de que un arreglo se aplique solo en una. Aquí la versión de Colab se GENERA
#  concatenando envolvimiento_core.py + envolvimiento_script.py, de modo que
#  siempre son el mismo código.
#
#  Uso:
#      python3 construir_colab.py
#  Produce: envolvimiento_colab.py  (autocontenido, pégalo en una celda)
#
#  Después de generarlo, verifica_equivalencia.py comprueba que las dos
#  versiones dan exactamente los mismos números.
# =============================================================================

from pathlib import Path

AQUI = Path(__file__).resolve().parent
CORE = AQUI / "envolvimiento_core.py"
SCRIPT = AQUI / "envolvimiento_script.py"
SALIDA = AQUI / "envolvimiento_colab.py"

ENCABEZADO = '''#@title Envolvimiento de membrana · BHE/SENACYT (Fase 2: tareas 2.4, 2.5, 2.6)
# =============================================================================
#  ARCHIVO GENERADO AUTOMÁTICAMENTE: NO LO EDITES A MANO.
#  Se genera con construir_colab.py a partir de:
#      envolvimiento_core.py    (la física)
#      envolvimiento_script.py  (la presentación)
#  Si necesitas cambiar algo, cámbialo allí y vuelve a generar este archivo.
#
#  AUTOCONTENIDO: pégalo en una celda de Colab y ejecútalo. Solo necesita numpy
#  y matplotlib, que Colab ya trae.
# =============================================================================

'''

PIE = '''

# =============================================================================
#  EJECUCIÓN EN COLAB
# =============================================================================
cabecera()
print()
test_limites()
tabla_compuertas()
tabla_sensibilidad()
lectura()

# Figuras en línea (en Colab no hace falta el backend Agg)
try:
    import matplotlib
    import matplotlib.pyplot as plt

    D = np.linspace(D0_nm, 6.0, 4000)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.2))

    for lip in LIPOSOMAS:
        G = energia_libre_J(D, lip["R_nm"], lip["zeta_mV"],
                            lip["peg_nm"], HAMAKER_REF) / KT_J
        axs[0].plot(D, G, label=lip["nombre"])
    axs[0].axhline(0, color="k", lw=0.6)
    axs[0].set_xlabel("separación D (nm)"); axs[0].set_ylabel("G(D)  [kT]")
    axs[0].set_title("Energía libre DLVO"); axs[0].set_ylim(-30, 15); axs[0].legend()

    w_grid = np.logspace(-5, -1, 400)
    for k in KAPPA_kT:
        axs[1].plot(w_grid * 1e6, [radio_critico_nm(w, k) for w in w_grid],
                    label=f"kappa = {k:.0f} kT")
    axs[1].set_xscale("log"); axs[1].set_yscale("log")
    axs[1].set_xlabel("w  (uN/m)"); axs[1].set_ylabel("R_min  (nm)")
    axs[1].set_title("Radio crítico  R_min = sqrt(2k/w)")
    axs[1].legend(); axs[1].grid(alpha=0.3, which="both")

    st = np.logspace(-3, 3, 600)
    for k in KAPPA_kT:
        axs[2].plot(st, barrera_envolvimiento_kT(st, k), label=f"kappa = {k:.0f} kT")
    axs[2].axhline(3.0, color="gray", ls="--", lw=0.8)
    axs[2].set_xscale("log"); axs[2].set_yscale("log")
    axs[2].set_xlabel("sigma_tilde"); axs[2].set_ylabel("barrera (kT)")
    axs[2].set_title("Barrera de envolvimiento (orden de magnitud)")
    axs[2].legend(); axs[2].grid(alpha=0.3, which="both")

    plt.tight_layout(); plt.show()
except Exception as _e:
    print(f"[aviso] no se pudieron dibujar las figuras: {_e}")
'''


def limpiar_script(texto: str) -> str:
    """Quita del script lo que no aplica en una celda de Colab."""
    fuera = []
    saltando_main = False
    for linea in texto.splitlines():
        s = linea.strip()
        # el core ya está pegado arriba: fuera los imports del módulo
        if s in ("import envolvimiento_core as E", "import argparse"):
            continue
        # el bloque main() y el guard no se usan en Colab
        if s.startswith("def main("):
            saltando_main = True
            continue
        if saltando_main:
            if linea and not linea.startswith((" ", "\t")):
                saltando_main = False
            else:
                continue
        if s.startswith('if __name__ == "__main__"') or s == "main()":
            continue
        # la función figuras() se sustituye por el bloque en línea del pie
        fuera.append(linea)
    texto = "\n".join(fuera)
    # como el core queda en el mismo espacio de nombres, "E." sobra
    texto = texto.replace("E.", "")
    return texto


def main():
    core = CORE.read_text(encoding="utf-8")
    script = limpiar_script(SCRIPT.read_text(encoding="utf-8"))
    # el core ya importa numpy; quitamos el import duplicado del script
    script = script.replace("import numpy as np\n", "", 1)

    SALIDA.write_text(ENCABEZADO + core + "\n\n" + script + PIE, encoding="utf-8")
    n = len(SALIDA.read_text(encoding="utf-8").splitlines())
    print(f"Generado: {SALIDA.name}  ({n} líneas)")
    print("Verifica con: python3 verifica_equivalencia.py")


if __name__ == "__main__":
    main()
