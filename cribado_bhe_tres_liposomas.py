#@title Cribado BHE · 3 liposomas de fingolimod x corona de proteínas x barrido de Hamaker
#
# =============================================================================
#  SUPERADO por envolvimiento/envolvimiento_core.py, marcado el 2026-08-16.
#  NO USAR PARA RESULTADOS PUBLICADOS. Se conserva solo como registro histórico
#  (aquí nació la idea del cribado). No alimenta la web ni las figuras del
#  pipeline (correr.sh / construir_web.py), que usan envolvimiento_core.py.
#
#  Dos motivos, encontrados en la auditoría de trazabilidad de Jhovan:
#   1. El van der Waals NO escala con el radio (PENDIENTES.md, tarea P5):
#      `_derivar()` calibra D_sat para que el pozo en el contacto valga
#      siempre 15 kT, PARA CUALQUIER R. El vdW esfera-plano debe escalar
#      proporcional a R. Cualquier conclusión sobre TAMAÑO que salga de este
#      script es un artefacto de esa calibración, no física.
#      envolvimiento_core.py corrige esto: usa solo el corte D0=0.165 nm.
#   2. ZETA_CORONA_mV=-7.0 (más abajo) está etiquetado como si fuera un valor
#      medido por Walter et al. 2021, y no lo es: Walter respalda el
#      MECANISMO (la corona invierte la carga hacia negativo), no esta cifra
#      concreta, que es un parámetro de simulación propio. Corregido el
#      comentario, pero el origen del número sigue siendo una calibración,
#      no una medida.
# =============================================================================
#
# AUTOCONTENIDO: pégalo en Colab y ejecuta. No necesita ningún otro archivo (solo numpy).
#
# ---------------------------------------------------------------------------------------------
# QUÉ HACE ESTE CÓDIGO
# ---------------------------------------------------------------------------------------------
# Responde a una sola pregunta: ¿se PEGA el liposoma a la pared interna del vaso cerebral?
#
# Ese es el primer paso, necesario pero no suficiente, para que el nanotransportador cruce la
# barrera hematoencefálica (BHE). Si no se adhiere a la superficie luminal del endotelio, la célula
# nunca lo va a internalizar. Lo que este código NO hace es simular el cruce completo (la célula
# tragándose la vesícula y soltándola al otro lado): eso es transporte activo y requiere biología,
# no solo energía libre.
#
# CÓMO LO CALCULA
# Se calcula la energía libre G(D) del sistema en función de D = separación entre la superficie del
# liposoma y la membrana, con teoría DLVO extendida. G(D) es la suma de cuatro fuerzas:
#   1. van der Waals      -> siempre ATRAE. Depende del tamaño y de la constante de Hamaker.
#   2. Doble capa eléctrica -> ATRAE si el liposoma y el endotelio tienen carga opuesta,
#                              REPELE si tienen el mismo signo. Aquí el endotelio es negativo
#                              (-11.4 mV) y tus tres liposomas son positivos -> atrae.
#   3. Hidratación        -> siempre REPELE, de muy corto alcance (las dos superficies están
#                            rodeadas de agua que cuesta desplazar).
#   4. Estérica del PEG   -> REPELE, solo si el liposoma está PEGilado. El cepillo de PEG actúa
#                            como un colchón que estorba el contacto.
#
# Con G(D) construida, la fracción adherida en el equilibrio sale EXACTA de la ley de Boltzmann:
#
#     f = ∫[contacto] exp(-G/kT) dD  /  ∫[todo el rango] exp(-G/kT) dD
#
# No hace falta simular partículas moviéndose: la distribución de equilibrio de ese sistema ES
# exp(-G/kT), así que integrarla da el mismo número, exacto y en milisegundos. (Una versión previa
# sí usaba dinámica Browniana y tardaba minutos en equilibrar los pozos profundos.)
#
# QUÉ SE BARRE Y POR QUÉ
#  · La constante de Hamaker no está medida para tu composición lipídica, así que en vez de fijar
#    un valor se prueban tres. Si la conclusión no cambia entre ellos, el resultado es robusto sin
#    necesidad del valor exacto.
#  · Cada diseño se evalúa con su zeta nominal Y con un zeta post-corona: al entrar en sangre, las
#    proteínas del plasma se adsorben sobre la nanopartícula y desplazan su potencial zeta hacia el
#    rango negativo INDEPENDIENTEMENTE de su carga original. El MECANISMO tiene respaldo en Walter
#    et al. 2021 (Tissue Barriers 9:1904773, Sección 3), pero la CIFRA de -7.0 mV es un PARÁMETRO DE
#    SIMULACIÓN propio de este script, no un valor medido que Walter reporte para liposomas con
#    corona proteica. No citar -7.0 mV como dato de Walter 2021 (auditoría de Jhovan, 2026-08-16).
#
# CONTROL DE FALSABILIDAD
# Cada caso se repite con la carga del liposoma invertida a fuertemente aniónica (-30 mV). Como en
# estos diseños la adhesión la impulsa la atracción electrostática, ese control debe dar adhesión
# mucho menor. Si no bajara, el modelo estaría roto.
# ---------------------------------------------------------------------------------------------
import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz

# ---- constantes físicas (SI) ----
EPS0 = 8.8541878128e-12   # F/m
KB   = 1.380649e-23       # J/K
QE   = 1.602176634e-19    # C
NM   = 1e-9               # m por nm


# =============================================================================
#  PARÁMETROS DE LA BHE Y DEL MODELO
# =============================================================================
ZETA_BHE_mV = -11.4   # endotelio cerebral humano hCMEC/D3 a 37 C.
                      # Santa-Maria et al. 2019, BBA Biomembranes 1861:1579, Figura 4A.
                      # Corroborado: -12.7 mV en Kincses et al. 2020, Lab Chip 20:3792, Fig. 6F.
                      # Rango defendible entre modelos celulares: -11 a -15 mV.

TEMP_K   = 310.0      # K (37 C)
EPS_R    = 74.0       # permitividad relativa del agua a 37 C (Malmberg-Maryott da 74.15)
DEBYE_nm = 0.8        # longitud de Debye a ~150 mM. Calculable:
                      # lambda_D = sqrt(eps0*eps_r*kB*T/(2*NA*e^2*I)) -> 0.78 nm a 150 mM y 37 C
D0_nm    = 0.165      # separación de contacto físico (convención DLVO)
D_MAX_nm = 10.0       # frontera lejana (más allá el gradiente es ~0)
W_ADH_nm = 1.0        # ventana desde el contacto que cuenta como "adherido"
B_HID_kT = 1.5        # amplitud de la repulsión de hidratación
W_HID_nm = 0.3        # alcance de la hidratación
B_PEG_kT = 8.0        # amplitud de la repulsión estérica del cepillo de PEG
POZO_VDW_OBJETIVO_kT = 15.0   # ver nota en _derivar()

KT_J = KB * TEMP_K


def _derivar(R_ext_nm, hamaker_J):
    """Longitud de saturación del van der Waals.

    El vdW puro (-A_H*R/6D) diverge cuando D->0, lo cual no es físico (a escala atómica el modelo
    continuo deja de valer) y da pozos de cientos de kT. Se regulariza suavizando el denominador
    (D -> D + D_sat), con D_sat calibrada para que el pozo en el contacto tenga una profundidad
    finita y razonable. Es una calibración explícita, no un valor medido.
    """
    R_m = R_ext_nm * NM
    objetivo_J = POZO_VDW_OBJETIVO_kT * KT_J
    D_sat_m = (hamaker_J * R_m) / (6.0 * objetivo_J) - (D0_nm * NM)
    return R_m, max(0.0, D_sat_m / NM)


def energia_libre(D_nm, R_ext_nm, zeta_nano_mV, peg_nm, hamaker_J, zeta_bhe_mV=ZETA_BHE_mV):
    """G(D) en unidades de kT. D = separación superficie-a-superficie en nm."""
    D = np.maximum(np.asarray(D_nm, dtype=float), D0_nm)
    R_m, D_sat_nm = _derivar(R_ext_nm, hamaker_J)
    D_m = D * NM

    # 1) van der Waals (atractivo), Derjaguin esfera-plano, con saturación de contacto
    g_vdw = (-hamaker_J * R_m / (6.0 * (D_m + D_sat_nm * NM))) / KT_J

    # 2) Doble capa eléctrica (aprox. de superposición lineal). El SIGNO sale solo:
    #    cargas opuestas -> gamma1*gamma2 < 0 -> término negativo (atractivo).
    kT_e = KT_J / QE                      # V
    g1 = np.tanh((zeta_nano_mV * 1e-3) / (4.0 * kT_e))
    g2 = np.tanh((zeta_bhe_mV  * 1e-3) / (4.0 * kT_e))
    pref_J = 64.0 * np.pi * EPS0 * EPS_R * R_m * (kT_e ** 2)
    g_edl = (pref_J * g1 * g2 * np.exp(-D_m / (DEBYE_nm * NM))) / KT_J

    # 3) Hidratación (repulsiva, corto alcance) y 4) cepillo de PEG (repulsivo, si lo hay)
    g_hid = B_HID_kT * np.exp(-(D - D0_nm) / W_HID_nm)
    g_peg = B_PEG_kT * np.exp(-(D - D0_nm) / peg_nm) if peg_nm > 0 else 0.0

    return g_vdw + g_edl + g_hid + g_peg


def fraccion_equilibrio(R_ext_nm, zeta_nano_mV, peg_nm, hamaker_J):
    """Fracción adherida en el equilibrio, exacta desde Boltzmann."""
    D = np.linspace(D0_nm, D_MAX_nm, 20000)
    G = energia_libre(D, R_ext_nm, zeta_nano_mV, peg_nm, hamaker_J)
    w = np.exp(-(G - G.min()))           # se resta el mínimo por estabilidad numérica
    dentro = D < (D0_nm + W_ADH_nm)
    return float(_trapz(w[dentro], D[dentro]) / _trapz(w, D))


# =============================================================================
#  FICHAS: valores de DISEÑO del estudio (zeta asignado en un rango, NO medido)
# =============================================================================
LIPOSOMAS = [
    dict(nombre="Convencional", R_ext_nm=20.0, zeta_mV=+5.0, peg_nm=0.0, frac_encaps=0.78),
    dict(nombre="Furtivo/PEG",  R_ext_nm=15.5, zeta_mV=+2.0, peg_nm=5.0, frac_encaps=0.69),
    dict(nombre="Catiónico",    R_ext_nm=17.5, zeta_mV=+6.7, peg_nm=0.0, frac_encaps=0.76),
]

HAMAKER_J      = [3.0e-21, 4.5e-21, 6.5e-21]   # sin valor primario verificado -> se barre
ZETA_CORONA_mV = -7.0                          # zeta tras adsorción de proteínas séricas: PARÁMETRO
                                                # DE SIMULACIÓN, no medido; mecanismo respaldado por
                                                # Walter et al. 2021 Sec. 3, la cifra NO
ZETA_CONTROL_mV = -30.0                        # control de falsabilidad


def barra(frac, ancho=20):
    n = int(round(frac * ancho))
    return "#" * n + "." * (ancho - n)


if __name__ == "__main__":
    print("=" * 78)
    print(" CRIBADO DE ADHESIÓN A LA BHE: 3 diseños de liposoma con fingolimod")
    print("=" * 78)
    print(f" Endotelio: zeta = {ZETA_BHE_mV} mV (hCMEC/D3 humano, 37 C; Santa-Maria 2019 Fig. 4A)")
    print(f" Escenario corona: zeta del liposoma -> {ZETA_CORONA_mV:+.1f} mV "
          f"(parámetro de simulación; mecanismo de Walter 2021 Sec. 3, cifra no medida por esa fuente)")
    print(f" Barrido de Hamaker: {[f'{h:.1e}' for h in HAMAKER_J]} J (sin valor primario verificado)")
    print(" Ninguno de los tres diseños lleva ligando dirigido.")
    print("=" * 78)

    resumen = []
    for lip in LIPOSOMAS:
        diam = 2 * lip["R_ext_nm"]
        cabe = "CABE" if diam <= 60 else ("LÍMITE" if diam <= 80 else "NO CABE")
        print(f"\n--- {lip['nombre']}  (R_ext {lip['R_ext_nm']:.1f} nm, "
              f"zeta ficha {lip['zeta_mV']:+.1f} mV, PEG {lip['peg_nm']:.1f} nm) ---")
        print(f"    Compuerta de tamaño: Ø {diam:.0f} nm vs caveolas 60-80 nm -> {cabe}")
        for etiqueta, z in (("nominal", lip["zeta_mV"]), ("corona ", ZETA_CORONA_mV)):
            for h in HAMAKER_J:
                f_adh  = fraccion_equilibrio(lip["R_ext_nm"], z, lip["peg_nm"], h)
                f_ctrl = fraccion_equilibrio(lip["R_ext_nm"], ZETA_CONTROL_mV, lip["peg_nm"], h)
                print(f"    {etiqueta} | A_H={h:.1e} J | adherido {f_adh*100:5.1f}% "
                      f"[{barra(f_adh)}] | control {f_ctrl*100:5.1f}%")
                resumen.append((lip["nombre"], etiqueta.strip(), f_adh))

    print("\n" + "=" * 78); print(" LECTURA DEL CRIBADO"); print("=" * 78)
    for nombre in [l["nombre"] for l in LIPOSOMAS]:
        nom = [r[2] for r in resumen if r[0] == nombre and r[1] == "nominal"]
        cor = [r[2] for r in resumen if r[0] == nombre and r[1] == "corona"]
        print(f" {nombre:14s} nominal {min(nom)*100:5.1f}-{max(nom)*100:5.1f}%   "
              f"corona {min(cor)*100:5.1f}-{max(cor)*100:5.1f}%")
    print("-" * 78)
    print(" El rango de cada celda es la sensibilidad al Hamaker. Si la conclusión (adhiere / no")
    print(" adhiere) NO cambia dentro del rango, el resultado es robusto pese a no tener el Hamaker")
    print(" medido, y eso es más fuerte que un valor puntual sin respaldo.")
    print("=" * 78)
    print(" LÍMITES QUE HAY QUE DECLARAR EN EL ARTÍCULO:")
    print(" 1. Los zeta son valores de DISEÑO asignados en un rango, no medidas. Esto es un cribado")
    print("    del espacio de diseño, no una predicción sobre formulaciones concretas.")
    print(" 2. El modelo predice ADHESIÓN luminal: paso NECESARIO pero NO suficiente para el cruce.")
    print("    Adhesión favorable != entrega cerebral.")
    print(" 3. El glicocálix mide 0.2-5 um (Walter 2021), de 6 a 160 veces más grueso que estos")
    print("    liposomas. El modelo trata la superficie como plano liso y NO modela la penetración")
    print("    de ese cepillo, que plausiblemente es el paso limitante real.")
    print(" 4. In vivo los catiónicos tienen aclaramiento rápido y alta captación hepática, y llegan")
    print("    PEOR al cerebro que los neutros/ligeramente negativos (Walter 2021). Si el catiónico")
    print("    gana en adhesión, decláralo explícitamente.")
    print("=" * 78)
