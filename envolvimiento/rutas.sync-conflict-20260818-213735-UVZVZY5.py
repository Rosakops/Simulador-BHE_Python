#!/usr/bin/env python3
# =============================================================================
#  SIMULADOR DE RUTAS COMPLETAS  —  proyecto BHE / SENACYT
#  Versión BINARIA (v1). Diseñada para crecer a cuantitativa sin rehacerse.
# =============================================================================
#
#  QUÉ HACE
#  --------
#  Toma un diseño de nanotransportador (tamaño, carga, PEG) y lo hace recorrer
#  cada ruta candidata hasta la mielina cerebral. Cada ruta es una secuencia
#  ordenada de COMPUERTAS. Si el diseño falla una sola, la ruta queda EXCLUIDA.
#
#  QUÉ AFIRMA Y QUÉ NO
#  -------------------
#  Este simulador está construido para EXCLUIR, no para predecir.
#    · "EXCLUIDA"      -> el diseño NO puede usar esa ruta. Afirmación FUERTE.
#    · "NO EXCLUIDA"   -> el diseño es CANDIDATO por esa ruta. Afirmación DÉBIL.
#    · "NO EVALUABLE"  -> falta información para decidir. NO es un aprobado.
#
#  La tercera categoría es deliberada y es lo que separa este simulador de uno
#  que solo sabe decir que sí. Una compuerta sin dato NO se da por superada.
#
#  CÓMO CRECE A CUANTITATIVO
#  -------------------------
#  Cada compuerta ya devuelve `valor`, `umbral` y `margen`. La versión
#  cuantitativa sustituirá el booleano por una probabilidad de superar la
#  compuerta y multiplicará a lo largo de la ruta. La estructura no cambia.
#
#  ANCLAJES EXPERIMENTALES
#  -----------------------
#  Cada compuerta implementada reproduce al menos un dato publicado. Ver
#  validar_contra_experimentos().
# =============================================================================

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

import envolvimiento_core as E
import glicocalix as G

PASA, FALLA, DESCONOCIDA = "PASA", "FALLA", "DESCONOCIDA"


# =============================================================================
#  ESTRUCTURAS
# =============================================================================

@dataclass
class Diseno:
    """Un nanotransportador candidato. Los tamaños son DIÁMETRO en nm.

    `clase` importa y no es decorativa: el suelo geométrico del transportador
    depende de su arquitectura. Un liposoma es una bicapa cerrada y no puede
    medir menos de dos espesores de bicapa; un dendrímero es una molécula
    ramificada maciza y no tiene ese suelo. Aplicar el límite del liposoma a un
    dendrímero lo excluiría por un motivo que no le afecta.
    """
    nombre: str
    diametro_nm: float
    zeta_mV: float
    peg_nm: float = 0.0
    # Fingolimod. DERIVADO, no medido: 1.683 nm es la dimensión máxima extremo a
    # extremo, promedio de 50 confórmeros ETKDGv3 + MMFF94. Ver
    # `tamano_farmaco.py`, que es reproducible y explica por qué se toma la
    # dimensión máxima y no el diámetro de esfera equivalente (0.858 nm).
    # Sustituye al «~1.0 nm sin fuente» que había hasta el 2026-08-13 [tarea C9].
    farmaco_diametro_nm: float = 1.683
    clase: str = "liposoma"            # liposoma | dendrimero | polimerico | micela
    nota: str = ""
    # Solo para clase="dendrimero" (tarea G.1a-bis). La ventana geométrica
    # suelo/techo está derivada de fuente primaria de PAMAM y SOLO de PAMAM.
    # PPI y carbosilano son otra química: heredarles esa ventana sería
    # transferencia de parámetro entre sistemas distintos, que es justo lo que
    # se declaró DESCONOCIDA en C.2 y D.3. Sin subquímica se asume PAMAM.
    subquimica: str = "pamam"          # pamam | ppi | carbosilano
    # True = el dato de entrada NO viene de medida ni de fuente primaria: es un
    # valor inventado dentro de un rango teórico plausible, para probar el
    # simulador. Nunca puede citarse como resultado.
    sintetico: bool = False
    # Solo para clase="polimerico" (tarea G.1b). Su suelo NO es un número fijo
    # sino una función de la masa molar de la cadena, así que el diseño tiene
    # que declararla. Sin ella la compuerta devuelve DESCONOCIDA.
    masa_molar_kDa: Optional[float] = None
    densidad_g_cm3: Optional[float] = None

    @property
    def radio_nm(self):
        return self.diametro_nm / 2.0


@dataclass
class Resultado:
    """Resultado de una compuerta.

    `advertencia` es la salvedad del respaldo, que no debe perderse. Su uso
    principal es el PASA: un PASA con advertencia NO es un PASA limpio y se
    marca en la salida con "✓!" y con punto hueco en la figura del recorrido.

    Una DESCONOCIDA también puede llevarla (B.3 desde el 2026-08-12), y ahí
    documenta POR QUÉ el dato disponible no sirve. En ese caso NO se pinta
    hueca: hueco significa aprobado con salvedad, que es lo contrario.
    """
    compuerta: str
    estado: str
    valor: Optional[float] = None
    umbral: Optional[float] = None
    unidad: str = ""
    margen: Optional[float] = None
    fuente: str = ""
    motivo: str = ""
    advertencia: str = ""


def _cmp(nombre, valor, umbral, unidad, fuente, mayor_es_mejor):
    """Compuerta numérica simple. margen > 0 significa que pasa con holgura."""
    margen = (valor - umbral) if mayor_es_mejor else (umbral - valor)
    return Resultado(nombre, PASA if margen > 0 else FALLA,
                     valor, umbral, unidad, margen, fuente)


# =============================================================================
#  COMPUERTAS IMPLEMENTADAS  (cada una con su anclaje)
# =============================================================================

# -----------------------------------------------------------------------------
#  LÍMITE GEOMÉTRICO DEL DENDRÍMERO PAMAM  (tarea G.1a, cerrada 2026-08-12)
#
#  Al contrario que el liposoma, el dendrímero tiene SUELO Y TECHO, y el techo
#  es el que manda. Ficha completa en
#  verificacion/verificacion_dendrimero_tarea_G_1a.md
#
#  SUELO = G3. Es la generación más pequeña que cumple a la vez las dos
#  condiciones del criterio ("existe Y lleva el fármaco dentro"):
#    · diámetro MEDIDO en fuente primaria — Prosa et al. 2001, Macromolecules
#      34:4897, Tabla 1: R = 18.8 Å por SAXS en metanol, o sea 3.76 nm.
#    · alojamiento demostrado de un fármaco hidrófobo — Devarakonda et al. 2004,
#      Int J Pharm 284:133, Tabla 2: complejo 1:1 con nifedipino, K = 287.6 M⁻¹
#      a pH 7, un orden de magnitud por encima de G1 (25.6) y G2 (52.7).
#  Por DEBAJO del suelo la compuerta NO falla: devuelve DESCONOCIDA. G1 y G2 sí
#  forman complejo 1:1, lo que falta es su diámetro medido (Prosa no bajó de G3).
#
#  TECHO = G10, la última generación completable. Tres líneas independientes:
#    · de Gennes & Hervet 1983, J. Physique Lett. 44:351: m1 = 2.88(ln P + 1.5),
#      con P = 7 da 9.9. Maiti lo escribe como 2.88·ln P + 4.4 y da 10.2.
#    · Maiti et al. 2004, Macromolecules 37:6236, Fig. 21: la energía de tensión
#      por monómero se dispara en G10-G11. Textual: "generation 11 is already
#      past the limiting generation".
#    · Maiti 2004, Fig. 22: el área disponible por amina terminal cae por debajo
#      de los 28.4 Å² que ocupa una amina aislada entre G10 y G11.
#  Valor: 14.16 nm (Maiti, Tabla 4, R_SAV = 70.80 Å). Se usa el MAYOR de los dos
#  disponibles a propósito, para no excluir de más; el medido por Prosa es
#  13.98 nm (R = 69.9 Å) y coincide con el calculado dentro del 1.3 %.
# -----------------------------------------------------------------------------

# Diámetros MEDIDOS por generación. Prosa et al. 2001, Macromolecules 34:4897,
# Tabla 1, columna R del ajuste a distribución gaussiana de esferas (SAXS en
# metanol, PAMAM de núcleo etilendiamina). Diámetro = 2R. Prosa no midió G0–G2.
DENDRIMERO_GENERACIONES_nm = {
    3: 3.76, 4: 4.60, 5: 5.64, 6: 7.26, 7: 8.38, 8: 10.04, 9: 11.84, 10: 13.98,
}

DENDRIMERO_SUELO_nm = 3.76          # G3  — Prosa 2001, Tabla 1 (medido, SAXS)
DENDRIMERO_TECHO_nm = 14.16         # G10 — Maiti 2004, Tabla 4 (R_SAV)
DENDRIMERO_TECHO_MEDIDO_nm = 13.98  # G10 — Prosa 2001, Tabla 1 (medido, SAXS)
DENDRIMERO_PRECISION = 0.05         # ±5 % de precisión global, Prosa 2001

_AV_CARGA = ("el único dato de estequiometría en fuente primaria (Devarakonda "
             "2004, Tabla 2, generaciones 0.5 a 3) es 1:1, UNA molécula por "
             "dendrímero. OJO CON TRES COSAS, comprobadas el 2026-08-13 sobre "
             "el cuerpo del artículo: (1) el fármaco es NIFEDIPINO, el "
             "fingolimod no aparece ni una vez; (2) el 1:1 NO se mide, se "
             "INFIERE de que la pendiente del diagrama de solubilidad tipo A_L "
             "es menor que la unidad (Higuchi y Connors), y describe un "
             "complejo de SOLUBILIZACIÓN, no una carga útil máxima; (3) llega "
             "a G3, no a las generaciones altas. Sin dato primario de carga "
             "para G4–G10 ni para fingolimod  [tarea G.2]")

_AV_DISOLVENTE = ("los tamaños están medidos en METANOL (Prosa 2001); en agua "
                  "las aminas terminales protonadas hinchan el dendrímero "
                  "(Maiti estima +15 % de radio en G6), así que el techo real "
                  "podría ser mayor  [tarea G.3]")


def g_transportador_fabricable(d: Diseno):
    """¿Es geométricamente posible un transportador de esa clase y ese tamaño?

    DEPENDE DE LA CLASE, y esto no es un detalle: el suelo del liposoma sale de
    que es una BICAPA CERRADA (d_ext = d_núcleo + 2·t_bicapa, con t medido en
    3.6–4.6 nm), y se exige un núcleo acuoso mínimo de 4 nm para que encapsule
    algo. Un dendrímero es una molécula ramificada maciza y no tiene ese suelo,
    pero sí tiene un TECHO que el liposoma no tiene.

    El polímero macizo tiene suelo pero no techo, y su suelo es función de la
    masa molar (ver `_polimerico_fabricable`). Cualquier otra clase devuelve
    DESCONOCIDA en vez de que se le aplique el límite de otra.
    """
    if d.clase == "liposoma":
        minimo = G.diametro_liposoma_minimo_nm(4.0, 4.0)
        return _cmp("Transportador fabricable", d.diametro_nm, minimo, "nm",
                    "Pan et al. 2008, PRL 100:198103, Fig. 3c", True)

    if d.clase == "micela":
        return _micela_fabricable(d)

    if d.clase == "dendrimero":
        return _dendrimero_fabricable(d)

    if d.clase == "polimerico":
        return _polimerico_fabricable(d)

    return Resultado("Transportador fabricable", DESCONOCIDA,
                     d.diametro_nm, None, "nm", None, "",
                     f"clase '{d.clase}' desconocida para esta compuerta")


# -----------------------------------------------------------------------------
#  POLÍMERO MACIZO  (tarea G.1b, cerrada 2026-08-12)
#  Ficha: verificacion/verificacion_polimero_micela_tarea_G_1b.md
#
#  POLÍMERO MACIZO. No tiene techo, y su suelo NO es un récord de fabricación
#  (que caduca) sino una identidad geométrica: una partícula maciza no puede ser
#  más pequeña que UNA SOLA CADENA del polímero colapsada sobre sí misma.
#      V = M / (ρ·N_A)        d = (6V/π)^(1/3)
#  Las dos magnitudes están medidas. Densidades de Parker et al. 2010, Biomed
#  Mater 5:055004, TABLA 2 (derivadas por los autores de su propia velocidad del
#  sonido e impedancia, ρ = Z/c). OJO: la Tabla 1 de ese mismo artículo da
#  densidades del FABRICANTE y no se usan.
#  El fármaco NO entra en el cálculo, y no por comodidad: una molécula de
#  fingolimod sobre una cadena de 53 kDa mueve el diámetro un 0.19 %. Así el
#  suelo no arrastra el tamaño del fármaco (tarea C9, ya derivado).
#
#  MICELA: fuera de alcance el 2026-08-12, REINCORPORADA el 2026-08-12 a
#  petición de Jhovan, cuando resultó que dos de las tres fichas de polímero
#  teórico que trajo (PEG-PLA y PCL-Pluronic-PCL) son micelas.
#
#  NO ES UNA CUARTA CLASE (decisión de Jhovan, 2026-08-12). El proyecto sigue
#  teniendo TRES clases: liposoma, dendrímero y polimérico. La micela es una
#  ARQUITECTURA dentro de la polimérica: sigue siendo polímero, y por eso vive
#  en la sección de polímeros de la web y de la hipótesis. `clase="micela"` solo
#  existe para que la compuerta del suelo elija la fórmula correcta.
#
#  Una micela ESTÁ hecha de polímero pero NO es un polímero macizo: es un
#  agregado autoensamblado de núcleo y corona, sostenido por el efecto
#  hidrófobo, con concentración micelar crítica por debajo de la cual se
#  deshace. Su suelo no es el glóbulo de UNA cadena colapsada sino el que fija
#  el número de agregación y la carga de fármaco, así que no puede pasar por
#  `_polimerico_fabricable`.
#
#  SUELO = 13.0 nm, micela CARGADA. Sochor et al. 2020, Langmuir 36:3494, SANS:
#  núcleo ≈30 Å + cáscara 1 ≈35 Å, cáscara 2 ≈0 a esa carga. Vacía mide 3.6 nm:
#  el fármaco no es un pasajero, multiplica el diámetro por 3.6.
#  Por DEBAJO del suelo la compuerta devuelve DESCONOCIDA, no FALLA (decisión de
#  la ficha G.1b §5): el dato sale de UN sistema, UN polímero y UN fármaco
#  distintos de los nuestros. Mismo criterio que C.2.
#  Ficha y las cuatro salvedades: verificacion/verificacion_polimero_micela_tarea_G_1b.md
# -----------------------------------------------------------------------------

_N_AVOGADRO = 6.02214076e23

# Parker et al. 2010, Biomed Mater 5:055004, Tabla 2. Densidad en g/cm³.
POLIMERO_DENSIDAD_g_cm3 = {
    "PLGA 85:15": 1.19, "PLA15": 1.14, "PLA24": 1.13, "PLA60": 1.22,
}
POLIMERO_DENSIDAD_POR_DEFECTO = 1.19   # PLGA 85:15, el único PLGA medido

def diametro_globulo_colapsado_nm(masa_molar_g_mol, densidad_g_cm3):
    """Diámetro de una sola cadena de polímero colapsada en esfera compacta."""
    volumen_nm3 = masa_molar_g_mol / (densidad_g_cm3 * _N_AVOGADRO) * 1e21
    return (6.0 * volumen_nm3 / np.pi) ** (1.0 / 3.0)


def _polimerico_fabricable(d: Diseno):
    fuente = ("Parker et al. 2010, Biomed Mater 5:055004, Tabla 2 (densidad) · "
              "glóbulo de cadena colapsada")
    if d.masa_molar_kDa is None:
        return Resultado("Transportador fabricable", DESCONOCIDA,
                         d.diametro_nm, None, "nm", None, fuente,
                         "el suelo del polímero macizo es función de la masa "
                         "molar y el diseño no la declara  [tarea G.1b]")

    rho = d.densidad_g_cm3 or POLIMERO_DENSIDAD_POR_DEFECTO
    minimo = diametro_globulo_colapsado_nm(d.masa_molar_kDa * 1000.0, rho)
    r = _cmp("Transportador fabricable", d.diametro_nm, minimo, "nm", fuente, True)
    if r.estado == PASA:
        r.advertencia = (
            "la densidad la derivan sus autores de impedancia y velocidad del "
            "sonido (ρ = Z/c) y ellos mismos la llaman ESTIMACIÓN, no medida "
            "directa; con el rango medido completo (1.13–1.22 g/cm³) el suelo "
            "varía un 2.5 %  [tarea G.6]")
    else:
        r.motivo = ("por debajo de una sola cadena colapsada: no cabe ni el "
                    "polímero, sin contar el fármaco")
    return r


MICELA_SUELO_CARGADA_nm = 13.0   # Sochor 2020, SANS, micela CARGADA 10/1
MICELA_VACIA_nm = 3.6            # Sochor 2020, la misma micela sin fármaco

_AV_MICELA = (
    "el suelo sale de UN sistema distinto del nuestro: tribloque de "
    "poli(2-oxazolina)/poli(2-oxazina), no PEG-PLA ni PEG-PCL, y cargado con "
    "curcumina (368 Da), no con fingolimod (307 Da). Además los 13.0 nm son "
    "una suma de núcleo + cáscara 1 hecha por nosotros, no una cifra publicada, "
    "y 10/1 es la carga más baja que ellos estudian, no la mínima posible "
    " [tareas G.4 y G.5]")


def _micela_fabricable(d: Diseno):
    """Suelo de la micela CARGADA. No tiene techo arquitectónico conocido."""
    fuente = ("Sochor et al. 2020, Langmuir 36:3494, SANS (radio de núcleo y "
              "cáscaras medidos frente a la carga) · Israelachvili et al. 1976, "
              "J Chem Soc Faraday Trans 2 72:1525 (forma de la restricción)")

    if d.diametro_nm < MICELA_SUELO_CARGADA_nm:
        return Resultado("Transportador fabricable", DESCONOCIDA, d.diametro_nm,
                         MICELA_SUELO_CARGADA_nm, "nm", None, fuente,
                         f"por debajo del suelo de la micela cargada "
                         f"({MICELA_SUELO_CARGADA_nm} nm). NO es un FALLA: el "
                         f"dato viene de otro polímero y otro fármaco, así que "
                         f"no se convierte una extrapolación en veredicto. La "
                         f"micela VACÍA sí mide {MICELA_VACIA_nm} nm, pero "
                         f"vacía no lleva fármaco  [tareas G.4 y G.5]")

    r = _cmp("Transportador fabricable", d.diametro_nm,
             MICELA_SUELO_CARGADA_nm, "nm", fuente, True)
    r.advertencia = _AV_MICELA
    return r


def _dendrimero_fabricable(d: Diseno):
    """Ventana geométrica del PAMAM: SUELO en G3 y TECHO en G10.

    Es una compuerta de DOS lados, así que no puede usar `_cmp`. El `umbral`
    que se reporta es el lado que está decidiendo, y `margen` la distancia a
    ese lado (positiva si pasa).
    """
    dn = d.diametro_nm

    if d.subquimica != "pamam":
        return Resultado("Transportador fabricable", DESCONOCIDA, dn, None,
                         "nm", None, "",
                         f"dendrímero de {d.subquimica.upper()}: la ventana "
                         "suelo/techo del código está derivada solo de PAMAM "
                         "(Prosa 2001, Maiti 2004) y no hay fuente primaria "
                         "del límite geométrico de esta química por generación "
                         " [tarea G.1a-bis]")

    fuente = ("Prosa et al. 2001, Macromolecules 34:4897, Tabla 1 (tamaños) · "
              "Maiti et al. 2004, Macromolecules 37:6236, Tabla 4 y Figs. 21-22 "
              "(techo) · de Gennes & Hervet 1983, J. Physique Lett. 44:351 "
              "(techo, teoría) · Devarakonda et al. 2004, Int J Pharm 284:133 "
              "(alojamiento del fármaco)")

    if dn > DENDRIMERO_TECHO_nm:
        return Resultado("Transportador fabricable", FALLA, dn,
                         DENDRIMERO_TECHO_nm, "nm",
                         DENDRIMERO_TECHO_nm - dn, fuente,
                         "por encima de G10, la última generación completable: "
                         "la tensión estérica impide el crecimiento")

    if dn < DENDRIMERO_SUELO_nm:
        return Resultado("Transportador fabricable", DESCONOCIDA, dn,
                         DENDRIMERO_SUELO_nm, "nm", None, fuente,
                         "por debajo de G3: G1 y G2 sí forman complejo 1:1 con "
                         "un fármaco hidrófobo, pero no hay diámetro medido en "
                         "fuente primaria (Prosa no bajó de G3)")

    margen = min(dn - DENDRIMERO_SUELO_nm, DENDRIMERO_TECHO_nm - dn)
    umbral = (DENDRIMERO_SUELO_nm
              if (dn - DENDRIMERO_SUELO_nm) < (DENDRIMERO_TECHO_nm - dn)
              else DENDRIMERO_TECHO_nm)

    av = [_AV_CARGA]
    if dn > DENDRIMERO_TECHO_MEDIDO_nm * (1.0 - DENDRIMERO_PRECISION):
        av.append(_AV_DISOLVENTE)
    return Resultado("Transportador fabricable", PASA, dn, umbral, "nm",
                     margen, fuente, "", " · ".join(av))


def g_glicocalix_tamiz(d: Diseno):
    """¿Atraviesa el tamiz de la matriz de fibras del glicocálix?"""
    dmax = 2.0 * G.radio_exclusion_nm()
    return _cmp("Tamiz del glicocálix", d.diametro_nm, dmax, "nm",
                "Weinbaum et al. 2003, PNAS 100:7988, Transport Model", False)


def g_envolvimiento(d: Diseno, kappa_kT=25.0, sigma_mNm=0.03, hamaker_J=4.5e-21):
    """¿Es lo bastante grande para que la membrana lo envuelva?"""
    w = E.w_adhesion(d.radio_nm, d.zeta_mV, d.peg_nm, hamaker_J)
    d_min = 2.0 * E.radio_critico_nm(w, kappa_kT)
    return _cmp("Envolvimiento de membrana", d.diametro_nm, d_min, "nm",
                "Deserno 2004, PRE 69:031903, Sec. III C", True)


def g_caveola(d: Diseno):
    """¿Cabe en una caveola?"""
    return _cmp("Compuerta de caveola", d.diametro_nm, E.DIAM_CAVEOLA_MAX_nm, "nm",
                "Bastiani & Parton 2010, J Cell Sci 123:3831", False)


# -----------------------------------------------------------------------------
#  DIFUSIÓN EN EL ESPACIO EXTRACELULAR — compuerta de DOS VARIABLES
# -----------------------------------------------------------------------------
#  Reescrita el 2026-08-10 (tarea V-5). Antes era solo de tamaño, con el 38 nm
#  de Thorne & Nicholson 2006. Nance et al. 2012 impugna ese número: compraron
#  el MISMO lote comercial de puntos cuánticos que usó Thorne (35 nm, −5.1 mV),
#  vieron que no difundían rápido, los recubrieron más con PEG (34 nm, −3.1 mV)
#  y entonces sí. Su conclusión: el rango de Thorne está SUBESTIMADO por
#  adhesión, y lo que decide no es solo el tamaño sino el SIGILO SUPERFICIAL.
#
#  Dato tajante de Nance: el 100 % de las partículas COOH quedaron inmovilizadas
#  o fuertemente obstaculizadas, INCLUIDAS LAS DE 40 nm. Pequeña pero adhesiva
#  no difunde. El tamaño por sí solo no decide.
#
#  Esta compuerta reproduce el diagrama de fases de su Figura 4, que es
#  literalmente "ζ-potential and size vs transport behavior".
#
#  LO QUE NO SE PUEDE EVALUAR: Nance da además un criterio de DENSIDAD de PEG
#  (~9 cadenas de 5 kDa por 100 nm², Γ/SA ≥ 2). El campo `peg_nm` del diseño es
#  un ESPESOR, no una densidad de injerto, y no hay conversión sin más datos.
#  Se usa ζ como sustituto, que es lo que hace el propio Nance en su Fig. 4.
# -----------------------------------------------------------------------------

# Eje de superficie (Nance 2012, Resultados, "High PEG surface density required")
ZETA_DIFUSIVO_mV = -4.0      # menos negativo que esto: difunden 5 de 5
ZETA_NO_DIFUSIVO_mV = -6.0   # más negativo que esto: difunden 0 de 6
                             # entre ambos: 1 de 2  -> zona gris

# Rama POSITIVA del eje de superficie. Nance 2012 no midió ni un ζ positivo
# (su rango es −2.5 a −52 mV), pero la literatura de vectores génicos no
# virales sí, y con la misma técnica (MPT en rebanada de cerebro de rata
# ex vivo). Los dos únicos puntos positivos medidos son:
#   Berry 2016      DNA-UPN   108 ± 13 nm   ζ = +10.0 ± 1.2 mV  -> <10 % difunde
#   Mastorakos 2016 PBAE-CP   120 ± 3.6 nm  ζ = +35.3 ± 1.6 mV  -> inmovilizada
# ZETA_ADHESIVO_POSITIVO_mV es el MÍNIMO positivo medido, no un umbral físico.
# Por debajo NO se extrapola: la compuerta devuelve DESCONOCIDA, mismo criterio
# que en G.1a-bis (dendrímero no PAMAM), C.2 y D.3.
ZETA_ADHESIVO_POSITIVO_mV = 10.0

# Eje de tamaño, PARA SUPERFICIE SIGILOSA (Nance 2012, Tabla 2 y Discusión)
D_SIGILOSO_PASA_nm = 114.0   # demostrado que penetra
D_SIGILOSO_FALLA_nm = 200.0  # demostrado que no

# --- DEPENDENCIA CON LA EDAD del tejido. McKenna et al. 2021, ACS Nano
# 15:8559, Fig. 2B, con PS-PEG de 40 nm (51 nm hidrodinámico) en corteza de
# rata, ajuste de Amsden. NO gobierna ninguna compuerta todavía: se declara
# para poder corregir datos tomados en tejido neonatal y para documentar que
# los umbrales de Nance (corteza humana ADULTA) y los de la escuela de Nance
# en Washington (rata P14) NO son directamente comparables. Tarea E.1.
EDAD_RAZON_DIFUSION = {"P14": 5.0, "P21": 9.0, "P28": 12.0,
                       "P35": 18.0, "P70": 34.0}          # D_ACSF / D_b,eff
EDAD_PORO_EFECTIVO_nm = {"P14": 76.8, "P70": 36.0}        # poro medio efectivo
EDAD_FACTOR_P14_A_P70 = (EDAD_RAZON_DIFUSION["P70"]
                         / EDAD_RAZON_DIFUSION["P14"])    # 6.8

# Escenario conservador previo, que se conserva para comparar (Thorne 2006)
D_THORNE_CONSERVADOR_nm = 38.0   # modelo de placas paralelas
D_THORNE_PERMISIVO_nm = 64.0     # modelo de poro cilíndrico

_F_NANCE = "Nance et al. 2012, Sci Transl Med 4:149ra119, Fig. 4 y Tabla 2"
_F_THORNE = "Thorne & Nicholson 2006, PNAS 103:5567"
_F_ZETA_POS = ("Berry et al. 2016, RSC Adv 6:41665, Tabla 1 y Fig. 3; "
               "Mastorakos et al. 2016, Small 12:678, Tabla 1 y Fig. 2")


def g_difusion_ecs(d: Diseno, que="transportador", escenario="nance"):
    """¿Puede difundir por el espacio extracelular del cerebro hasta la mielina?

    escenario="nance"   -> criterio de dos variables (tamaño Y superficie).
    escenario="thorne"  -> criterio antiguo, solo tamaño, 38 nm. Se conserva
                           para poder enseñar en cuánto difieren.
    """
    nombre = f"Difusión en espacio extracelular ({que})"

    # --- el fármaco liberado NO se juzga por tamaño.
    #
    #     REESCRITA el 2026-08-17 (tarea F.1, decisión de Jhovan). Antes esta
    #     rama comparaba el diámetro del fármaco (1.683 nm) contra el umbral de
    #     Thorne (38 nm) y devolvía PASA trivialmente: el tamaño no puede ser
    #     limitante con veinte veces de margen. Medía el eje equivocado y su
    #     propia advertencia lo decía.
    #
    #     El eje correcto es si la especie activa puede MOVERSE por el espacio
    #     extracelular sin ayuda. La evidencia validada (filas 55-61) dice que
    #     no, y lo dice con medidas, no con huecos:
    #
    #     1. La especie que exporta SPNS2 es FTY720-FOSFATO, no el fingolimod
    #        neutro. Es un esfingolípido fosforilado, cargado.
    #     2. Foster 2007 declara que el FTY720-P no cruza por sí mismo y que su
    #        transporte depende de portador; es asunción declarada de los
    #        autores, no medida directa.
    #     3. Bucki 2010 (Fig. 1 A-E) MIDE que la gelsolina -- el único portador
    #        candidato propuesto para el LCR -- se une al FTY720-P de forma
    #        débil o nula, mientras que sí interactúa con el S1P. Esto no es un
    #        hueco: descarta activamente al candidato.
    #     4. Foster 2007 (Tabla 3) mide 30-80x menos FTY720-P en LCR que en
    #        tejido/plasma: casi no hay fracción libre en el compartimento
    #        extracelular, coherente con que no viaja solo.
    #
    #     El veredicto es DESCONOCIDA, y conviene ser preciso sobre por qué no
    #     es ninguna de las otras dos:
    #
    #     - NO es PASA. Ese era el resultado del eje equivocado. Que la molécula
    #       quepa por el poro no dice nada sobre si llega.
    #     - NO es FALLA. Se escribió FALLA en una primera versión de esta
    #       reescritura (2026-08-17) y era un SALTO: Bucki descarta a la
    #       gelsolina, no a todos los portadores. Albúmina y apoM son portadores
    #       conocidos del S1P y NO están medidos para FTY720-P, ni a favor ni en
    #       contra. El 1 % de BSA del medio de Hisano 2011 tampoco sirve de
    #       prueba: ese ensayo mide exportación DESDE la célula, no difusión POR
    #       el espacio extracelular. Y Mishima/Kurano 2018 (Biosci Rep 38(5),
    #       doi 10.1042/BSR20181288) mide que el DH-S1P, análogo del S1P, NO
    #       hereda su biología de portador -- el S1P se une a HDL vía apoM y el
    #       DH-S1P no. Ser análogo del S1P no basta para heredarle los
    #       portadores, así que tampoco se puede asumir lo contrario.
    #
    #     Lo que falta para cerrar la compuerta en un sentido u otro: una medida
    #     de unión de FTY720-P a albúmina/apoM, y la concentración de esos
    #     portadores en el ECS cerebral (la BHE excluye albúmina, pero ese
    #     número no está verificado en este proyecto).
    if que != "transportador":
        return Resultado(
            nombre, DESCONOCIDA, None, None, None, None,
            "Foster et al. 2007, JPET 323(2):469, Tabla 3 y Discussion · "
            "Bucki et al. 2010, Am J Physiol Cell Physiol 299(6):C1516, Fig. 1 A-E · "
            "Mishima y Kurano et al. 2018, Biosci Rep 38(5):BSR20181288",
            "la especie activa es FTY720-FOSFATO, un esfingolípido cargado: "
            "Foster 2007 declara que su transporte depende de portador y mide "
            "30-80x menos FTY720-P en LCR que en tejido, y Bucki 2010 mide que "
            "la gelsolina NO lo une, a diferencia del S1P. Pero eso descarta UN "
            "portador, no todos: albúmina y apoM no están medidos para esta "
            "molécula, y Mishima 2018 muestra que un análogo del S1P puede no "
            "heredar sus portadores. Sin esa medida no se puede decir ni que "
            "llega ni que no llega. Tarea F.1")

    if escenario == "thorne":
        return _cmp(nombre, d.diametro_nm, D_THORNE_CONSERVADOR_nm, "nm",
                    _F_THORNE, False)

    z = d.zeta_mV

    # ------------------------------------------------ eje 1: superficie
    if z >= ZETA_ADHESIVO_POSITIVO_mV:
        return Resultado(nombre, FALLA, z, ZETA_ADHESIVO_POSITIVO_mV, "mV",
                         ZETA_ADHESIVO_POSITIVO_mV - z, _F_ZETA_POS,
                         "superficie catiónica adhesiva: con ζ +10.0 mV menos "
                         "del 10 % de la población difundió en el parénquima "
                         "(Berry) y con ζ +35.3 mV quedó inmovilizada "
                         "(Mastorakos), las dos por MPT en cerebro de rata "
                         "ex vivo. SALVEDAD: los dos son polímero/ADN, no "
                         "liposomas, y su ζ está medido en NaCl 10 mM pH 7.0, "
                         "no en aCSF; en aCSF ambos pierden estabilidad "
                         "coloidal, así que adhesión y agregación no están "
                         "separadas experimentalmente")

    if z > 0.0:
        return Resultado(nombre, DESCONOCIDA, z, ZETA_ADHESIVO_POSITIVO_mV, "mV",
                         None, _F_ZETA_POS,
                         f"ζ positivo pero por DEBAJO del mínimo medido "
                         f"(+{ZETA_ADHESIVO_POSITIVO_mV:.1f} mV, Berry 2016). "
                         "Entre 0 y +10 mV no hay ni un dato: Nance solo midió "
                         "de −2.5 a −52 mV, y los dos únicos puntos positivos "
                         "que existen están en +10.0 y +35.3 mV. No se "
                         "extrapola. La física apunta a MÁS adhesión, porque la "
                         "matriz extracelular y las superficies celulares son "
                         "negativas, pero eso es expectativa, no medida")

    if z <= ZETA_NO_DIFUSIVO_mV:
        return Resultado(nombre, FALLA, z, ZETA_NO_DIFUSIVO_mV, "mV",
                         z - ZETA_NO_DIFUSIVO_mV, _F_NANCE,
                         "superficie adhesiva: 0 de 6 formulaciones con ζ más "
                         "negativo que −6 mV difundieron, a cualquier tamaño")

    if z < ZETA_DIFUSIVO_mV:
        return Resultado(nombre, DESCONOCIDA, z, ZETA_DIFUSIVO_mV, "mV",
                         None, _F_NANCE,
                         "zona gris de superficie: entre −4 y −6 mV difundió 1 "
                         "de 2 formulaciones")

    # ------------------------------------------------ eje 2: tamaño
    if d.diametro_nm <= D_SIGILOSO_PASA_nm:
        # OJO: no se usa _cmp aquí. El umbral es INCLUSIVO (Nance demuestra que
        # la de 114 nm penetra), y _cmp exige margen > 0 estricto, con lo que una
        # partícula de exactamente 114 nm saldría FALLA.
        r = Resultado(nombre, PASA, d.diametro_nm, D_SIGILOSO_PASA_nm, "nm",
                      D_SIGILOSO_PASA_nm - d.diametro_nm, _F_NANCE)
        r.advertencia = (
            "sostenido sobre Nance, que mide sobre todo EX VIVO y en corteza "
            "humana de cirugía de epilepsia. Con el criterio antiguo de Thorne "
            f"({D_THORNE_CONSERVADOR_nm:.0f} nm) este diseño "
            f"{'pasaría igual' if d.diametro_nm < D_THORNE_CONSERVADOR_nm else 'FALLARÍA'}")
        return r

    if d.diametro_nm >= D_SIGILOSO_FALLA_nm:
        return Resultado(nombre, FALLA, d.diametro_nm, D_SIGILOSO_FALLA_nm, "nm",
                         D_SIGILOSO_FALLA_nm - d.diametro_nm, _F_NANCE,
                         "las de 200 nm con PEG denso no se dispersaron ni "
                         "ex vivo ni in vivo")

    return Resultado(nombre, DESCONOCIDA, d.diametro_nm, D_SIGILOSO_PASA_nm, "nm",
                     None, _F_NANCE,
                     f"zona gris de tamaño: entre {D_SIGILOSO_PASA_nm:.0f} y "
                     f"{D_SIGILOSO_FALLA_nm:.0f} nm no hay dato COMPARABLE, aun "
                     "con superficie sigilosa. El único punto que cae dentro es "
                     "Curtis 2019 (PS-PEG, 163.2 nm de media en intensidad, "
                     "ζ −6.2 mV, Deff 0.22 µm²/s, razón 18) y NO decide, por dos "
                     f"motivos: es corteza de rata P14, la edad más permeable, y "
                     f"corregido a adulto por el factor {EDAD_FACTOR_P14_A_P70:.1f} "
                     "de McKenna 2021 la razón sube a ~124, que cae entre las que "
                     "difunden (36) y la que no (1600); y su diámetro es media en "
                     "INTENSIDAD, métrica que no es la de Nance")


# -- Tarea B.3 -------------------------------------------------------------
#  Abierta el 2026-08-10 como PASA con salvedad; REABIERTA y devuelta a
#  DESCONOCIDA el 2026-08-12 (decisión de Jhovan). Ficha:
#  verificacion/verificacion_transito_tarea_B_3.md
#
#  La compuerta compara DOS TIEMPOS. El transportador solo sirve si la célula
#  llega a la lesión ANTES de que el liposoma haya soltado el fármaco en sangre.
#
#  POR QUÉ SE CAYÓ EL PASA. La versión de 2026-08-10 usaba los 20 h de Yona
#  2013 como tiempo de tránsito, pero ese número es la semivida del monocito
#  Ly6C+ EN CIRCULACIÓN: mide SALIR DE LA SANGRE, no LLEGAR A LA LESIÓN. Se
#  sale también al bazo, al hígado y a la médula. Tong et al. 2016 sí mide lo
#  que la compuerta pregunta —infusión IV de monocitos cargados con
#  nanopartícula y recuento en cerebro inflamado— y da un PICO a las 48 h, con
#  células ya detectables en el primer punto de muestreo, 24 h.
#
#  Los dos intervalos SE SOLAPAN y por eso no hay decisión posible:
#    tránsito  24–48 h   (Tong 2016; resolución de muestreo 24 h)
#    descarga  > 24 h    (Mao 2014; es COTA INFERIOR, sin techo medido)

T_TRANSITO_PRIMERA_DETECCION_h = 24.0  # Tong 2016, S1 Text: primer muestreo
T_TRANSITO_PICO_h = 48.0               # Tong 2016, S1 Text: "peaked at 48h"
T_LIBERACION_COTA_INFERIOR_h = 24.0    # Mao 2014: a las 24 h retiene >50 %

# Conservado solo como registro de la versión retirada; NO gobierna nada.
T_MEDIO_MONOCITO_h = 20.0   # Yona 2013: vida media EN CIRCULACIÓN del Ly6C+


def g_transito_vs_liberacion(d: Diseno):
    """¿Llega la célula a la lesión antes de que el fármaco se suelte? (B.3)

    Tiempo de tránsito — Tong et al. 2016, S1 Text, textual: "Brain tissues
    were collected at 1, 2, 3 and 7 days following MDM transfers into LPS ICI
    treated mice. The number of recruited donor-derived cells peaked at 48h and
    decreased afterwards". Las células ya se detectan en el primer punto de
    muestreo (24 h), así que el tránsito real cae entre 24 y 48 h.

    Tiempo de descarga — Mao et al. 2014, textual: "The release over a time
    period of 12 h was <40%, while at 24 h the liposomes still retained over
    50% FTY720 in both media" (PBS y 10 % de suero, 37 °C). Es una COTA
    INFERIOR de la semivida de liberación: >24 h, sin techo.

    24–48 h contra >24 h sin techo: los intervalos se solapan. La compuerta NO
    se puede decidir con los datos publicados y devuelve DESCONOCIDA.
    """
    return Resultado(
        "Tránsito del monocito frente a cinética de liberación",
        DESCONOCIDA,
        None, None, "h", None,
        "Tong et al. 2016, PLoS ONE 11:e0154022 · Mao et al. 2014, Nanomedicine 10:393",
        motivo=(
            f"tránsito {T_TRANSITO_PRIMERA_DETECCION_h:.0f}–{T_TRANSITO_PICO_h:.0f} h "
            f"frente a descarga >{T_LIBERACION_COTA_INFERIOR_h:.0f} h sin techo "
            "medido: los intervalos se solapan  [tarea B.3]"),
        advertencia=(
            f"tránsito {T_TRANSITO_PRIMERA_DETECCION_h:.0f}–{T_TRANSITO_PICO_h:.0f} h "
            f"(Tong 2016) frente a descarga >{T_LIBERACION_COTA_INFERIOR_h:.0f} h sin "
            "techo medido (Mao 2014): los intervalos SE SOLAPAN y no hay decisión. "
            "El tránsito de Tong es LPS intracraneal en ratón, no EAE, y el "
            "cargamento es un SPION de 38-40 nm, no un liposoma de 700 nm; la "
            "medida equivalente en EAE a resolución de horas no existe. La "
            "descarga de Mao es IN VITRO en PBS y suero, con un liposoma de "
            "157 nm. Los 20 h de Yona 2013 que gobernaban esta compuerta hasta "
            "el 2026-08-12 miden SALIDA DE CIRCULACIÓN, no llegada a la lesión, "
            "y quedan RETIRADOS del veredicto"))


def g_salida_farmaco(d: Diseno):
    """¿Puede el fármaco salir de la célula que lo transportó? (tarea B.5)

    Implementada el 2026-08-10. Antes devolvía DESCONOCIDA por falta de fuente
    primaria; ya se leyó.

    Hisano et al. 2011 demuestra la maquinaria completa PARA ESTE FÁRMACO:
      · el fingolimod es un PROFÁRMACO; la esfingosina quinasa lo fosforila
        dentro de la célula a FTY720-fosfato, que es la forma activa;
      · el transportador SPNS2 EXPORTA el FTY720-P fuera de la célula, por la
        misma vía que la S1P.
    Textual del resumen: "human SPNS2 can transport several S1P analogues,
    including FTY720-P... This is the first identification of an FTY720-P
    transporter in cells".

    Es una compuerta DEPENDIENTE DEL FÁRMACO, no del diseño de la partícula.
    Para otro fármaco habría que rehacerla.
    """
    return Resultado("Salida del fármaco de la célula transportadora", PASA,
                     fuente="Hisano et al. 2011, J Biol Chem 286:1758",
                     motivo="SphK fosforila dentro, SPNS2 exporta fuera",
                     advertencia=(
                         "demostrado en células CHO TRANSFECTADAS. NO se ha "
                         "demostrado que un macrófago que fagocitó un liposoma "
                         "exporte el fármaco en la lesión. El salto es una "
                         "INFERENCIA. Y SEGUNDA SALVEDAD, del 2026-08-13: el "
                         "medio en el que se demostró la exportación es F-12 con "
                         "1 % de BSA, 10 mM de glicerofosfato sódico, 5 mM de "
                         "fluoruro sódico y 1 mM de semicarbazida, es decir, con "
                         "portador proteico y con la degradación bloqueada. El "
                         "espacio extracelular de una lesión no es ese medio"))


def g_captacion_fagocitica(d: Diseno):
    """¿Lo capta un macrófago con eficacia suficiente?

    Compuerta INVERSA a las anteriores: aquí grande es mejor.
    Anclaje: Muselman et al. 2026 comparó 150, 550 y 700 nm en EAE. Los de
    150 nm captaban mal; los de 550 y 700 nm funcionaron. El umbral está entre
    150 y 550 nm y NO está determinado con precisión, así que se declara una
    ZONA GRIS y en ella la compuerta devuelve DESCONOCIDA.

    CORREGIDO el 2026-08-13 (auditoría de Jhovan): se cambió el umbral alto de
    550 a 500 nm, creyendo que el artículo comparaba 150/500/700 nm.
    REVERTIDO el 2026-08-17 (C10, lectura completa del artículo vía PMC): el
    diseño experimental real (Métodos y pie de Figura 2) fabricó y marcó con
    FITC tres tamaños, 150/550/700 nm — ese es el dato que define el umbral.
    El "500" que motivó la corrección de agosto solo aparece una vez, en la
    prosa de Resultados, resumiendo cuáles liposomas targetearon mejor; no hay
    un cuarto grupo de 500 nm descrito en Métodos. Se lee como errata de
    redacción de los propios autores, no como un grupo experimental real. Sin
    verificar contra las gráficas del PDF (lectura vía texto extraído de PMC).
    """
    bajo, alto = 150.0, 550.0
    fuente = "Muselman et al. 2026, Front Immunol 16:1657131"
    if d.diametro_nm <= bajo:
        return Resultado("Captación fagocítica", FALLA, d.diametro_nm, alto, "nm",
                         d.diametro_nm - alto, fuente,
                         "por debajo o igual a 150 nm la captación fue pobre")
    if d.diametro_nm >= alto:
        return Resultado("Captación fagocítica", PASA, d.diametro_nm, alto, "nm",
                         d.diametro_nm - alto, fuente)
    return Resultado("Captación fagocítica", DESCONOCIDA, d.diametro_nm, alto, "nm",
                     None, fuente,
                     "zona gris entre 150 y 550 nm: no hay dato que la resuelva")


# =============================================================================
#  COMPUERTAS SIN DATO  (devuelven DESCONOCIDA a propósito)
#  Cada una corresponde a una tarea abierta del cronograma v4.
# =============================================================================

def _sin_dato(nombre, motivo, tarea):
    def f(d: Diseno):
        return Resultado(nombre, DESCONOCIDA, motivo=f"{motivo}  [tarea {tarea}]")
    return f


g_union_glicocalix = _sin_dato(
    "Unión al glicocálix e internalización",
    "no se sabe si una partícula unida a las fibras externas alcanza la "
    "membrana. Cheng 2016 separa los dos sucesos y mide que la unión NO implica "
    "llegada, pero es endotelio de almohadilla grasa de rata en cultivo "
    "estático y con partícula ANIÓNICA, y la ruta C necesita una catiónica. "
    "Además el trayecto real es MÁS largo de lo que se creía: 540-726 nm de "
    "glicocálix cerebral (Shi y Larsen 2025) frente a los 150-400 nm que se "
    "venían usando", "C.2")

g_acceso_receptor = _sin_dato(
    "Acceso al receptor bajo el glicocálix",
    "en endotelio CEREBRAL el glicocálix mide 540 nm (Shi 2025) o 726 nm "
    "(Larsen 2025), no los 150-400 nm de rana y hámster que se venían usando, "
    "y cubre el 93 % de la superficie, no el 75 % de cultivo estático. La vía "
    "de escape que se había propuesto para D.3 (huecos sin cubrir) queda "
    "reducida al 7 %. La contradicción con la transcitosis por receptor, que "
    "SÍ funciona en la BHE, se AGRAVA en vez de resolverse",
    "D.3")

g_transcitosis = _sin_dato(
    "Transcitosis completa",
    "transporte activo dependiente de ATP, fuera del alcance de la termodinámica de equilibrio",
    "fuera de alcance")

# Añadida el 2026-08-13, decisión de Jhovan. Hasta hoy el modelo preguntaba si
# el transportador LLEGA y nunca si entrega FÁRMACO SUFICIENTE, con lo que un
# diseño que llegase con una sola molécula habría salido igual de bien que uno
# que llegase con miles. Era un hueco INVISIBLE, y la lógica de este proyecto es
# que los huecos se declaren. Va en las CUATRO rutas porque la pregunta no
# depende del mecanismo de entrada.
# ACTUALIZADA el 2026-08-17. La versión del 2026-08-13 decía que faltaban TRES
# números y que "ninguno está en el proyecto". Eso ya NO es cierto: dos de los
# tres entraron, y conviene que la compuerta lo diga, porque un hueco declarado
# más grande de lo que es también desinforma.
#
#   (1) CARGA DEL LIPOSOMA — RESUELTO. Mouzoura et al. 2025 (Int J Nanomedicine
#       20:239, doi 10.2147/IJN.S494512) da relación molar fármaco:lípido 1:8
#       con eficiencia de carga 94-97.2 %, con FTY720 directo. No es la
#       eficiencia de encapsulación de Mao 2014, que era lo que había antes y
#       no servía.
#   (2) UMBRAL EN PARÉNQUIMA — RESUELTO. Foster et al. 2007 (JPET 323(2):469,
#       Tabla 3) mide 398 ± 186 ng/g de FTY720-P en cerebro a una dosis que FUE
#       terapéuticamente eficaz en EAE (0.3 mg/kg, días 11-33). Eso da la
#       concentración diana contra la que comparar.
#   (3) CARGA DEL DENDRÍMERO CON FINGOLIMOD — sigue sin existir. El 1:1 de
#       Devarakonda 2004 es con nifedipino y es estequiometría de complejo
#       inferida de una pendiente, no carga útil medida. Congelado: el
#       dendrímero quedó fuera del foco de simulación el 2026-08-17.
#
#   (4) FRACCIÓN DE DOSIS ENTREGADA AL PARÉNQUIMA — RESUELTO el 2026-08-17i.
#       Era el último eslabón. Johnsen et al. 2019 (J Control Release 295:237,
#       Fig. 5D) es el único estudio encontrado que reúne las TRES condiciones
#       a la vez, y por eso se usa este y no otro:
#         - dosis INTRAVENOSA (no perfusión in situ, no inyección intracerebral),
#         - la sangre se separa DOS veces: perfusión transcardíaca Y ADEMÁS
#           depleción capilar con dextrano (Triguero 1990), que aparta la
#           fracción vascular del parénquima,
#         - marcador ELEMENTAL por ICP-MS, no fluoróforo ni radiomarca que se
#           despegue del vehículo.
#       Ratón Balb/c, n = 7-8, 2.5 h. Liposoma de control "stealth": 132.6 ±
#       0.5 nm, ζ = −17.1 ± 2.5 mV, PDI 0.032 — PEGilado y SIN ligando, es
#       decir, la RUTA PASIVA, y cae dentro del rango del catálogo.
#
# POR QUÉ ESTA COMPUERTA YA NO ES DESCONOCIDA
#
# Con (2) y (4) la pregunta se vuelve aritmética: si el parénquima recibe
# 0.0145 %ID/g por vía pasiva, ¿qué dosis hay que inyectar para llegar a los
# 398 ng/g que Foster midió a una dosis EFICAZ de 0.3 mg/kg? Salen ~110-137
# mg/kg, entre dos y tres órdenes de magnitud por encima. La compuerta FALLA.
#
# El resultado se expresa como FACTOR DE EXCESO DE DOSIS (cuántas veces la
# dosis eficaz de Foster haría falta) y no como concentración, porque el factor
# es adimensional y sobrevive a los supuestos de masa corporal.
#
# LA SALVEDAD QUE NO SE PUEDE PERDER, y que va dentro del Resultado:
# el 0.0145 %ID/g es de OXALIPLATINO encapsulado, no de fingolimod. Es
# transferencia de parámetro entre contextos, la misma clase de salto ya
# fichado en la salvedad B.3(a). Se acepta aquí porque el orden de magnitud
# (>100x) es demasiado grande para que lo revierta un cambio de fármaco, pero
# NO se presenta como medida de fingolimod. Además: ratón sano, no EAE (sin
# factor de apertura de la BHE), un solo punto temporal, y los propios autores
# admiten que su señal parenquimal podría ser arrastre del método de depleción.
#
# Lo que esto NO dice: no dice que el liposoma no llegue nunca al cerebro. Dice
# que por la ruta pasiva no alcanza la concentración que Foster asoció a
# eficacia, ni de lejos, y que encapsular un fármaco tan lipófilo lo empeora
# frente a administrarlo libre.
#
# El DENDRÍMERO sigue sin su carga con fingolimod (el 1:1 de Devarakonda es con
# nifedipino), así que para él la compuerta se queda DESCONOCIDA: no se le puede
# aplicar la aritmética sin el numerador.

# Fracción de la dosis que llega al parénquima, por gramo de tejido.
#
# RUTA PASIVA — grupo mPEG (sin ligando) de Johnsen 2019.
# NO se usa la lectura del eje de la Fig. 5D (≈0.0145 %ID/g), porque ese número
# no está tabulado ni escrito en ninguna parte del artículo y dependería de leer
# una barra a ojo. Se usa una COTA SUPERIOR DERIVADA DEL TEXTO, que es
# verificable palabra por palabra:
#   el cuerpo dice que el grupo 0.6e3 Abs/µm² alcanzó "approximately 0.04 %ID/g"
#   y que eso fue "more than a doubling of the platinum concentration ...
#   compared to the mPEG control".
#   Si 0.6e3 ≈ 0.04 y es MÁS del doble del mPEG  =>  mPEG < 0.02 %ID/g.
# Se toma 0.02: es cota superior (favorece que la compuerta pase) y no cuelga de
# ninguna lectura visual. El factor sale 295x en vez de 407x: menos espectacular
# y mucho más defendible.
_F_PARENQUIMA_PASIVA = 2.0e-4           # <0.02 %ID/g, COTA derivada del texto
# Mejor densidad de anticuerpo del mismo estudio (0.6e3 Abs/µm²), ruta activa.
# SE USA EL VALOR DEL TEXTO (0.04), NO EL DE LA FIGURA (0.036). El cuerpo del
# artículo dice "approximately 0.04%ID/g" y la lectura del eje de la Fig. 5D da
# 0.036. Ante la discrepancia se toma el MÁS ALTO, que es el que más favorece
# que la compuerta pase, por la misma regla de cota superior de más abajo.
_F_PARENQUIMA_ACTIVA = 4.0e-4           # 0.04 %ID/g  (texto; figura da 0.036)
# Foster 2007: concentración medida a una dosis que fue eficaz en EAE.
# VERIFICADO EN FUENTE PRIMARIA el 2026-08-17i (PDF leído, tabla de niveles
# valle). Al verificarlo aparecieron DOS cosas que ni el handoff ni este archivo
# decían, y las dos importan:
#   (a) son NIVELES VALLE ("Trough Level Brain"), no pico ni media. Es el mínimo
#       sostenido entre dosis. Como umbral terapéutico es defendible, pero no es
#       "la concentración que alcanzó el fármaco".
#   (b) son valle tras dosificación DIARIA de los días 11 a 33 = 23 días. NO es
#       una dosis única. Comparar eso contra la entrega de un solo pinchazo de
#       liposoma (Johnsen mide a 2.5 h de UNA dosis) penaliza al liposoma de
#       forma artificial.
# Por eso la compuerta reporta una BANDA y no un número (decisión de Jhovan,
# 2026-08-17i): el extremo alto compara dosis única contra dosis diaria, el
# extremo bajo compara dosis acumulada contra dosis acumulada. La verdad está
# entre los dos y para estrecharla haría falta la vida media del liposoma en
# parénquima, que no está publicada.
_UMBRAL_PARENQUIMA_NG_G = 398.0         # ± 186 ng/g FTY720-P (neutro: 313 ± 149)
_DOSIS_EFICAZ_MG_KG = 0.3               # mg/kg/DÍA, RATA DA con EAE, n=6 cerebro
# Factor de acumulación en CEREBRO por dosificación repetida. MEDIDO por Foster,
# no supuesto: "The accumulation factor (i.e., ratio between the tissue
# concentration after multiple versus single dosing) was 3.2 for brain and 3.7
# for spinal cord at 24 h after the seventh dose compared with the single dose"
# (Results, subsección de autorradiografía — NO Discussion). CORRIGE un error de Kiel del 2026-08-17i, que había asumido
# acumulación LINEAL sobre los 23 días de dosificación (×23) sin buscar si el
# artículo la daba medida. La daba. El supuesto lineal inflaba el margen 7 veces
# a favor del liposoma, y además era incoherente con el propio dato que usaba:
# si no hubiera eliminación entre dosis, no existiría un "nivel valle".
#
# TERCERA TRANSFERENCIA ENTRE CONTEXTOS, declarada (detectada al auditar por
# tercera vez, 2026-08-17i): este 3.2 NO sale del experimento del umbral. Sale
# del estudio de autorradiografía del mismo artículo, que es OTRO montaje:
#     umbral 398 ng/g -> rata DA hembra, CON EAE, 0.3 mg/kg/d, 23 dosis, LC-MS/MS
#     factor 3.2      -> rata LE/CR WIGA macho SANA, 7.5 mg/kg/d (25x más),
#                        7 dosis, radiactividad 14C total
# Se usaron juntos porque están en el mismo artículo. No son el mismo experimento.
# Dirección del sesgo DESCONOCIDA: a 25x la dosis puede haber saturación de
# eliminación (inflaría el 3.2), pero con 7 dosis en vez de 23 puede no haberse
# alcanzado el estado estacionario (lo desinflaría). Sin la vida media no se
# acota. Es el supuesto más frágil de la cadena junto al salto de fármaco.
_ACUMULACION_CEREBRO = 3.2              # medido a la 7ª dosis, no a los 23 días
# Supuesto DECLARADO, no medido, y no está en ninguno de los dos artículos.
# Es RATÓN porque la fracción entregada es de ratón (Johnsen); el umbral que se
# compara es de rata. Esa mezcla va declarada en la advertencia del Resultado.
_MASA_RATON_KG = 0.0225                 # 20-25 g


def _g_carga_util(d: Diseno, fraccion, procedencia):
    """G.2 — ¿la carga útil entregada alcanza el umbral terapéutico?

    DESCONOCIDA. REABIERTA el 2026-08-17l tras una auditoría independiente.

    HISTORIA, porque este archivo no oculta sus retractaciones:
      · 2026-08-13 se añade la compuerta, DESCONOCIDA por falta de tres números.
      · 2026-08-17i se cierra con FALLA (46x-147x) al entrar Johnsen 2019.
      · 2026-08-17j se recalcula a 2.2x-3.9x metiendo un cociente de eficiencias.
      · 2026-08-17k una auditoría independiente tumba la cifra con 16 hallazgos.
      · 2026-08-17l se RETIRA la cadena entera y la compuerta vuelve a DESCONOCIDA.

    POR QUÉ SE RETIRA (los tres motivos que la matan, en orden de gravedad):

    (1) EL COCIENTE NO CANCELA EL UMBRAL. Un cociente de %ID/g cancela la DOSIS
        ADMINISTRADA, no la concentración cerebral eficaz. Solo se cancelaría si
        los dos fármacos exigieran la misma concentración en cerebro y toleraran
        la misma dosis sistémica, y eso no está establecido para ninguno de los
        dos. Lo que la comparación mide legítimamente es «por unidad de dosis
        administrada, qué vehículo deposita más masa por gramo de cerebro», que
        NO autoriza a afirmar nada sobre eficacia. Este es un fallo LÓGICO y no
        se arregla con mejores números.

    (2) LA BANDA DEL 17j ERA UN ARTEFACTO ARITMÉTICO. El factor 756 aparecía en
        numerador y denominador y se cancelaba: 756/A con A = 3.2x(756/224) da
        224/3.2 = 70.0 ng/g exactamente, o sea la fila de 7 días de la Tabla 2
        dividida por 3.2. Los datos de 23 dosis no aportaban NADA. La cota
        inferior de la banda era además un error de contabilidad, por mezclar
        numerador de la Tabla 3 con denominador de la Tabla 2.

    (3) LOS DOS TRAMOS DE ACUMULACIÓN SON INCOMPATIBLES. Componer el 3.2 (tramo
        1->7 dosis, SUBLINEAL) con 3.375 (tramo 7->23, LINEAL) invierte la
        cinética: los factores de acumulación desaceleran hacia el estado
        estacionario, no aceleran. Y los propios datos de Foster (Tabla 2, ng/g
        por dosis administrada a 0.3 mg/kg: 32.0 a 7 d, 36.1 a 12 d, 32.9 a
        23 d) describen acumulación LINEAL en el número de dosis, o sea A ~ 23,
        no 10.8. Con A = 23 la ruta activa quedaba en 1.0x-1.8x: paridad dentro
        del ruido, sin poder de decisión.

    SESGO SISTEMÁTICO, que es lo más serio: los cuatro supuestos de mayor
    apalancamiento numérico (factor de acumulación, ancla temporal de la
    enfermedad, masa mínima de la rata, elección del fármaco de referencia)
    apuntaban TODOS en la dirección que sostenía la conclusión, y dos de ellos
    no estaban declarados. No era un error aislado.

    QUÉ HARÍA FALTA PARA CERRARLA (y por qué no se inventa mientras tanto):
      (a) concentración cerebral EFICAZ de fingolimod, como umbral propio, no
          prestada de un nivel valle de otro régimen;
      (b) fracción de dosis entregada al parénquima por un liposoma cargado con
          FINGOLIMOD tras administración IV, medida con separación de la
          fracción vascular;
      (c) idealmente ambas en el mismo modelo animal y estado de enfermedad.
    Mientras falte cualquiera de las tres, la salida es DESCONOCIDA. Un hueco
    declarado NO es un aprobado: ver la cabecera del archivo, NO EVALUABLE.

    Las constantes de la cadena se conservan abajo SOLO como registro de lo que
    se retiró. NINGUNA se usa ya para decidir.
    """
    return Resultado(
        "Carga útil suficiente", DESCONOCIDA,
        fuente="retirada 2026-08-17l tras auditoría independiente",
        motivo="RETRACTACIÓN. La cadena Foster 2007 + Johnsen 2019 que cerraba "
               "esta compuerta con FALLA se retira por tres motivos: (1) el "
               "cociente de %ID/g cancela la DOSIS ADMINISTRADA, no el umbral "
               "de concentración eficaz, así que no autoriza una afirmación de "
               "eficacia; (2) la banda del 17j era un artefacto: el 756 se "
               "cancelaba algebraicamente y la cota inferior mezclaba Tabla 2 "
               "con Tabla 3; (3) los dos tramos de acumulación (3.2 sublineal "
               "y 3.375 lineal) describen cinéticas incompatibles, y los "
               "propios datos de Foster dan acumulación ~23, con la que el "
               "resultado quedaba en paridad dentro del ruido. Falta el umbral "
               "eficaz de FINGOLIMOD y la fracción entregada de un liposoma "
               "cargado con FINGOLIMOD: sin esos dos números medidos en el "
               "mismo contexto, la pregunta no es decidible  [tarea G.2]",
        advertencia="NO es un aprobado. Es un hueco declarado: la ruta queda NO "
                    "EVALUABLE, no NO EXCLUIDA. Para la ruta D · mediada por "
                    "receptor el veredicto NO depende de esta compuerta: esa "
                    "ruta queda EXCLUIDA por la compuerta de transcitosis por "
                    "TfR, que es un argumento de mecanismo y no de dosis")


def _advertencia_retirada_17l():
    """Texto de la salvedad de la versión RETIRADA de G.2.

    Se conserva literal para que quede registro de qué se afirmaba y con qué
    salvedades, y para que la retractación sea auditable. NO se usa.
    """
    fraccion, procedencia, factor, factor_acum = 0.0, "", 0.0, 0.0
    return (
        f"BANDA {factor_acum:.0f}×-{factor:.0f}×: el extremo bajo descuenta el "
        f"factor de acumulación por dosis repetida MEDIDO por Foster en cerebro "
        f"({_ACUMULACION_CEREBRO}×, Results/autorradiografía — OJO: medido en OTRO montaje del "
        "mismo artículo, rata sana a 25× la dosis y 7 dosis en vez de 23, con "
        "dirección de sesgo desconocida); el alto compara dosis única "
        "contra dosis diaria sin descontar nada. El veredicto se decide con el "
        "extremo BAJO. El umbral de Foster es un nivel VALLE tras 23 días de "
        "dosificación, no un pico ni una dosis única. ASIMETRÍA DE MÉTODO "
        "declarada: el umbral es de CEREBRO EN BLOQUE sin perfundir (Foster "
        "homogeneiza medio cerebro directo, nivel 1-2 de Yokel 2020) mientras la "
        "fracción entregada es de parénquima con doble separación de sangre; "
        "numéricamente es despreciable (plasma 7.35 ng/ml × ~1.4 % ≈ 0.1 ng/g "
        "frente a 398) y va en contra del liposoma, no a su favor. "
        f"Fracción entregada {fraccion*100:.4f} %ID/g ({procedencia}). "
        "OXALIPLATINO encapsulado, no de fingolimod: transferencia de parámetro "
        "entre contextos, misma clase de salto que la salvedad B.3(a). OJO: con la "
        "banda corregida el margen del extremo bajo YA NO es de dos órdenes de "
        "magnitud, así que este salto de fármaco SÍ puede discutir el resultado y "
        "deja de ser una salvedad menor. NO es una medida de fingolimod. Además: ratón sano (sin "
        "factor de apertura de la BHE por EAE), un solo punto temporal (2.5 h), "
        "masa de ratón supuesta (22.5 g, no está en los artículos), y los propios "
        "autores admiten que su señal parenquimal podría ser arrastre del método "
        "de depleción capilar. MEZCLA DE ESPECIES, declarada: el umbral y la "
        "dosis eficaz son de RATA DA con EAE (Foster), la fracción entregada es "
        "de RATÓN Balb/c SANO (Johnsen); el cociente asume que la dosis eficaz "
        "en mg/kg se transfiere entre especies, que el escalado alométrico NO "
        "garantiza. MEZCLA DE ESPECIE QUÍMICA, declarada: el umbral es de "
        "FTY720-P (metabolito fosforilado) y el liposoma entregaría FTY720; con "
        "el neutro de Foster (313 ng/g) el factor bajaría ~21 %. Las dos "
        "fracciones son COTAS SUPERIORES tomadas del TEXTO de Johnsen, no "
        "lecturas de figura. El veredicto FALLA aguanta las cuatro cosas a la vez (peor "
        "caso combinado: 87×); la CIFRA exacta no debe citarse como precisa  "
        "[tarea G.2 — RETIRADA el 2026-08-17l]")


# =============================================================================
#  TRANSCITOSIS POR TfR DE ALTA AFINIDAD  (ruta D, añadida el 2026-08-17l)
#
#  Esta compuerta sustituye a G.2 como fundamento del veredicto de la ruta D, y
#  la diferencia de fondo es que NO es un argumento de dosis sino de MECANISMO.
#  No pregunta cuánto llega: pregunta si la vía por la que se supone que llega
#  existe. Por eso no depende de ninguna transferencia de parámetro entre
#  especies, fármacos ni métodos analíticos, que es lo que tumbó a G.2.
#
#  LO QUE DICE LA FUENTE, verificado en el cuerpo del artículo el 2026-08-17l
#  (PDF de Johnsen 2019 leído entero, no el resumen):
#
#   (a) «the utility of this antibody for transport across the BBB is disputed»
#       y RI7 «resembles that of the anti-rat TfR antibody, OX26, which is known
#       to be confined to the brain capillaries [30,31], unless its affinity is
#       reduced [32,33]».
#   (b) «The transport capacity of the RI7 antibody also reflects the accepted
#       pathway of iron uptake into the brain, which DOES NOT INCLUDE
#       TRANSCYTOSIS OF THE TfR from the luminal to the abluminal surface
#       [9,10,34]. Accordingly, we expected only little transport of the AuNPs
#       or liposomes into the brain parenchyma.»
#   (c) medido: «the brain uptake of RI7 nanoparticles was < 1/30 of what was
#       previously shown for the antibody itself».
#
#  EL PESO ESTÁ EN (a) Y (b), NO EN (c), y esto es deliberado. El «< 1/30» NO es
#  un control interno: compara la medida de Johnsen contra 1.6 %ID/g tomado de
#  la referencia [23], que es OTRO estudio, y el propio Johnsen dice en el mismo
#  párrafo que ese valor está disputado («Others, however, failed to prove the
#  increased transport across the BBB»). Usarlo como fundamento principal sería
#  repetir exactamente el error de G.2: apoyar un veredicto en una transferencia
#  entre artículos. Se cita como consistente, no como prueba.
#
#  TAMPOCO se apoya en la frase del «proceso pasivo», y también es deliberado.
#  Johnsen escribe «We SPECULATE that the measured quantities ... is a result of
#  a passive process», y antes reconoce que «Whether this is an active process
#  mediated by the TfR CANNOT BE DECIPHERED FROM THESE DATA». Un veredicto no
#  puede montarse sobre una frase en la que los autores declaran que sus datos
#  no deciden. Se registra en la advertencia, no en el motivo.
#
#  LA CONDICIONAL ES PARTE DEL RESULTADO, no una salvedad decorativa: el propio
#  texto dice «unless its affinity is reduced». La afinidad del ligando es una
#  palanca de DISEÑO. Esta compuerta excluye el diseño de alta afinidad, no la
#  idea de dirigirse al TfR. Si algún día el catálogo incorpora ligandos de
#  afinidad reducida, la compuerta debe reevaluarse, no heredarse.
# =============================================================================

def g_transcitosis_tfr(d: Diseno):
    """Ruta D — ¿el TfR de alta afinidad transcita de luminal a abluminal?"""
    return Resultado(
        "Transcitosis por TfR de alta afinidad", FALLA,
        fuente="Johnsen 2019 (J Control Release 295:237) § Discusión, y las "
               "fuentes primarias que cita: vía del hierro [9,10,34], "
               "confinamiento capilar de OX26 [30,31], rescate por afinidad "
               "reducida [32,33]",
        motivo="la vía aceptada de captación de hierro al cerebro NO incluye "
               "transcitosis del TfR de la superficie luminal a la abluminal; "
               "el anticuerpo de alta afinidad queda confinado al endotelio "
               "capilar, como OX26. Es un argumento de MECANISMO del receptor, "
               "no de dosis entregada: no depende de ninguna transferencia de "
               "parámetro entre especies, fármacos ni métodos  [tarea D.4]",
        advertencia="CONDICIONAL A LA AFINIDAD DEL LIGANDO: el propio texto dice "
                    "'unless its affinity is reduced' [32,33]. Esta compuerta "
                    "excluye el diseño de ALTA afinidad (RI7, OX26), no la "
                    "estrategia de dirigirse al TfR en general; con un ligando "
                    "de afinidad reducida debe REEVALUARSE, nunca heredarse. "
                    "NO se apoya en el '< 1/30', que compara contra 1.6 %ID/g "
                    "de OTRO estudio [23] y que el propio Johnsen declara "
                    "disputado; se cita solo como consistente. TAMPOCO se apoya "
                    "en la frase del 'proceso pasivo', que Johnsen marca como "
                    "especulación y sobre la que escribe que 'whether this is "
                    "an active process mediated by the TfR cannot be deciphered "
                    "from these data'. Un solo artículo como fuente secundaria; "
                    "las primarias [9,10,34,30-33] NO se han leído enteras "
                    "todavía y verificarlas es la primera tarea pendiente")


# La fracción entregada NO es la misma en todas las rutas, y usar la pasiva en
# todas inflaría el FALLA. Regla: dato propio donde lo hay, cota superior donde
# no. Con eso, un FALLA lo es incluso en el mejor caso documentado.
#   A pasiva   -> grupo mPEG de Johnsen (sin ligando). Es SU dato.
#   D receptor -> grupo anti-TfR de Johnsen. Es SU dato, y es 2.5x más alto.
#   B y C      -> Johnsen no midió macrófago ni superficie catiónica. Se les
#                 presta la cota superior (la del ligando), que las favorece.
#
# RETIRADAS el 2026-08-17l: las tres devuelven ya DESCONOCIDA, porque
# `_g_carga_util` ignora la fracción. Se conservan los tres nombres para no
# romper el resto del código y para que el registro de qué fracción se prestaba
# a cada ruta no se pierda.
def g_carga_util_pasiva(d: Diseno):
    return _g_carga_util(d, _F_PARENQUIMA_PASIVA,
                         "COTA del texto de Johnsen: el grupo 0.6e3 ≈0.04 %ID/g "
                         "es 'more than a doubling' del mPEG, luego mPEG <0.02. "
                         "Sin lectura de figura")


def g_carga_util_cota(d: Diseno):
    return _g_carga_util(d, _F_PARENQUIMA_ACTIVA,
                         "COTA SUPERIOR: Johnsen anti-TfR 0.6e3 Abs/µm², la "
                         "entrega más alta que midió. Esta ruta no tiene dato "
                         "propio, así que se le presta la más favorable")


# Alias para el código y las pruebas que piden "la" compuerta sin más contexto:
# se queda con la versión conservadora (cota superior), nunca con la que más
# ayuda a fallar.
g_carga_util = g_carga_util_cota


# =============================================================================
#  RUTAS
# =============================================================================

# g_carga_util va la ÚLTIMA de cada ruta, después de la difusión: primero se
# pregunta si el fármaco llega al sitio y solo entonces si llega bastante.
RUTAS = {
    "A · pasiva (adhesión y envolvimiento)": [
        g_transportador_fabricable,
        g_glicocalix_tamiz,
        g_envolvimiento,
        g_caveola,
        g_transcitosis,
        lambda d: g_difusion_ecs(d, "transportador"),
        g_carga_util_pasiva,          # única ruta con dato propio sin ligando
    ],
    "B · celular (macrófago de Troya)": [
        g_transportador_fabricable,
        g_captacion_fagocitica,
        g_transito_vs_liberacion,
        g_salida_farmaco,
        lambda d: g_difusion_ecs(d, "fármaco liberado"),
        g_carga_util,
    ],
    "C · adsortiva (carga positiva)": [
        g_transportador_fabricable,
        g_union_glicocalix,
        g_transcitosis,
        lambda d: g_difusion_ecs(d, "transportador"),
        g_carga_util,
    ],
    # g_transcitosis_tfr va ANTES de la difusión y de la carga útil: si el
    # mecanismo de entrada no existe, lo que pase después es irrelevante. Es la
    # compuerta que EXCLUYE esta ruta desde el 2026-08-17l, y sustituye a G.2
    # como fundamento del veredicto.
    "D · mediada por receptor": [
        g_transportador_fabricable,
        g_acceso_receptor,
        g_transcitosis_tfr,
        g_transcitosis,
        lambda d: g_difusion_ecs(d, "transportador"),
        g_carga_util,
    ],
}


def evaluar_ruta(diseno: Diseno, compuertas):
    """Evalúa una ruta completa. Devuelve (veredicto, lista de resultados)."""
    res = [c(diseno) for c in compuertas]
    if any(r.estado == FALLA for r in res):
        return "EXCLUIDA", res
    if any(r.estado == DESCONOCIDA for r in res):
        return "NO EVALUABLE", res
    return "NO EXCLUIDA", res


def evaluar(diseno: Diseno):
    """Evalúa el diseño en todas las rutas."""
    return {nombre: evaluar_ruta(diseno, comps) for nombre, comps in RUTAS.items()}


# =============================================================================
#  VALIDACIÓN CONTRA DATOS PUBLICADOS
#  Aquí es donde el simulador deja de ser coherente consigo mismo y pasa a estar
#  contrastado. Cada prueba enfrenta una compuerta a un experimento real.
# =============================================================================

def validar_contra_experimentos(verbose=True):
    """Dos bloques que NO valen lo mismo. Leer el aviso de abajo antes de citarlo.

    CORRECCIÓN de una sobreafirmación anterior: las pruebas de CONSISTENCIA
    reproducen los umbrales con los que se construyeron las propias compuertas.
    Son circulares por construcción. No son validación. La validación real del
    proyecto está en verificacion/resultados_validacion.md, contra literatura
    que no se usó para calibrar nada.
    """
    ok = []

    def chequeo(nombre, cond, detalle=""):
        ok.append(bool(cond))
        if verbose:
            print(f"  [{'OK ' if cond else 'FALLA'}] {nombre}{'  ' + detalle if detalle else ''}")

    def titulo(t, aviso=""):
        if verbose:
            print("=" * 78)
            print(f" {t}")
            if aviso:
                print(f" {aviso}")
            print("=" * 78)

    titulo("BLOQUE 1 · CONSISTENCIA",
           "Reproduce la calibración. CIRCULAR por construcción: NO es validación.")

    # -- Muselman 2026: 150 nm se captan mal, 700 nm funcionan.
    r150 = g_captacion_fagocitica(Diseno("test", 150.0, 0.0))
    r700 = g_captacion_fagocitica(Diseno("test", 700.0, 0.0))
    chequeo("C1 captación fagocítica a 150 nm debe FALLAR", r150.estado == FALLA,
            f"da {r150.estado} (Muselman: captación pobre)")
    chequeo("C2 captación fagocítica a 700 nm debe PASAR", r700.estado == PASA,
            f"da {r700.estado} (Muselman: funcionó)")

    # -- Nance 2012, Tabla 2: la de 114 nm con ζ −2.5 mV penetra.
    r114 = g_difusion_ecs(Diseno("Nance 114 PEG", 114.0, -2.5))
    chequeo("C3 difusión: 114 nm con ζ −2.5 mV debe PASAR", r114.estado == PASA,
            f"da {r114.estado} (Nance Tabla 2, calificada +++)")

    # -- Nance 2012, Tabla 1: la de 198 nm con PEG (ζ −7.8 mV) NO se dispersa.
    r198 = g_difusion_ecs(Diseno("Nance 200 PEG", 198.0, -7.8))
    chequeo("C4 difusión: 198 nm con ζ −7.8 mV debe FALLAR", r198.estado == FALLA,
            f"da {r198.estado} (Nance: 1600x más lenta que en ACSF)")

    # -- Weinbaum 2003: la albúmina (7 nm) atraviesa el tamiz.
    r_alb = g_glicocalix_tamiz(Diseno("albúmina", 7.0, 0.0))
    chequeo("C5 la albúmina (7 nm) pasa el tamiz del glicocálix",
            r_alb.estado == PASA)

    # -- el fármaco liberado no se juzga por tamaño y no hay dato para juzgarlo
    #    por portador: la compuerta debe salir DESCONOCIDA.
    #    REESCRITO el 2026-08-17 con la compuerta (tarea F.1): antes exigía PASA
    #    (eje equivocado, tamaño); una primera versión de la reescritura exigió
    #    FALLA, que era un salto de "la gelsolina no lo une" a "no tiene
    #    portador". Ver el comentario de g_difusion_ecs.
    r_far = g_difusion_ecs(Diseno("test", 700.0, 0.0), "fármaco liberado")
    chequeo("C6 la difusión del FTY720-P liberado sale DESCONOCIDA (falta medida "
            "de portador)",
            r_far.estado == DESCONOCIDA,
            f"da {r_far.estado} (Bucki 2010 descarta la gelsolina, pero albúmina "
            "y apoM no están medidas para esta molécula)")

    titulo("BLOQUE 2 · LA REESCRITURA DE LA COMPUERTA (tarea V-5)",
           "Casos donde el criterio nuevo y el viejo dan resultados DISTINTOS.")

    # -- El caso que justifica la reescritura. Nance Tabla 1: la partícula
    #    nominal de "40 nm" con PEG mide 69 nm y ζ −2.8 mV, y es de las que más
    #    rápido difunden (solo 37x más lenta que en ACSF). El criterio viejo de
    #    Thorne (38 nm) la habría dado por EXCLUIDA. Es un dato medido.
    d69 = Diseno("Nance 40 PEG (mide 69 nm)", 69.0, -2.8)
    nuevo = g_difusion_ecs(d69)
    viejo = g_difusion_ecs(d69, escenario="thorne")
    chequeo("R1 la de 69 nm y ζ −2.8 mV PASA con el criterio nuevo",
            nuevo.estado == PASA, f"da {nuevo.estado}")
    chequeo("R2 ...y FALLABA con el criterio viejo de solo tamaño",
            viejo.estado == FALLA,
            f"da {viejo.estado}  <- esto es lo que corrige la reescritura")

    # -- El otro lado: pequeña pero adhesiva NO difunde. Nance Tabla 1, la COOH
    #    nominal de "40 nm" mide 57 nm con ζ −38.6 mV y quedó inmovilizada.
    #    El criterio viejo la juzgaba solo por tamaño.
    d57 = Diseno("Nance 40 COOH (mide 57 nm)", 57.0, -38.6)
    chequeo("R3 pequeña pero adhesiva (57 nm, ζ −38.6 mV) debe FALLAR",
            g_difusion_ecs(d57).estado == FALLA,
            "(Nance: el 100% de las COOH inmovilizadas, incluidas las de 40 nm)")

    # -- Conservadurismo declarado, NO es un fallo. La de 106 nm con ζ −4.4 mV
    #    cae en la banda −4 a −6 mV, donde Nance reporta 1 de 2. Difundía, pero
    #    la compuerta devuelve DESCONOCIDA a propósito: no se aprueba una
    #    superficie cuya estadística publicada es una moneda al aire.
    r106 = g_difusion_ecs(Diseno("Nance 100 PEG (mide 106 nm)", 106.0, -4.4))
    chequeo("R4 la de 106 nm y ζ −4.4 mV cae en zona gris de superficie",
            r106.estado == DESCONOCIDA,
            "(conservador a propósito: Nance da 1 de 2 en esa banda)")

    titulo("BLOQUE 3 · FALSABILIDAD",
           "Que la lógica no pueda repartir aprobados gratis.")

    # -- una ruta con compuerta desconocida NO puede salir aprobada.
    #    OJO: antes esta prueba usaba la ruta B, pero al cerrarse la tarea B.3
    #    la ruta B dejó de tener incógnitas y la prueba perdió sentido. Se pasa
    #    a la ruta C, que sigue con dos compuertas sin dato. El diseño es
    #    sigiloso y de tamaño cómodo para que NINGUNA compuerta falle: lo único
    #    que puede impedir el aprobado es la existencia de DESCONOCIDAS.
    #    REESCRITA el 2026-08-17i, y OTRA VEZ el 2026-08-17l. Al retirarse la
    #    cadena que cerraba G.2, la carga útil vuelve a DESCONOCIDA y la ruta C
    #    vuelve a salir NO EVALUABLE, que es justo lo que esta prueba quiere
    #    comprobar. Ya no hace falta quitar la compuerta para aislar el
    #    mecanismo: se evalúa la ruta C entera.
    v, res_c = evaluar_ruta(Diseno("sigiloso", 50.0, -2.0),
                            RUTAS["C · adsortiva (carga positiva)"])
    chequeo("F1 una DESCONOCIDA impide el aprobado",
            v == "NO EVALUABLE" and not any(r.estado == FALLA for r in res_c),
            f"da '{v}'")
    # F1b REESCRITA el 2026-08-17l. Antes comprobaba que G.2 excluía la ruta C.
    # Ese resultado se retractó. Lo que SÍ excluye una ruta hoy es un argumento
    # de mecanismo, y solo en la ruta D: el mismo diseño que sale NO EVALUABLE
    # por la C tiene que salir EXCLUIDO por la D.
    v_d, res_d = evaluar_ruta(Diseno("sigiloso", 50.0, -2.0),
                              RUTAS["D · mediada por receptor"])
    chequeo("F1b el mismo diseño queda EXCLUIDO por la ruta D",
            v_d == "EXCLUIDA"
            and any(r.compuerta == "Transcitosis por TfR de alta afinidad"
                    and r.estado == FALLA for r in res_d),
            f"da '{v_d}' (por mecanismo del receptor, no por dosis)")

    # -- un diseño imposible queda excluido en todas las rutas.
    veredictos = {k: vv for k, (vv, _) in evaluar(Diseno("imposible", 3.0, 0.0)).items()}
    chequeo("F2 un liposoma de 3 nm queda EXCLUIDO en todas las rutas",
            all(vv == "EXCLUIDA" for vv in veredictos.values()))

    # -- ζ POSITIVO debe salir DESCONOCIDA, nunca PASA. Es el caso de los tres
    #    diseños del estudio, y es una extrapolación fuera del rango medido.
    r_pos = g_difusion_ecs(Diseno("catiónico", 35.0, +6.7))
    chequeo("F3 un diseño de ζ POSITIVO no puede salir aprobado en difusión",
            r_pos.estado == DESCONOCIDA, f"da {r_pos.estado}")

    # -- el suelo geométrico del liposoma NO puede aplicarse a otras clases.
    #    Un dendrímero de 5 nm es más pequeño que dos bicapas, pero eso no lo
    #    excluye: no tiene bicapa. Antes salía FALLA por un motivo ajeno.
    peq = 5.0
    chequeo("F5 un liposoma de 5 nm FALLA por su suelo geométrico",
            g_transportador_fabricable(Diseno("l", peq, 0.0)).estado == FALLA)
    chequeo("F6 un dendrímero de 5 nm NO se excluye por el suelo del liposoma",
            g_transportador_fabricable(
                Diseno("d", peq, 0.0, clase="dendrimero")).estado == PASA,
            "(G4–G5, dentro de la ventana G3–G10)")

    # -- G.1a: la ventana del dendrímero tiene DOS lados y el que manda es el
    #    techo. G10 es la última generación completable.
    def _dend(x):
        return g_transportador_fabricable(Diseno("d", x, 0.0, clase="dendrimero"))

    chequeo("F7 un dendrímero de 30 nm FALLA por su techo geométrico",
            _dend(30.0).estado == FALLA,
            f"(techo G10 = {DENDRIMERO_TECHO_nm} nm)")
    chequeo("F8 un dendrímero por debajo de G3 sale DESCONOCIDA, no FALLA",
            _dend(2.0).estado == DESCONOCIDA,
            "(G1–G2 sí complejan, falta su diámetro medido)")
    chequeo("F9 el techo del dendrímero NO alcanza la ventana de envolvimiento",
            DENDRIMERO_TECHO_nm < g_envolvimiento(Diseno("x", 20.0, 0.0)).umbral,
            f"(margen {g_envolvimiento(Diseno('x', 20.0, 0.0)).umbral - DENDRIMERO_TECHO_nm:.2f} nm — EMPATE TÉCNICO)")
    chequeo("F10 el dendrímero cabe por el glicocálix y el liposoma no",
            _dend(8.0).estado == PASA
            and g_glicocalix_tamiz(Diseno("d", 8.0, 0.0)).estado == PASA
            and g_transportador_fabricable(Diseno("l", 8.0, 0.0)).estado == FALLA)
    chequeo("F11 un dendrímero fabricable siempre avisa de la carga 1:1",
            "1:1" in _dend(8.0).advertencia)

    # -- G.1b: polímero macizo.
    def _poli(x, mw=53.0):
        return g_transportador_fabricable(
            Diseno("p", x, 0.0, clase="polimerico", masa_molar_kDa=mw))

    chequeo("F12 un polímero sin masa molar declarada sale DESCONOCIDA",
            g_transportador_fabricable(
                Diseno("p", 50.0, 0.0, clase="polimerico")).estado == DESCONOCIDA,
            "(sin M no hay suelo que calcular)")
    chequeo("F13 el suelo del polímero CRECE con la masa molar",
            _poli(50.0, 15.0).umbral < _poli(50.0, 60.0).umbral,
            f"({_poli(50.,15.).umbral:.2f} nm a 15 kDa · "
            f"{_poli(50.,60.).umbral:.2f} nm a 60 kDa)")
    chequeo("F14 un polímero más pequeño que una cadena colapsada FALLA",
            _poli(2.0).estado == FALLA,
            f"(suelo {_poli(2.0).umbral:.2f} nm para 53 kDa)")
    chequeo("F15 el suelo del polímero es INSENSIBLE al fármaco (<2 % con 10 moléculas)",
            abs(diametro_globulo_colapsado_nm(53000 + 10 * 307.48, 1.19)
                / diametro_globulo_colapsado_nm(53000, 1.19) - 1.0) < 0.02)
    _mic = lambda x: g_transportador_fabricable(
        Diseno("m", x, 0.0, clase="micela"))
    chequeo("F16 una micela por encima del suelo cargado PASA, con salvedad",
            _mic(20.0).estado == PASA and bool(_mic(20.0).advertencia),
            f"(suelo {MICELA_SUELO_CARGADA_nm} nm, Sochor 2020)")
    chequeo("F17 una micela por debajo del suelo sale DESCONOCIDA, no FALLA",
            _mic(5.0).estado == DESCONOCIDA,
            "(el suelo viene de otro polímero y otro fármaco: G.4 y G.5)")
    chequeo("F18 una clase que no existe NO se evalúa por defecto",
            g_transportador_fabricable(
                Diseno("x", 20.0, 0.0, clase="inventada")).estado == DESCONOCIDA)

    # -- B.3 reabierta el 2026-08-12. El tránsito de Tong (24-48 h) y la
    #    descarga de Mao (>24 h, sin techo) se SOLAPAN, así que la compuerta no
    #    se puede decidir. Antes devolvía PASA apoyada en los 20 h de Yona, que
    #    miden salida de circulación y no llegada a la lesión.
    chequeo("F19 el tránsito del monocito sale DESCONOCIDA: los intervalos se solapan",
            g_transito_vs_liberacion(Diseno("x", 700.0, 0.0)).estado == DESCONOCIDA,
            f"(tránsito {T_TRANSITO_PRIMERA_DETECCION_h:.0f}-{T_TRANSITO_PICO_h:.0f} h "
            f"frente a descarga >{T_LIBERACION_COTA_INFERIOR_h:.0f} h)")
    # TRIPWIRE, no invariante. Fija el estado del modelo mientras B.3 siga
    # abierta. El día que B.3 se cierre con dato, o que entre al catálogo un
    # diseño que sobreviva, ESTA PRUEBA DEBE FALLAR: es la señal de que hay que
    # actualizarla a mano, no un error.
    chequeo("F20 con B.3 abierta, NINGÚN diseño del catálogo sale NO EXCLUIDO",
            all(v != "NO EXCLUIDA"
                for d in CATALOGO for v, _ in evaluar(d).values()),
            "(tripwire: si falla, B.3 se cerró o hay un diseño nuevo vivo; "
            "revisar y actualizar la prueba)")

    # -- rama positiva del eje de superficie, cerrada el 2026-08-13 con Berry
    #    2016 y Mastorakos 2016. El umbral es el MÍNIMO medido, no un umbral
    #    físico: por debajo NO se extrapola.
    _zeta = lambda z: g_difusion_ecs(Diseno("x", 50.0, z))
    chequeo("F21 un ζ igual o mayor que el mínimo positivo medido FALLA",
            _zeta(ZETA_ADHESIVO_POSITIVO_mV).estado == FALLA
            and _zeta(35.3).estado == FALLA,
            f"(+{ZETA_ADHESIVO_POSITIVO_mV:.1f} mV Berry, +35.3 mV Mastorakos)")
    chequeo("F22 un ζ positivo por debajo del mínimo medido sale DESCONOCIDA",
            all(_zeta(z).estado == DESCONOCIDA for z in (0.5, 2.0, 6.7, 9.9)),
            "(entre 0 y +10 mV no hay ni un dato: no se extrapola)")
    chequeo("F23 la rama positiva NO toca la negativa ni el ζ neutro",
            _zeta(-0.24).estado == PASA and _zeta(-28.33).estado == FALLA,
            "(Chow 2025 y Gong 2022 conservan su veredicto)")

    # -- zona gris de TAMAÑO. Sigue abierta tras leer Curtis 2019 y McKenna
    #    2021: el único punto dentro de la banda es de rata P14 y con una
    #    métrica de diámetro distinta. TRIPWIRE: si algún día se cierra, esta
    #    prueba DEBE fallar para obligar a revisarla.
    chequeo("F24 la zona gris de tamaño 114-200 nm sigue sin decidir",
            g_difusion_ecs(Diseno("x", 134.0, -0.24)).estado == DESCONOCIDA,
            "(tripwire: Chow 2025; si falla, la banda se cerró con dato)")
    chequeo("F25 el factor de edad P14->P70 se deriva de McKenna, no a mano",
            abs(EDAD_FACTOR_P14_A_P70 - 34.0 / 5.0) < 1e-9,
            f"(x{EDAD_FACTOR_P14_A_P70:.1f}; poro {EDAD_PORO_EFECTIVO_nm['P14']} "
            f"-> {EDAD_PORO_EFECTIVO_nm['P70']} nm)")

    # -- carga útil, añadida el 2026-08-13. Sigue DESCONOCIDA, pero ya NO por
    #    tres números: quedan resueltos la carga del liposoma (Mouzoura 2025) y
    #    el umbral en parénquima (Foster 2007). Falta UN eslabón, el mismo de la
    #    transcitosis. Razón actualizada el 2026-08-17h.
    chequeo("F26 la carga útil está en las CUATRO rutas",
            all(any(c(CATALOGO[0]).compuerta == "Carga útil suficiente"
                    for c in comps) for comps in RUTAS.values()))
    # F27 CAMBIÓ DOS VECES. El 2026-08-17i pasó de DESCONOCIDA a FALLA al entrar
    # Johnsen 2019. El 2026-08-17l vuelve a DESCONOCIDA: la auditoría retiró la
    # cadena entera. Esta prueba es ahora un TRIPWIRE DE RETRACTACIÓN — si
    # alguien vuelve a cerrar G.2 sin los dos números que faltan, salta.
    chequeo("F27 la carga útil es DESCONOCIDA para TODO el catálogo",
            all(g_carga_util(d).estado == DESCONOCIDA for d in CATALOGO),
            "(retirada el 17l: falta el umbral eficaz de fingolimod y la "
            "fracción entregada de un liposoma cargado con fingolimod)")
    chequeo("F27b la retractación queda escrita en el motivo",
            "RETRACTACIÓN" in g_carga_util(CATALOGO[0]).motivo,
            "(quien lea la salida tiene que ver que hubo un resultado retirado, "
            "no un hueco que nunca se intentó llenar)")
    # Una DESCONOCIDA no puede leerse como aprobado. La advertencia lo dice, y
    # esta prueba impide que alguien la borre por parecer redundante.
    chequeo("F27c la DESCONOCIDA declara que NO es un aprobado",
            "NO es un aprobado" in g_carga_util(CATALOGO[0]).advertencia)
    # Tripwire anti-sesgo, reescrito: la fracción entregada ya NO decide nada,
    # así que las tres variantes tienen que dar exactamente lo mismo. Si alguien
    # reconecta una fracción a la decisión, esta prueba salta.
    chequeo("F27d la fracción entregada ya no decide el veredicto",
            (g_carga_util_cota(CATALOGO[0]).estado
             == g_carga_util_pasiva(CATALOGO[0]).estado == DESCONOCIDA)
            and g_carga_util_cota(CATALOGO[0]).valor is None,
            "(las tres variantes dan DESCONOCIDA sin valor numérico)")
    # Lo que SÍ decide hoy, y por un camino que no es de dosis.
    chequeo("F27e la ruta D queda EXCLUIDA por mecanismo del receptor",
            g_transcitosis_tfr(CATALOGO[0]).estado == FALLA,
            "(la vía del hierro no incluye transcitosis del TfR; RI7 y OX26 "
            "quedan confinados al endotelio capilar)")
    # La condicional NO es decorativa: es la palanca de diseño que deja la
    # compuerta abierta a un ligando de afinidad reducida. Si se pierde, el
    # veredicto se vuelve más fuerte de lo que la fuente permite.
    chequeo("F27f el FALLA del TfR conserva su condicional de afinidad",
            "afinidad" in g_transcitosis_tfr(CATALOGO[0]).advertencia.lower())
    # Y el tripwire que impide repetir el error de G.2: el veredicto del TfR NO
    # puede apoyarse en el '< 1/30', que es una comparación contra otro estudio.
    chequeo("F27g el veredicto del TfR no se apoya en el '< 1/30'",
            "1/30" not in g_transcitosis_tfr(CATALOGO[0]).motivo,
            "(esa cifra compara contra 1.6 %ID/g de la ref. [23], otro estudio, "
            "y el propio Johnsen la declara disputada)")

    # -- las advertencias no se pueden perder por el camino.
    chequeo("F4 un PASA con salvedad conserva su advertencia",
            bool(g_salida_farmaco(Diseno("x", 100.0, 0.0)).advertencia)
            and bool(r114.advertencia))

    if verbose:
        print("-" * 78)
        print(f" RESULTADO: {sum(ok)}/{len(ok)} pruebas superadas")
        print(" Recordatorio: el bloque 1 es consistencia, no validación.")
        print("=" * 78)
    return all(ok)


# =============================================================================
#  INFORME
# =============================================================================

# Liposomas TEÓRICOS: formulaciones propuestas, ningún número medido.
CATALOGO_TEORICO = [
    Diseno("Diseño furtivo (teórico)", 31.0, +2.0, 5.0, sintetico=True),
    Diseno("Diseño convencional (teórico)", 40.0, +5.0, 0.0, sintetico=True),
    Diseno("Diseño catiónico (teórico)", 35.0, +6.7, 0.0, sintetico=True),
]

# Liposomas REALES: Ø y ζ medidos y publicados, con la referencia en la nota.
CATALOGO_REAL = [
    Diseno("Mao 2014 (real)", 157.5, +3.99, 5.0, nota="Nanomedicine 10:393"),
    Diseno("Gong 2022 (real)", 145.0, -28.33, 5.0, nota="Nanophotonics 11:5133"),
    Diseno("Chow 2025 (real)", 134.0, -0.24, 0.0, nota="Drug Deliv Transl Res 15:2022"),
    Diseno("Muselman 2026 (real)", 700.0, 0.0, 5.0, nota="Front Immunol 16:1657131"),
]

# El catálogo completo se mantiene: es el que usan el informe de texto, las
# pruebas y la equivalencia con el Colab. El orden es el de siempre.
CATALOGO = CATALOGO_TEORICO + CATALOGO_REAL

_SIMB = {PASA: "✓", FALLA: "✗", DESCONOCIDA: "?"}
_CORTO = {"EXCLUIDA": "NO", "NO EXCLUIDA": "SÍ", "NO EVALUABLE": "??"}


def _quien_lo_mata(res):
    """Devuelve el nombre de la primera compuerta que falla, o None."""
    for r in res:
        if r.estado == FALLA:
            return r.compuerta
    return None


def _lectura_corta(veredictos, resultados, nombres):
    """La frase más informativa para un diseño, en una línea.

    Antes esto imprimía "la primera compuerta que lo mata en alguna ruta", y
    salía SIEMPRE la de la ruta A, o sea el tamiz del glicocálix, para los siete
    diseños. Además decía "se cae por el tamiz" de Muselman, cuyo interés es
    justo el contrario: que la ruta B NO lo excluye. Ahora se prioriza el mejor
    desenlace, que es el que de verdad informa.
    """
    vivas = [n for n in nombres if veredictos[n] == "NO EXCLUIDA"]
    if vivas:
        return "CANDIDATO por " + ", ".join(v.split(" ·")[0] for v in vivas)

    grises = [n for n in nombres if veredictos[n] == "NO EVALUABLE"]
    if grises:
        return "faltan datos en " + ", ".join(g.split(" ·")[0] for g in grises)

    # excluido en todas: el motivo útil es qué compuerta lo mata más veces
    from collections import Counter
    matadores = Counter(m for n in nombres
                        if (m := _quien_lo_mata(resultados[n])))
    peor, veces = matadores.most_common(1)[0]
    peor = _ETIQUETA.get(peor, peor)
    return f"excluido en todas · {peor} ({veces}/{len(nombres)})"


def resumen(catalogo=None, pistas=True):
    """SALIDA SIMPLE. Una línea por diseño. Es la vista por defecto."""
    catalogo = catalogo or CATALOGO
    ancho = max(len(d.nombre) for d in catalogo)
    print()
    print("  ¿PUEDE ESTE DISEÑO LLEGAR A LA MIELINA?")
    print("  " + "-" * (ancho + 63))
    print(f"  {'':{ancho}s}  {'A':^4s} {'B':^4s} {'C':^4s} {'D':^4s}   lectura")
    print("  " + "-" * (ancho + 63))
    nombres = list(RUTAS.keys())
    for d in catalogo:
        r = evaluar(d)
        veredictos = {n: r[n][0] for n in nombres}
        resultados = {n: r[n][1] for n in nombres}
        marcas = " ".join(f"{_CORTO[veredictos[n]]:^4s}" for n in nombres)
        print(f"  {d.nombre:{ancho}s}  {marcas}   "
              f"{_lectura_corta(veredictos, resultados, nombres)}")
    print("  " + "-" * (ancho + 63))
    print("  A pasiva · B macrófago · C adsortiva · D receptor")
    print("  NO = excluido (firme)   SÍ = candidato (débil)   ?? = faltan datos")
    if not pistas:
        print()
        return
    print()
    print("  Para el desglose completo:  python3 rutas.py --detalle")
    print("  Para las figuras:           python3 rutas.py --figuras")
    print()


# Nombres cortos SOLO para imprimir. Los de las compuertas son descriptivos a
# propósito, pero algunos pasan de 50 caracteres y descuadraban la columna.
_ETIQUETA = {
    "Tránsito del monocito frente a cinética de liberación":
        "Tránsito del monocito vs liberación",
    "Salida del fármaco de la célula transportadora":
        "Salida del fármaco de la célula",
    "Difusión en espacio extracelular (transportador)":
        "Difusión extracelular · transportador",
    "Difusión en espacio extracelular (fármaco liberado)":
        "Difusión extracelular · fármaco",
}


def informe(catalogo=None, detalle=False, tabla=True):
    """Solo resultados. La prosa vive en notas().

    tabla=False omite la tabla de veredictos: la usa correr.sh, que ya imprime
    la versión compacta de resumen() y no tiene por qué repetir lo mismo en
    dos formatos distintos.
    """
    catalogo = catalogo or CATALOGO
    nombres = list(RUTAS.keys())

    # Las salvedades se numeran y se listan UNA vez al final. Antes se repetían
    # enteras dentro de cada diseño: el mismo párrafo siete veces.
    notas_txt, orden = {}, []

    def _nota(t):
        if t not in notas_txt:
            orden.append(t)
            notas_txt[t] = len(orden)
        return notas_txt[t]

    if tabla:
        anc = max(len(d.nombre) for d in catalogo)
        print()
        print(" VEREDICTO POR RUTA")
        print(" " + "-" * (anc + 56))
        print(f" {'diseño':{anc}s} " + " ".join(f"{n.split(' ·')[0]:>13s}" for n in nombres))
        print(" " + "-" * (anc + 56))
        for d in catalogo:
            r = evaluar(d)
            print(f" {d.nombre:{anc}s} " + " ".join(f"{r[n][0]:>13s}" for n in nombres))
        print(" " + "-" * (anc + 56))
        print(" EXCLUIDA = no puede    NO EXCLUIDA = candidato    NO EVALUABLE = faltan datos")

    if detalle:
        print()
        print(" DESGLOSE COMPUERTA A COMPUERTA")

    if detalle:
        for d in catalogo:
            print()
            print(f" {d.nombre}   Ø {d.diametro_nm:g} nm · ζ {d.zeta_mV:+.2f} mV · PEG {d.peg_nm:g} nm")
            print(" " + "-" * 74)
            for nombre, (v, res) in evaluar(d).items():
                print(f"   {nombre.split(' (')[0]:<38s} {v}")
                for x in res:
                    txt = _ETIQUETA.get(x.compuerta, x.compuerta)
                    if x.valor is not None and x.umbral is not None:
                        dato = f"{x.valor:.1f} / {x.umbral:.1f} {x.unidad}".strip()
                    else:
                        dato = ""
                    ref = f" [{_nota(x.advertencia)}]" if x.advertencia else ""
                    print(f"     {_SIMB[x.estado]} {txt:<40s} "
                          f"{dato:>16s}{ref}".rstrip())

    if orden:
        print()
        print(" SALVEDADES")
        print(" " + "-" * 74)
        import textwrap as _tw
        for t in orden:
            n = notas_txt[t]
            lineas = _tw.wrap(t, 70)
            print(f"  [{n}] {lineas[0]}")
            for ln in lineas[1:]:
                print(f"      {ln}")
    print()


def notas():
    """Toda la prosa: lectura de los resultados y qué cambió en el modelo.

    Vive aparte porque el informe debe imprimir resultados, no explicaciones.
    Se llega con:  sh correr.sh notas
    """
    print("=" * 79)
    print(" LECTURA")
    print("=" * 79)
    print(" 1. NINGUNA ruta sale NO EXCLUIDA. La B lo estuvo entre el 2026-08-10 y")
    print("    el 2026-08-12, para un diseño de 700 nm, pero la tarea B.3 volvió a")
    print("    abrirse: el tránsito del monocito al cerebro inflamado (24-48 h,")
    print("    Tong 2016) se solapa con la descarga del liposoma (>24 h sin techo")
    print("    medido, Mao 2014), así que no hay decisión. Sus otras dos últimas")
    print("    compuertas siguen llevando salvedad.")
    print(" 2. Las rutas A, C y D incluyen la transcitosis completa, que es transporte")
    print("    activo dependiente de ATP y queda fuera del alcance del método. Mientras")
    print("    siga ahí, ninguna de las tres puede salir de 'NO EVALUABLE'.")
    print(" 3. La ruta B es la única cuyas incógnitas son cerrables con literatura.")
    print("    B.5, salida del fármaco de la célula, se cerró con Hisano 2011. B.3,")
    print("    tránsito del monocito, sigue ABIERTA: hacen falta la cinética de")
    print("    reclutamiento a una lesión de EAE a resolución de horas y la semivida")
    print("    de descarga de un liposoma PEGilado de ~700 nm. Ninguna existe hoy.")
    print(" 4. Los tres diseños teóricos quedan EXCLUIDOS de la ruta A por el tamiz")
    print("    del glicocálix, y de la B por captación fagocítica insuficiente.")
    print(" 5. La ruta B es además la única donde lo que tiene que difundir por el")
    print("    espacio extracelular es el FÁRMACO (1.68 nm) y no el transportador.")
    print("    Por eso Muselman 2026, de 700 nm, sale EXCLUIDO de C y D pero no de B:")
    print("    una partícula de 700 nm no atraviesa el parénquima, su carga sí.")
    print(" 6. LA CARGA ÚTIL es incógnita en las CUATRO rutas, y se añadió el")
    print("    2026-08-13. Hasta entonces el modelo preguntaba si el transportador")
    print("    LLEGA y nunca si entrega FÁRMACO SUFICIENTE: un diseño que llegara con")
    print("    una sola molécula habría salido igual de bien que uno con miles. Era un")
    print("    hueco INVISIBLE. Al 2026-08-17 faltaba UN solo eslabón, no tres:")
    print("    la FRACCIÓN de la dosis que llega al parénquima, que es el MISMO")
    print("    hueco de la transcitosis y no uno independiente. Ya están resueltos")
    print("    la carga del liposoma (1:8 molar, eficiencia 94-97.2 %, Mouzoura 2025,")
    print("    = 4.3-4.4 % p/p de fingolimod) y el umbral en parénquima (398 ± 186")
    print("    ng/g de FTY720-P a dosis eficaz en EAE, Foster 2007 Tabla 3). Para el")
    print("    dendrímero falta además su carga con fingolimod (el 1:1 de Devarakonda")
    print("    2004 es con nifedipino y es un complejo de solubilización), pero el")
    print("    dendrímero quedó fuera del foco de simulación el 2026-08-17.")
    print()
    print("=" * 79)
    print(" LA COMPUERTA DE DIFUSIÓN CAMBIÓ (tarea V-5, 2026-08-10)")
    print("=" * 79)
    print(" Antes era solo de TAMAÑO, con el 38 nm de Thorne & Nicholson 2006.")
    print(" Ahora es de DOS VARIABLES, tamaño y superficie, siguiendo el diagrama de")
    print(" fases de Nance et al. 2012, que impugna ese 38 nm por subestimación:")
    print(" compraron el mismo lote de puntos cuánticos que usó Thorne, no difundían,")
    print(" los recubrieron más con PEG y entonces sí.")
    print()
    print("   eje 1: superficie (ζ)              eje 2: tamaño, solo si es sigilosa")
    _sup = [("ζ > 0", "DESCONOCIDA"),
            (f"0 >= ζ >= {ZETA_DIFUSIVO_mV:.0f}", "sigilosa -> eje 2"),
            (f"{ZETA_DIFUSIVO_mV:.0f} > ζ > {ZETA_NO_DIFUSIVO_mV:.0f}", "DESCONOCIDA"),
            (f"ζ <= {ZETA_NO_DIFUSIVO_mV:.0f}", "FALLA")]
    _tam = [(f"Ø <= {D_SIGILOSO_PASA_nm:.0f} nm", "PASA"),
            (f"{D_SIGILOSO_PASA_nm:.0f} < Ø < {D_SIGILOSO_FALLA_nm:.0f} nm", "DESCONOCIDA"),
            (f"Ø >= {D_SIGILOSO_FALLA_nm:.0f} nm", "FALLA"),
            ("", "")]
    for (sa, sb), (ta, tb) in zip(_sup, _tam):
        izq = f"{sa:>14s} -> {sb:<20s}" if sa else " " * 38
        der = f"{ta:>14s} -> {tb}" if ta else ""
        print(f"   {izq}{der}")
    print()
    print(" CONSECUENCIA INCÓMODA PARA ESTE PROYECTO. Los tres diseños tienen ζ")
    print(" POSITIVO (+2.0, +5.0, +6.7 mV) y el barrido de Nance va de −2.5 a −52 mV.")
    print(" Sí hay dato con carga positiva, pero MÁS ARRIBA: Berry 2016 mide +10.0 mV")
    print(" (menos del 10 % de la población difunde) y Mastorakos 2016 mide +35.3 mV")
    print(" (inmovilizada), las dos por MPT en cerebro de rata ex vivo. Con ζ ≥ +10 mV")
    print(" la compuerta FALLA. Entre 0 y +10 mV no hay ni un dato y NO se extrapola:")
    print(" devuelve DESCONOCIDA. Los tres diseños caen justo en ese hueco.")
    print(" SALVEDAD VIVA: Berry y Mastorakos son polímero/ADN, no liposomas, y su ζ")
    print(" está medido en NaCl 10 mM, no en aCSF; en aCSF ambos pierden estabilidad")
    print(" coloidal, así que adhesión y agregación no están separadas.")
    print()
    print(" Y hay una TENSIÓN DE DISEÑO que conviene tener delante: la carga positiva")
    print(" es lo que propone la ruta C para ENTRAR (adsorción al glicocálix aniónico)")
    print(" y es lo que estorbaría para DIFUNDIR una vez dentro. Dos tramos del")
    print(" recorrido con requisitos de signo opuesto. Refuerza la idea de que el")
    print(" transportador no tenga que hacer los dos.")
    print()
    print(" En cuánto difieren los dos criterios, para partículas sigilosas (ζ = −2 mV):")
    for d_test in (31.0, 40.0, 69.0, 106.0, 114.0, 150.0, 200.0):
        dd = Diseno("t", d_test, -2.0)
        n = g_difusion_ecs(dd).estado
        t = g_difusion_ecs(dd, escenario="thorne").estado
        marca = "  <- difieren" if n != t else ""
        print(f"    Ø {d_test:6.1f} nm ->  Nance: {n:11s} | Thorne (38 nm): {t:11s}{marca}")
    print("=" * 79)


# =============================================================================
#  FIGURAS
# =============================================================================

def _ventana_fabricable(catalogo):
    """La barra de 'Transportador fabricable' de la figura de ventanas.

    ANTES estaba escrita a mano con el suelo del LIPOSOMA (12 nm), así que la
    figura no valía para dendrímero ni para polímero y esas clases se quedaban
    sin ella. Ahora sale de la propia compuerta, para la clase que tenga el
    catálogo que se esté dibujando.

    Devuelve (etiqueta, lo, hi, indefinida). `indefinida` es True cuando algún
    diseño del catálogo no tiene ventana verificada (p. ej. un dendrímero que no
    es PAMAM, o un polímero sin masa molar declarada): entonces la barra se
    pinta rayada y no debe leerse como un rango cerrado.
    """
    catalogo = catalogo or CATALOGO
    clases = {d.clase for d in catalogo}
    clase = clases.pop() if len(clases) == 1 else None
    indefinida = any(g_transportador_fabricable(d).estado == DESCONOCIDA
                     for d in catalogo)
    TOPE = 2000.0

    if clase == "liposoma":
        return ("Transportador fabricable (liposoma)",
                G.diametro_liposoma_minimo_nm(4.0, 4.0), TOPE, indefinida)

    if clase == "dendrimero":
        # La ventana solo está verificada para PAMAM: si el catálogo trae otra
        # química, la etiqueta no debe decir PAMAM.
        etq = ("Transportador fabricable (dendrímero)" if indefinida
               else "Transportador fabricable (dendrímero PAMAM)")
        return (etq, DENDRIMERO_SUELO_nm, DENDRIMERO_TECHO_nm, indefinida)

    if clase == "micela":
        return ("Transportador fabricable (micela cargada)",
                MICELA_SUELO_CARGADA_nm, TOPE, indefinida)

    if clase == "polimerico":
        suelos = [diametro_globulo_colapsado_nm(
                      d.masa_molar_kDa * 1000.0,
                      d.densidad_g_cm3 or POLIMERO_DENSIDAD_POR_DEFECTO)
                  for d in catalogo if d.masa_molar_kDa]
        if suelos:
            return ("Transportador fabricable (polímero macizo)",
                    min(suelos), TOPE, indefinida or len(suelos) != len(catalogo))
        return ("Transportador fabricable (polímero macizo)", 1.0, TOPE, True)

    return ("Transportador fabricable", 1.0, TOPE, True)


def figuras(prefijo="rutas", catalogo=None, incluir_ventanas=True):
    """Genera las tres figuras del simulador: ventanas, matriz y recorrido.

    Las tres son GENÉRICAS: valen para cualquier catálogo y para cualquier
    clase. `incluir_ventanas=False` se conserva por compatibilidad, pero ya no
    hace falta usarlo para las clases que no son liposoma.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    catalogo = catalogo or CATALOGO
    VERDE, ROJO, GRIS = "#2e7d32", "#c62828", "#9e9e9e"

    if incluir_ventanas:
        # ---------------------------------------------------------------- FIGURA 1
        #  Ventanas de tamaño. Cada barra es el rango PERMITIDO por una compuerta.
        #  Es la figura que muestra de un vistazo por qué la ventana está vacía.
        etq_fab, lo_fab, hi_fab, fab_indef = _ventana_fabricable(catalogo)
        _u = _umbrales_modelo()
        ventanas = [
            (etq_fab, lo_fab, hi_fab, GRIS if fab_indef else VERDE),
            ("Tamiz del glicocálix",   1.0, _u["glicocalix"], VERDE),
            ("Envolvimiento de membrana", _u["envolvimiento"], 2000.0, VERDE),
            ("Compuerta de caveola",              1.0,   80.0, VERDE),
            ("Difusión extracelular (portador),\n"
             "  solo si ζ entre −4 y 0 mV",       1.0, D_SIGILOSO_PASA_nm, VERDE),
            ("Captación por macrófago",         550.0, 2000.0, VERDE),
            # La barra "Difusión extracelular (fármaco)" se RETIRÓ el 2026-08-17
            # (tarea F.1). Dibujaba una ventana de 1 a 38 nm (Thorne) que daba a
            # entender que el fármaco liberado difunde si es pequeño. Ya no es
            # así: la compuerta pasó a juzgar dependencia de portador, no
            # tamaño, y FALLA en todo el rango de diámetros. Dejar la barra
            # verde habría contradicho al propio modelo. Ahora esta compuerta
            # aparece en el pie, entre las que no dependen del tamaño.
        ]
        # Las etiquetas van como marcas del eje Y, NO como texto encima de la barra.
        # Puestas encima quedaban a media altura entre dos barras y no se sabía a
        # cuál correspondía cada una.
        fig, ax = plt.subplots(figsize=(12.5, 5.6))
        etiquetas_y = []
        for i, (nombre, lo, hi, col) in enumerate(ventanas):
            y = len(ventanas) - 1 - i
            rayada = (i == 0 and fab_indef)
            ax.barh(y, hi - lo, left=lo, height=0.62, color=col, alpha=0.65,
                    edgecolor="black", linewidth=0.5,
                    hatch=("//" if rayada else None))
            etiquetas_y.append((y, nombre + ("\n  (sin ventana verificada)"
                                             if rayada else "")))
        # Las dos zonas grises van sobre la barra que les toca. El índice se
        # BUSCA por nombre, no se escribe a mano: al retirar una barra el
        # 2026-08-17 los índices literales (1 y 2) habrían quedado apuntando a
        # la barra equivocada sin que ninguna prueba lo detectara.
        _y_de = {nom: len(ventanas) - 1 - i
                 for i, (nom, _lo, _hi, _c) in enumerate(ventanas)}
        # zona gris de la captación fagocítica
        ax.barh(_y_de["Captación por macrófago"], 550 - 150, left=150,
                height=0.55, color=GRIS, alpha=0.45,
                hatch="//", edgecolor="black", linewidth=0.5)
        # zona gris de la difusión extracelular del portador: 114 a 200 nm
        ax.barh(_y_de["Difusión extracelular (portador),\n"
                      "  solo si ζ entre −4 y 0 mV"],
                D_SIGILOSO_FALLA_nm - D_SIGILOSO_PASA_nm, left=D_SIGILOSO_PASA_nm,
                height=0.55, color=GRIS, alpha=0.45, hatch="//", edgecolor="black",
                linewidth=0.5)

        # Las líneas de los diseños. AVISO IMPORTANTE: esta figura es de UNA sola
        # variable (el diámetro), pero la compuerta de difusión pasó a ser de DOS
        # (tamaño y superficie). Para un diseño de ζ positivo esa barra NO aplica:
        # la compuerta devuelve DESCONOCIDA sin llegar a mirar el tamaño. Si se
        # pintaran todas las líneas igual, esta figura diría que los tres diseños
        # teóricos "caben" en la ventana de difusión, que es justo lo contrario de
        # lo que concluye el modelo. Por eso se distinguen.
        n_indef = por_zeta = por_tamano = 0
        for d in catalogo:
            rdif = g_difusion_ecs(d)
            indef = rdif.estado == DESCONOCIDA
            n_indef += indef
            # No todas las DESCONOCIDA de esta compuerta son por ζ. La de ζ
            # ni llega a mirar el tamaño; la de la zona gris 114-200 nm es justo
            # lo contrario, ζ vale y lo que no decide es el tamaño. Decirlas
            # iguales en el pie era falso.
            # La causa se lee de la UNIDAD del propio Resultado ("mV" -> salió
            # por la rama de superficie, "nm" -> por la de tamaño), no se
            # reconstruye a mano la condición: así no puede desincronizarse de
            # g_difusion_ecs cuando esa función cambie. Ya pasó al añadir la
            # rama positiva el 2026-08-13.
            if indef:
                if rdif.unidad == "mV":
                    por_zeta += 1
                else:
                    por_tamano += 1
            ax.axvline(d.diametro_nm, color=(GRIS if indef else ROJO),
                       ls=(":" if indef else "--"), lw=1.0, alpha=0.75)
            # dentro de los ejes: fuera se montaba sobre el título
            ax.text(d.diametro_nm, len(ventanas) - 0.62,
                    f"{d.diametro_nm:.0f}", rotation=90, fontsize=7,
                    color=(GRIS if indef else ROJO), va="bottom", ha="right")

        # AVISO OBLIGATORIO. Esta figura solo dibuja las compuertas que TIENEN
        # ventana de tamaño. Las que no dependen del diámetro no aparecen, y sin
        # decirlo la figura se lee como un veredicto: un diseño cuya línea cae en
        # verde en todas las barras de una ruta parece que la supera, cuando
        # puede estar bloqueado por una compuerta que aquí no se ve. Los nombres
        # se derivan de RUTAS, no se escriben a mano.
        _con_ventana = {
            "Transportador fabricable", "Tamiz del glicocálix",
            "Envolvimiento de membrana", "Compuerta de caveola",
            "Captación fagocítica",
            "Difusión en espacio extracelular (transportador)"}
        # "Difusión en espacio extracelular (fármaco liberado)" salió de este
        # conjunto el 2026-08-17 (tarea F.1): ya no tiene ventana de tamaño, así
        # que pasa a listarse en el pie como compuerta que no depende del
        # diámetro.
        _sin_ventana, _vistas = [], set()
        for _comps in RUTAS.values():
            for _c in _comps:
                _n = _c(catalogo[0]).compuerta
                if _n not in _con_ventana and _n not in _vistas:
                    _vistas.add(_n)
                    _sin_ventana.append(_ETIQUETA.get(_n, _n))

        pie = []
        if por_zeta:
            # OJO con la redacción: la BARRA verde dice cuándo la compuerta
            # PASA (ζ entre −4 y 0 mV) y este pie dice cuándo sale DESCONOCIDA.
            # Son cosas distintas y si se enuncian como un solo rango la figura
            # se lee contradictoria. Hay DOS bandas sin dato, una a cada lado.
            # El signo menos va con el tipográfico "−", como en el resto de la
            # figura; un f-string de un número negativo saca el ASCII "-" y las
            # dos formas juntas en el mismo pie quedan desiguales.
            pie.append(f"Punteado gris por ζ ({por_zeta}): el ζ cae en una "
                       f"banda SIN DATO (de −{abs(ZETA_NO_DIFUSIVO_mV):.0f} a "
                       f"−{abs(ZETA_DIFUSIVO_mV):.0f} mV, o de 0 a "
                       f"+{ZETA_ADHESIVO_POSITIVO_mV:.0f} mV); la compuerta de "
                       "difusión ni mira el tamaño.")
        if por_tamano:
            pie.append(f"Punteado gris por TAMAÑO ({por_tamano}): ζ válido, pero el "
                       f"diámetro cae en la franja rayada de "
                       f"{D_SIGILOSO_PASA_nm:.0f}–{D_SIGILOSO_FALLA_nm:.0f} nm.")
        if _sin_ventana:
            pie.append("NO son todas las compuertas: estas no dependen del tamaño y "
                       "no se dibujan → " + " · ".join(_sin_ventana) + ".")
            pie.append("Que una línea caiga en verde en todas las barras de una ruta "
                       "NO significa que la supere. El veredicto está en la matriz.")
        if pie:
            ax.text(0.5, -0.20, "\n".join(pie), transform=ax.transAxes,
                    ha="center", va="top", fontsize=7.5, color="#555555")
        # La marca azul del "fingolimod libre (1.68 nm)" se RETIRÓ el 2026-08-17
        # (tarea F.1) junto con su barra. Señalaba dónde caía el fármaco dentro
        # de una ventana de tamaño que ya no existe; sin la barra, la marca
        # sugería que el tamaño del fármaco sigue decidiendo algo. No decide.

        ax.set_xscale("log")
        ax.set_xlim(1, 2000)
        ax.set_ylim(-0.8, len(ventanas) - 0.2)
        ax.set_yticks([y for y, _ in etiquetas_y])
        ax.set_yticklabels([n for _, n in etiquetas_y], fontsize=9)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("diámetro (nm)")
        # "de cada compuerta" era falso: solo salen las que dependen del tamaño.
        ax.set_title("Ventanas de tamaño · solo las compuertas que dependen del "
                     "diámetro\n"
                     "barra = rango permitido · línea vertical = diseño evaluado")
        ax.grid(axis="x", alpha=0.3, which="both")
        # La entrada azul "permitido (fármaco liberado)" salió el 2026-08-17
        # (tarea F.1): ya no hay ninguna barra azul que explicar.
        ax.legend(handles=[
            Patch(facecolor=VERDE, alpha=0.65, label="permitido (portador)"),
            Patch(facecolor=GRIS, alpha=0.45, hatch="//", label="zona sin dato"),
        ], loc="lower left", fontsize=8)
        fig.tight_layout(rect=[0, 0.10, 1, 1])   # hueco para la nota de abajo
        fig.savefig(f"{prefijo}_ventanas.png", dpi=160, bbox_inches="tight")
        plt.close(fig)

    # ---------------------------------------------------------------- FIGURA 2
    #  Matriz de veredictos: diseños por rutas.
    nombres = list(RUTAS.keys())
    col = {"EXCLUIDA": ROJO, "NO EXCLUIDA": VERDE, "NO EVALUABLE": GRIS}
    fig, ax = plt.subplots(figsize=(9, 0.62 * len(catalogo) + 2.4))
    for i, d in enumerate(catalogo):
        r = evaluar(d)
        for j, n in enumerate(nombres):
            v = r[n][0]
            ax.add_patch(plt.Rectangle((j, len(catalogo) - 1 - i), 0.94, 0.9,
                                       facecolor=col[v], alpha=0.75,
                                       edgecolor="white"))
            # Las salvedades solo se marcan cuando la ruta NO está excluida:
            # ahí es donde restan fuerza a la afirmación. En una casilla "NO"
            # son irrelevantes, porque la ruta queda descartada igualmente, y
            # ponerlas solo añadía ruido.
            n_salv = (sum(1 for x in r[n][1] if x.advertencia)
                      if v == "NO EXCLUIDA" else 0)
            etiqueta = _CORTO[v] + "!" * min(n_salv, 3)
            ax.text(j + 0.47, len(catalogo) - 1 - i + 0.45, etiqueta,
                    ha="center", va="center", fontsize=11, color="white",
                    fontweight="bold")
    ax.set_xlim(0, len(nombres)); ax.set_ylim(0, len(catalogo))
    ax.set_xticks([j + 0.47 for j in range(len(nombres))])
    ax.set_xticklabels([n.split(" (")[0] for n in nombres], fontsize=8)
    ax.set_yticks([len(catalogo) - 1 - i + 0.45 for i in range(len(catalogo))])
    ax.set_yticklabels([d.nombre for d in catalogo], fontsize=8)
    ax.set_title("¿Puede este diseño llegar a la mielina por esta ruta?")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.legend(handles=[
        Patch(facecolor=ROJO, alpha=0.75, label="NO · excluido (firme)"),
        Patch(facecolor=VERDE, alpha=0.75, label="SÍ · candidato (débil)"),
        Patch(facecolor=GRIS, alpha=0.75, label="?? · faltan datos"),
    ], loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=8,
        frameon=False)
    ax.text(0.5, -0.155, "Cada '!' de un SÍ = una compuerta que pasa pero con "
            "salvedad declarada. SÍ!!! se apoya en tres.",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
            color="#555555")
    fig.tight_layout()
    fig.savefig(f"{prefijo}_matriz.png", dpi=160)
    plt.close(fig)

    # ---------------------------------------------------------------- FIGURA 3
    #  Recorrido: en qué compuerta muere cada diseño, ruta por ruta.
    # Alta y ancha: con 7 diseños y hasta 7 compuertas por ruta, con el tamaño
    # anterior (11 x 2.9 por ruta) las etiquetas se solapaban y no se leía nada.
    fig, axes = plt.subplots(len(RUTAS), 1, figsize=(14, 3.6 * len(RUTAS)),
                             sharex=False)
    for ax, (nombre, comps) in zip(axes, RUTAS.items()):
        etiquetas = [c(catalogo[0]).compuerta for c in comps]
        for i, d in enumerate(catalogo):
            _, res = evaluar_ruta(d, comps)
            for j, r in enumerate(res):
                c = {PASA: VERDE, FALLA: ROJO, DESCONOCIDA: GRIS}[r.estado]
                # Un PASA CON SALVEDAD no puede pintarse igual que un PASA
                # limpio: en el texto se distingue con '✓!'. Se marca como
                # ROSQUILLA (mismo color, centro hueco). Antes llevaba un anillo
                # NEGRO, pero en la ruta B casi todas las compuertas tienen
                # salvedad y el anillo aparecía en todo: dejaba de distinguir
                # nada y solo ensuciaba la figura.
                # OJO: solo un PASA se pinta hueco. Una DESCONOCIDA también
                # puede llevar advertencia (B.3 desde el 2026-08-12), y
                # pintarla hueca la haría pasar por un aprobado con salvedad,
                # que es justo lo contrario de lo que dice.
                if r.advertencia and r.estado == PASA:
                    ax.scatter(j, len(catalogo) - 1 - i, s=105, color=c,
                               edgecolor="white", linewidth=1.0, zorder=3)
                    ax.scatter(j, len(catalogo) - 1 - i, s=30, color="white",
                               zorder=4)
                else:
                    ax.scatter(j, len(catalogo) - 1 - i, s=95, color=c,
                               edgecolor="white", zorder=3)
            ax.plot(range(len(res)), [len(catalogo) - 1 - i] * len(res),
                    color="#cccccc", lw=1, zorder=1)
        ax.set_title(nombre, fontsize=11, loc="left", fontweight="bold")
        import textwrap
        ax.set_xticks(range(len(etiquetas)))
        ax.set_xticklabels(["\n".join(textwrap.wrap(e, 22)) for e in etiquetas],
                           fontsize=8.5)
        ax.set_yticks(range(len(catalogo)))
        ax.set_yticklabels([d.nombre for d in reversed(catalogo)], fontsize=8.5)
        ax.set_ylim(-0.7, len(catalogo) - 0.3)
        # margen simétrico: con -0.6 a la derecha la etiqueta de la última
        # compuerta se salía del eje y salía cortada al guardar.
        ax.set_xlim(-0.45, len(etiquetas) - 1 + 0.45)
        ax.grid(axis="x", alpha=0.25)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    fig.suptitle("Recorrido compuerta a compuerta\n"
                 "verde = pasa · rojo = falla · gris = sin dato\n"
                 "punto hueco = pasa PERO con salvedad declarada", fontsize=12)
    fig.tight_layout(rect=[0, 0.01, 1, 0.96])
    fig.savefig(f"{prefijo}_recorrido.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"\n  Figuras guardadas:")
    if incluir_ventanas:
        print(f"    {prefijo}_ventanas.png    ventanas de tamaño de cada compuerta")
    print(f"    {prefijo}_matriz.png      qué diseño puede usar qué ruta")
    print(f"    {prefijo}_recorrido.png   dónde se cae cada diseño, paso a paso")


# =============================================================================
#  DENDRÍMERO — informe y figura propios  (tarea G.1a)
#
#  Van aparte a propósito. El CATÁLOGO y las tres figuras de arriba son de
#  liposomas y ya están citadas en el informe del equipo; meter el dendrímero
#  dentro cambiaría tablas que ya están hechas. Aquí se ve solo esta clase.
# =============================================================================

def _umbrales_modelo():
    """Los tres umbrales de tamaño, leídos del propio código, no a mano."""
    return dict(
        glicocalix=2.0 * G.radio_exclusion_nm(),
        envolvimiento=g_envolvimiento(Diseno("_", 20.0, 0.0)).umbral,
        fagocitosis=550.0,
    )


def informe_dendrimero():
    """Resultado de G.1a en texto, sin prosa."""
    u = _umbrales_modelo()
    techo, suelo = DENDRIMERO_TECHO_nm, DENDRIMERO_SUELO_nm
    banda = (DENDRIMERO_TECHO_MEDIDO_nm * (1 - DENDRIMERO_PRECISION),
             DENDRIMERO_TECHO_MEDIDO_nm * (1 + DENDRIMERO_PRECISION))

    print("=" * 79)
    print(" DENDRÍMERO PAMAM · ventana geométrica (tarea G.1a)")
    print("=" * 79)
    print(f"  suelo  G3   {suelo:6.2f} nm   Prosa 2001 (medido) + Devarakonda 2004 (aloja fármaco)")
    print(f"  techo  G10  {techo:6.2f} nm   Maiti 2004 (calc.) · medido por Prosa: "
          f"{DENDRIMERO_TECHO_MEDIDO_nm:.2f} nm")
    print()
    print("  GENERACIÓN   Ø (nm)   fabricable   ¿pasa el glicocálix?   ¿se envuelve?")
    print("  " + "-" * 74)
    for g, dn in DENDRIMERO_GENERACIONES_nm.items():
        r = g_transportador_fabricable(Diseno(f"G{g}", dn, 0.0, clase="dendrimero"))
        gl = "sí" if dn <= u["glicocalix"] else "no"
        en = "sí" if dn >= u["envolvimiento"] else "no"
        print(f"  G{g:<10} {dn:6.2f}   {r.estado:11s}  {gl:^20s}   {en:^12s}")
    print()
    print("  UMBRALES DEL MODELO (leídos del código)")
    print(f"    tamiz del glicocálix        d <= {u['glicocalix']:7.3f} nm   (Weinbaum)")
    print(f"    envolvimiento de membrana   d >= {u['envolvimiento']:7.3f} nm   (Deserno)")
    print(f"    captación fagocítica        d >= {u['fagocitosis']:7.1f} nm   (ruta B)")
    print()
    print(f"  RESULTADO  el techo del dendrímero queda {u['envolvimiento'] - techo:.2f} nm "
          "por debajo del envolvimiento.")
    print("             Ningún dendrímero puede ser envuelto, ni captado por macrófago.")
    print()
    print("  !! EMPATE TÉCNICO. Prosa declara ±5 % de precisión global. Sobre sus")
    print(f"     {DENDRIMERO_TECHO_MEDIDO_nm:.2f} nm de G10 eso da {banda[0]:.2f}–{banda[1]:.2f} nm, "
          f"y {banda[1]:.2f} > {u['envolvimiento']:.3f}.")
    print("     Firme en valores centrales, marginal en el extremo. No es un caso cerrado.")
    print()
    print("  Salvedades y fuentes: verificacion/verificacion_dendrimero_tarea_G_1a.md")
    print("=" * 79)


# Constantes de estabilidad K(1:1) del complejo dendrímero–nifedipino.
# Devarakonda et al. 2004, Int J Pharm 284:133, Tabla 2. Unidades M⁻¹.
# Generación entera = superficie de aminas (–NH2); media = superficie de éster
# (–COOCH3). G0 no aparece en su tabla: no formó complejo medible.
DEVARAKONDA_K = {
    "amina":  {4: {}, 7: {1: (25.6, 1.1), 2: (52.7, 2.1), 3: (287.6, 6.1)},
               10: {1: (18.0, 1.3), 2: (27.7, 1.6), 3: (187.1, 2.3)}},
    "ester":  {4: {0.5: (23.4, 1.5), 1.5: (79.6, 2.2), 2.5: (116.3, 3.3)},
               7: {0.5: (42.5, 1.3), 1.5: (91.2, 2.5), 2.5: (338.7, 5.6)},
               10: {0.5: (44.8, 1.1), 1.5: (95.7, 1.2), 2.5: (338.5, 4.5)}},
}


def catalogo_dendrimero():
    """Las ocho generaciones con diámetro MEDIDO, como diseños evaluables.

    OJO CON ζ: un PAMAM de aminas terminales es catiónico a pH fisiológico, pero
    NO tenemos un valor medido en fuente primaria para esta clase. Se pone un
    ζ positivo simbólico y no un número concreto: +1.0 mV cae en el hueco sin
    dato (0 a +10 mV) y la compuerta devuelve DESCONOCIDA.
    OJO, esto SÍ depende de la magnitud desde el 2026-08-13: si el ζ medido
    resultara ≥ +10 mV (Berry 2016), la compuerta pasaría a FALLA y el veredicto
    del dendrímero cambiaría. Un PAMAM-NH2 de generación alta puede perfectamente
    estar por encima de +10 mV, así que este valor simbólico NO es inocuo:
    es la tarea G.3 y hay que medirlo antes de dar por bueno el resultado.
    """
    return [Diseno(f"PAMAM G{g}", dn, +1.0, clase="dendrimero",
                   nota="ζ>0 sin valor medido")
            for g, dn in DENDRIMERO_GENERACIONES_nm.items()]


# -----------------------------------------------------------------------------
#  DENDRÍMEROS TEÓRICOS DEL SIMULADOR DE ACOPLE  (entrada de Jhovan, 2026-08-12)
#
#  Tres fichas generadas con el simulador de acople fármaco–nanotransportador
#  (dinámica Browniana sobre PMF radial derivado del logP). Fármaco: fingolimod,
#  logP 4.16.
#
#  LOS VALORES SON SINTÉTICOS. Palabras de Jhovan: "la mayoría de los datos son
#  inventados dentro de un rango teórico real". Entran al catálogo marcados con
#  sintetico=True para ver DÓNDE CHOCAN con las compuertas, y NO pueden citarse
#  como resultado ni cerrar ninguna tarea de verificación.
#
#  El ζ tampoco cierra el hueco del ζ del dendrímero: sigue sin valor medido.
# -----------------------------------------------------------------------------

DENDRIMEROS_TEORICOS = [
    ("PAMAM G4 (ficha)",   15.0,  +9.9, "pamam",
     "Ø de ficha; Prosa 2001 mide 4.60 nm para G4"),
    ("PPI G4 (ficha)",      4.6,  +7.0, "ppi",       "Ø de ficha"),
    ("Carbosilano G2",      2.63, +0.3, "carbosilano", "Ø de ficha"),
]


def catalogo_dendrimeros_teoricos():
    """Los tres dendrímeros de las fichas de acople, tal cual, sin corregir."""
    return [Diseno(nom, dn, z, clase="dendrimero", subquimica=sq,
                   sintetico=True, nota=f"DATO SINTÉTICO · {nota}")
            for nom, dn, z, sq, nota in DENDRIMEROS_TEORICOS]


# -----------------------------------------------------------------------------
#  POLÍMEROS Y MICELAS TEÓRICOS DEL SIMULADOR DE ACOPLE  (Jhovan, 2026-08-12)
#
#  Tres fichas más del simulador de acople, con fingolimod (logP 4.16). VALORES
#  SINTÉTICOS, igual que los dendrímeros: no son medidas.
#
#  Dos de las tres son MICELAS, no polímeros macizos, y van con clase="micela".
#
#  MASA MOLAR DEL PLGA: la ficha no la declara. Se DERIVA de sus propios datos
#  (504 cadenas, 232 869 monómeros → 462 monómeros por cadena) suponiendo la
#  razón 85:15, que es la única con densidad medida por Parker 2010:
#      462 × (0.85·72.06 + 0.15·58.04) = 462 × 69.96 = 32.3 kDa
#  La razón casi no importa: con 50:50 el suelo baja de 4.42 a 4.31 nm (2.5 %) y
#  ningún veredicto cambia. Queda marcado como DERIVADO, no como dato de ficha.
# -----------------------------------------------------------------------------

PLGA_MONOMEROS_POR_CADENA = 232869 / 504      # de la propia ficha
PLGA_M_MONOMERO_85_15 = 0.85 * 72.06 + 0.15 * 58.04    # g/mol
PLGA_MW_DERIVADA_kDa = PLGA_MONOMEROS_POR_CADENA * PLGA_M_MONOMERO_85_15 / 1000.0

POLIMEROS_TEORICOS = [
    ("Nanopartícula PLGA (ficha)", 33.3, +4.7, "polimerico",
     PLGA_MW_DERIVADA_kDa, 1.19,
     f"Mw {PLGA_MW_DERIVADA_kDa:.1f} kDa DERIVADA de la ficha (462 monómeros/"
     f"cadena, razón 85:15 supuesta), no declarada"),
    ("Micela PEG-PLA (ficha)", 11.0, +8.4, "micela", None, None,
     "micela, no polímero macizo · núcleo PLA 3.5 nm + corona PEG 2.0 nm"),
    ("Micela PCL-Pluronic-PCL", 13.4, 0.0, "micela", None, None,
     "micela, no polímero macizo · núcleo PCL 4.2 nm + corona PEO 2.5 nm"),
]


def catalogo_polimeros_teoricos():
    """Las tres fichas de polímero/micela, tal cual, sin corregir."""
    return [Diseno(nom, dn, z, clase=cl, masa_molar_kDa=mw,
                   densidad_g_cm3=rho, sintetico=True,
                   nota=f"DATO SINTÉTICO · {nota}")
            for nom, dn, z, cl, mw, rho, nota in POLIMEROS_TEORICOS]


def informe_polimeros_teoricos():
    """Los tres polímeros/micelas teóricos contra las compuertas geométricas."""
    disenos = catalogo_polimeros_teoricos()
    puertas = [g_transportador_fabricable, g_glicocalix_tamiz, g_envolvimiento]

    print()
    print("=" * 78)
    print("  POLÍMEROS Y MICELAS TEÓRICOS DEL ACOPLE — DATOS SINTÉTICOS")
    print("=" * 78)
    print("  Valores inventados dentro de un rango teórico. NO son medidas y no")
    print("  cierran ninguna tarea de verificación.")
    print()
    print(f"  {'transportador':28s} {'Ø nm':>7s} {'ζ mV':>6s} {'clase':>11s}  "
          f"{'fabric.':>11s} {'glicoc.':>11s} {'envolv.':>11s}")
    print("  " + "-" * 92)
    for d in disenos:
        est = [p(d).estado for p in puertas]
        print(f"  {d.nombre:28s} {d.diametro_nm:7.2f} {d.zeta_mV:+6.1f} "
              f"{d.clase:>11s}  {est[0]:>11s} {est[1]:>11s} {est[2]:>11s}")

    print()
    print("  POR QUÉ:")
    for d in disenos:
        print(f"  · {d.nombre}")
        for p in puertas:
            r = p(d)
            det = ""
            if r.valor is not None and r.umbral is not None:
                det = f"  ({r.valor:g} vs {r.umbral:.2f} {r.unidad})"
            print(f"      {r.compuerta:26s} {r.estado:12s}{det}")
            if r.motivo:
                print(f"          {r.motivo}")
            if r.advertencia:
                print(f"          salvedad: {r.advertencia}")
    print()
    print(f"  Mw del PLGA: {PLGA_MW_DERIVADA_kDa:.1f} kDa DERIVADA "
          f"({PLGA_MONOMEROS_POR_CADENA:.0f} monómeros/cadena × "
          f"{PLGA_M_MONOMERO_85_15:.2f} g/mol, razón 85:15 supuesta).")
    print("  ζ: ninguno de los tres entra en el modelo. Siguen sin valor medido.")
    print()


def figura_polimeros_teoricos(nombre="polimeros_teoricos.png"):
    """Los tres sobre las dos ventanas, con el suelo de su propia clase."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    u = _umbrales_modelo()
    VERDE, ROJO, GRIS, AZUL = "#2e7d32", "#c62828", "#9e9e9e", "#1565c0"
    COLOR = {PASA: VERDE, FALLA: ROJO, DESCONOCIDA: GRIS}
    disenos = catalogo_polimeros_teoricos()

    fig, ax = plt.subplots(figsize=(12.5, 4.8))
    ax.axvspan(0, u["glicocalix"], color=AZUL, alpha=0.09)
    ax.axvspan(u["envolvimiento"], 40, color=VERDE, alpha=0.09)
    ax.axvline(u["glicocalix"], color=AZUL, ls="--", lw=1.4)
    ax.axvline(u["envolvimiento"], color=VERDE, ls="--", lw=1.4)
    ax.axvline(MICELA_SUELO_CARGADA_nm, color=GRIS, ls=":", lw=1.3)

    for i, d in enumerate(disenos):
        y = len(disenos) - 1 - i
        est = g_transportador_fabricable(d).estado
        # Mismo convenio que la figura de recorrido: un PASA CON SALVEDAD se
        # pinta HUECO. Aquí se pintaba macizo y las dos figuras del mismo
        # diseño se contradecían.
        _r = g_transportador_fabricable(d)
        _hueco = bool(_r.advertencia) and est == PASA
        ax.plot([d.diametro_nm], [y], "o", ms=13, color=COLOR[est],
                mec=(COLOR[est] if _hueco else "black"),
                mew=(2.2 if _hueco else 0.8), zorder=5,
                mfc=("white" if _hueco else COLOR[est]))
        ax.text(d.diametro_nm, y + 0.22, f"{d.diametro_nm:g} nm",
                ha="center", fontsize=9)

    y_top = len(disenos) - 0.45
    for x, txt, col in ((u["glicocalix"], "tamiz glicocálix", AZUL),
                        (MICELA_SUELO_CARGADA_nm, "suelo micela cargada", GRIS),
                        (u["envolvimiento"], "envolvimiento", VERDE)):
        ax.text(x, y_top, f" {txt}", rotation=90, va="top", ha="left",
                fontsize=8, color=col)

    ax.set_yticks(range(len(disenos)))
    ax.set_yticklabels([d.nombre for d in reversed(disenos)], fontsize=9.5)
    ax.set_xlabel("diámetro (nm)")
    ax.set_xlim(0, 40)
    ax.set_ylim(-0.7, len(disenos) - 0.3)
    ax.set_title("Polímeros y micelas teóricos (DATOS SINTÉTICOS) contra las dos "
                 "ventanas\n"
                 f"glicocálix ≤ {u['glicocalix']:.2f} nm · envolvimiento ≥ "
                 f"{u['envolvimiento']:.2f} nm · suelo de la micela cargada "
                 f"{MICELA_SUELO_CARGADA_nm} nm", fontsize=11)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mec="black", ms=9,
                              color=COLOR[e], label=f"fabricable: {e}")
                       for e in (PASA, FALLA, DESCONOCIDA)]
                      + [Line2D([], [], marker="o", ls="", ms=9, mew=2.2,
                                mec=VERDE, color=VERDE, mfc="white",
                                label="fabricable: PASA con salvedad")],
              loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(nombre, dpi=150)
    plt.close(fig)
    print(f"  figura escrita: {nombre}")


def figuras_polimeros_teoricos():
    """El set común más la figura propia de esta sub-página."""
    import contextlib as _c
    import io as _io

    figura_polimeros_teoricos()
    with _c.redirect_stdout(_io.StringIO()):
        figuras(prefijo="pol_teoricos", catalogo=catalogo_polimeros_teoricos())
    for f in ("ventanas", "matriz", "recorrido"):
        print(f"  figura escrita: pol_teoricos_{f}.png")


def informe_dendrimeros_teoricos():
    """Los tres dendrímeros teóricos contra las compuertas geométricas."""
    disenos = catalogo_dendrimeros_teoricos()
    puertas = [g_transportador_fabricable, g_glicocalix_tamiz, g_envolvimiento]

    print()
    print("=" * 78)
    print("  DENDRÍMEROS TEÓRICOS DEL SIMULADOR DE ACOPLE — DATOS SINTÉTICOS")
    print("=" * 78)
    print("  Valores inventados dentro de un rango teórico. NO son medidas y no")
    print("  cierran ninguna tarea de verificación.")
    print()
    print(f"  {'transportador':22s} {'Ø nm':>7s}  {'ζ mV':>6s}  "
          f"{'fabricable':>12s} {'glicocálix':>12s} {'envolvim.':>11s}")
    print("  " + "-" * 74)
    for d in disenos:
        est = [p(d).estado for p in puertas]
        print(f"  {d.nombre:22s} {d.diametro_nm:7.2f}  {d.zeta_mV:+6.1f}  "
              f"{est[0]:>12s} {est[1]:>12s} {est[2]:>11s}")

    print()
    print("  POR QUÉ:")
    for d in disenos:
        print(f"  · {d.nombre}")
        for p in puertas:
            r = p(d)
            det = ""
            if r.valor is not None and r.umbral is not None:
                det = f"  ({r.valor:g} vs {r.umbral:.2f} {r.unidad})"
            print(f"      {r.compuerta:26s} {r.estado:12s}{det}")
            if r.motivo:
                print(f"          {r.motivo}")
    print()
    print("  ζ: ninguno de los tres entra en el modelo. Siguen sin valor medido")
    print("     en fuente primaria (tarea abierta).")
    print("  Carga útil: el simulador de acople muestrea UNA molécula con N")
    print("     réplicas independientes, así que no da estequiometría [tarea G.2].")
    print()


def figura_dendrimeros_teoricos(nombre="dendrimeros_teoricos.png"):
    """Los tres dendrímeros teóricos sobre las dos ventanas del modelo."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    u = _umbrales_modelo()
    disenos = catalogo_dendrimeros_teoricos()
    VERDE, ROJO, GRIS = "#2e7d32", "#c62828", "#9e9e9e"
    COLOR = {PASA: VERDE, FALLA: ROJO, DESCONOCIDA: GRIS}

    fig, ax = plt.subplots(figsize=(12.5, 4.6))

    ax.axvspan(DENDRIMERO_SUELO_nm, u["glicocalix"], color=VERDE, alpha=0.10)
    ax.axvline(DENDRIMERO_SUELO_nm, color=GRIS, ls=":", lw=1.2)
    ax.axvline(DENDRIMERO_TECHO_nm, color=GRIS, ls=":", lw=1.2)
    ax.axvline(u["glicocalix"], color="#1565c0", ls="--", lw=1.4)

    for i, d in enumerate(disenos):
        y = len(disenos) - 1 - i
        est = g_transportador_fabricable(d).estado
        # Mismo convenio que la figura de recorrido: un PASA CON SALVEDAD se
        # pinta HUECO. Aquí se pintaba macizo y las dos figuras del mismo
        # diseño se contradecían.
        _r = g_transportador_fabricable(d)
        _hueco = bool(_r.advertencia) and est == PASA
        ax.plot([d.diametro_nm], [y], "o", ms=13, color=COLOR[est],
                mec=(COLOR[est] if _hueco else "black"),
                mew=(2.2 if _hueco else 0.8), zorder=5,
                mfc=("white" if _hueco else COLOR[est]))
        ax.text(d.diametro_nm, y + 0.24, f"{d.diametro_nm:g} nm",
                ha="center", fontsize=9)

    y_top = len(disenos) - 0.45
    for x, txt, col in ((DENDRIMERO_SUELO_nm, "suelo PAMAM G3", GRIS),
                        (u["glicocalix"], "tamiz glicocálix", "#1565c0"),
                        (DENDRIMERO_TECHO_nm, "techo PAMAM G10", GRIS)):
        ax.text(x, y_top, f" {txt}", rotation=90, va="top", ha="left",
                fontsize=8, color=col)

    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mec="black", ms=9,
                              color=COLOR[e], label=f"fabricable: {e}")
                       for e in (PASA, FALLA, DESCONOCIDA)]
                      + [Line2D([], [], marker="o", ls="", ms=9, mew=2.2,
                                mec=VERDE, color=VERDE, mfc="white",
                                label="fabricable: PASA con salvedad")],
              loc="lower right", fontsize=8, framealpha=0.9)

    ax.set_yticks(range(len(disenos)))
    ax.set_yticklabels([d.nombre for d in reversed(disenos)], fontsize=10)
    ax.set_xlabel("diámetro (nm)")
    ax.set_xlim(0, 18)
    ax.set_ylim(-0.7, len(disenos) - 0.3)
    ax.set_title("Dendrímeros teóricos (DATOS SINTÉTICOS) contra la ventana "
                 "geométrica del modelo\n"
                 f"suelo PAMAM {DENDRIMERO_SUELO_nm} nm · tamiz del glicocálix "
                 f"{u['glicocalix']:.2f} nm · techo PAMAM "
                 f"{DENDRIMERO_TECHO_nm} nm",
                 fontsize=11)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(nombre, dpi=150)
    plt.close(fig)
    print(f"  figura escrita: {nombre}")


def figura_teoricos_vs_medidos(nombre="dend_teoricos_vs_medidos.png"):
    """Las tres fichas sobre la escala de generaciones PAMAM MEDIDA.

    Sirve para ver de un golpe que el Ø de la ficha del PAMAM G4 (15 nm) no cae
    donde caería un G4 medido (4.60 nm) sino por encima del G10.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    u = _umbrales_modelo()
    VERDE, ROJO, GRIS, AZUL = "#2e7d32", "#c62828", "#9e9e9e", "#1565c0"
    COLOR = {PASA: VERDE, FALLA: ROJO, DESCONOCIDA: GRIS}
    disenos = catalogo_dendrimeros_teoricos()

    fig, ax = plt.subplots(figsize=(12.5, 5.0))

    ax.axvspan(DENDRIMERO_SUELO_nm, u["glicocalix"], color=VERDE, alpha=0.09)
    ax.axvline(u["glicocalix"], color=AZUL, ls="--", lw=1.4)
    ax.axvline(DENDRIMERO_TECHO_nm, color=GRIS, ls=":", lw=1.3)

    # escala de referencia: las 8 generaciones con Ø medido (Prosa 2001)
    for g, dn in DENDRIMERO_GENERACIONES_nm.items():
        ax.plot([dn], [0], "|", ms=18, color="#37474f", mew=1.6, zorder=4)
        ax.text(dn, -0.30, f"G{g}", ha="center", va="top", fontsize=8,
                color="#37474f")

    for i, d in enumerate(disenos):
        y = i + 1
        est = g_transportador_fabricable(d).estado
        # Mismo convenio que la figura de recorrido: un PASA CON SALVEDAD se
        # pinta HUECO. Aquí se pintaba macizo y las dos figuras del mismo
        # diseño se contradecían.
        _r = g_transportador_fabricable(d)
        _hueco = bool(_r.advertencia) and est == PASA
        ax.plot([d.diametro_nm], [y], "o", ms=13, color=COLOR[est],
                mec=(COLOR[est] if _hueco else "black"),
                mew=(2.2 if _hueco else 0.8), zorder=5,
                mfc=("white" if _hueco else COLOR[est]))
        ax.text(d.diametro_nm, y + 0.20, f"{d.diametro_nm:g} nm",
                ha="center", fontsize=9)

    # el PAMAM de la ficha frente al PAMAM G4 medido
    d_ficha = disenos[0].diametro_nm
    d_medido = DENDRIMERO_GENERACIONES_nm[4]
    ax.annotate("", xy=(d_medido, 1), xytext=(d_ficha, 1),
                arrowprops=dict(arrowstyle="->", color=ROJO, lw=1.4,
                                linestyle="--"))
    ax.text((d_ficha + d_medido) / 2, 1.28,
            f"G4 medido = {d_medido:.2f} nm (Prosa 2001)", ha="center",
            fontsize=8.5, color=ROJO)

    ax.set_yticks([0] + [i + 1 for i in range(len(disenos))])
    ax.set_yticklabels(["PAMAM medido\n(Prosa 2001)"]
                       + [d.nombre for d in disenos], fontsize=9)
    ax.set_xlabel("diámetro (nm)")
    ax.set_xlim(0, 18)
    ax.set_ylim(-0.85, len(disenos) + 0.75)
    ax.set_title("Fichas teóricas frente a la escala de generaciones PAMAM medida\n"
                 f"banda verde = suelo {DENDRIMERO_SUELO_nm} nm a tamiz del "
                 f"glicocálix {u['glicocalix']:.2f} nm · línea punteada = techo "
                 f"G10 {DENDRIMERO_TECHO_nm} nm", fontsize=11)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mec="black", ms=9,
                              color=COLOR[e], label=f"fabricable: {e}")
                       for e in (PASA, FALLA, DESCONOCIDA)]
                      + [Line2D([], [], marker="o", ls="", ms=9, mew=2.2,
                                mec=VERDE, color=VERDE, mfc="white",
                                label="fabricable: PASA con salvedad")],
              loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(nombre, dpi=150)
    plt.close(fig)
    print(f"  figura escrita: {nombre}")


def figura_teoricos_pmf(nombre="dend_teoricos_pmf.png"):
    """El modelo de acople de las fichas, con sus propias fórmulas.

    Panel A: el pozo hidrófobo SOLO depende del logP del fármaco, así que vale
    lo mismo para los tres transportadores. Panel B: la barrera electrostática
    sí depende del ζ, pero es dos órdenes menor que el pozo.

    Fórmulas tomadas literalmente de las fichas de Jhovan:
      ΔG_partición = −2.303 · k_BT · logP
      B_elec       = z · 2ζ / (k_BT/e),  con k_BT/e = 26.714 mV a 310 K
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    LOGP_FARMACO, LOGP_CONTROL, Z_FARMACO = 4.16, -1.5, 1.0
    KT_mV = 1.380649e-23 * 310.0 / 1.602176634e-19 * 1000.0   # 26.714 mV
    KT_kcal = 1.380649e-23 * 310.0 / 4184.0 * 6.02214076e23   # 0.616 kcal/mol
    VERDE, ROJO, GRIS, AZUL = "#2e7d32", "#c62828", "#9e9e9e", "#1565c0"

    disenos = catalogo_dendrimeros_teoricos()
    dg_farmaco = -2.303 * LOGP_FARMACO

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    # ---- Panel A: el pozo no distingue transportador -----------------------
    lp = np.linspace(-2.0, 6.0, 200)
    ax1.plot(lp, -2.303 * lp, color=AZUL, lw=2.0)
    ax1.axhline(0.0, color=GRIS, lw=1.0)
    for x, col, txt in ((LOGP_FARMACO, VERDE, "fingolimod"),
                        (LOGP_CONTROL, ROJO, "control hidrofílico")):
        y = -2.303 * x
        ax1.plot([x], [y], "o", ms=10, color=col, mec="black", mew=0.8, zorder=5)
        ax1.annotate(f"{txt}\nlogP {x:g} → {y:+.2f} k_BT",
                     (x, y), textcoords="offset points", xytext=(8, -4),
                     fontsize=8.5, color=col, va="top")
    ax1.set_xlabel("logP del fármaco")
    ax1.set_ylabel("ΔG partición agua→interior (k_BT)")
    ax1.set_title("El pozo depende SOLO del logP del fármaco\n"
                  f"los tres transportadores comparten ΔG = {dg_farmaco:.2f} "
                  f"k_BT = {dg_farmaco*KT_kcal:.2f} kcal/mol", fontsize=10.5)
    ax1.grid(alpha=0.25)

    # ---- Panel B: la barrera electrostática, a escala del pozo -------------
    nombres = [d.nombre for d in disenos]
    barreras = [Z_FARMACO * 2.0 * d.zeta_mV / KT_mV for d in disenos]
    y = np.arange(len(disenos))
    ax2.barh(y, barreras, height=0.5, color=ROJO, alpha=0.75,
             edgecolor="black", linewidth=0.5)
    for i, (b, d) in enumerate(zip(barreras, disenos)):
        ax2.text(b + 0.15, i, f"{b:+.3f} k_BT   (ζ {d.zeta_mV:+.1f} mV)",
                 va="center", fontsize=8.5)
    ax2.axvline(0.0, color=GRIS, lw=1.0)
    ax2.axvline(-dg_farmaco, color=VERDE, lw=1.8)
    ax2.text(-dg_farmaco - 0.15, -0.42,
             f"profundidad del pozo: {-dg_farmaco:.2f} k_BT",
             color=VERDE, fontsize=8.5, va="bottom", ha="right")
    ax2.set_ylim(-0.75, len(disenos) - 0.25)
    ax2.set_yticks(y)
    ax2.set_yticklabels(nombres, fontsize=9)
    ax2.set_xlim(0, -dg_farmaco * 1.12)
    ax2.set_xlabel("barrera electrostática (k_BT)")
    ax2.set_title("La barrera de superficie sí depende del ζ,\n"
                  "pero es ≪ pozo: ralentiza, no revierte", fontsize=10.5)
    ax2.grid(axis="x", alpha=0.25)

    fig.suptitle("Modelo de acople de las fichas — DATOS SINTÉTICOS",
                 fontsize=11.5, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(nombre, dpi=150)
    plt.close(fig)
    print(f"  figura escrita: {nombre}")


def figuras_dendrimeros_teoricos():
    """Las cinco figuras de la sub-página de dendrímeros teóricos."""
    import contextlib as _c
    import io as _io

    figura_dendrimeros_teoricos()
    figura_teoricos_vs_medidos()
    figura_teoricos_pmf()
    with _c.redirect_stdout(_io.StringIO()):
        figuras(prefijo="dend_teoricos",
                catalogo=catalogo_dendrimeros_teoricos())
    for f in ("ventanas", "matriz", "recorrido"):
        print(f"  figura escrita: dend_teoricos_{f}.png")


def figura_dendrimero(nombre="dendrimero_ventana.png"):
    """Dos paneles: la ventana completa y el zoom del techo (el empate técnico)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    u = _umbrales_modelo()
    techo, suelo = DENDRIMERO_TECHO_nm, DENDRIMERO_SUELO_nm
    lo_b = DENDRIMERO_TECHO_MEDIDO_nm * (1 - DENDRIMERO_PRECISION)
    hi_b = DENDRIMERO_TECHO_MEDIDO_nm * (1 + DENDRIMERO_PRECISION)
    VERDE, ROJO, AZUL, GRIS = "#2e7d32", "#c62828", "#1565c0", "#9e9e9e"

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.5, 8.2))

    # ---- PANEL 1: la ventana entera, escala log ----------------------------
    ax1.barh(2, techo - suelo, left=suelo, height=0.5, color=VERDE, alpha=0.7,
             edgecolor="black", linewidth=0.6)
    ax1.text(suelo * 1.02, 2.42, f"dendrímero PAMAM: {suelo:.2f} – {techo:.2f} nm  (G3 – G10)",
             fontsize=9, color=VERDE, va="bottom")
    ax1.barh(1, u["glicocalix"] - 1.0, left=1.0, height=0.5, color=AZUL,
             alpha=0.55, edgecolor="black", linewidth=0.6)
    ax1.barh(0, 2000 - u["envolvimiento"], left=u["envolvimiento"], height=0.5,
             color=AZUL, alpha=0.55, edgecolor="black", linewidth=0.6)
    ax1.barh(-1, 2000 - u["fagocitosis"], left=u["fagocitosis"], height=0.5,
             color=AZUL, alpha=0.55, edgecolor="black", linewidth=0.6)

    for g, dn in DENDRIMERO_GENERACIONES_nm.items():
        cabe = dn <= u["glicocalix"]
        ax1.plot([dn], [2], marker="|", ms=13, mew=1.6,
                 color=("black" if cabe else ROJO))
        ax1.text(dn, 1.72, f"G{g}", fontsize=7, ha="center", va="top",
                 color=("black" if cabe else ROJO))

    ax1.set_xscale("log")
    ax1.set_xlim(1, 2000)
    ax1.set_ylim(-1.7, 3.0)
    ax1.set_yticks([2, 1, 0, -1])
    ax1.set_yticklabels(["existe el dendrímero", "tamiz del glicocálix",
                         "envolvimiento", "captación fagocítica"], fontsize=9)
    ax1.tick_params(axis="y", length=0)
    ax1.set_xlabel("diámetro (nm)")
    ax1.set_title("Dendrímero PAMAM frente a las ventanas del modelo\n"
                  "verde = lo que la arquitectura permite · azul = lo que cada compuerta exige")
    ax1.grid(axis="x", alpha=0.3, which="both")
    # loc="lower left": arriba a la izquierda se montaba sobre la etiqueta verde.
    ax1.legend(handles=[
        Patch(facecolor=VERDE, alpha=0.7, label="existe (G3–G10)"),
        Patch(facecolor=AZUL, alpha=0.55, label="permitido por la compuerta"),
    ], loc="lower left", fontsize=8)

    # ---- PANEL 2: el zoom del techo, escala lineal --------------------------
    ax2.axvspan(lo_b, hi_b, facecolor=GRIS, alpha=0.18, hatch="\\\\",
                edgecolor="#bdbdbd", linewidth=0.4,
                label=f"±5 % de Prosa sobre G10  ({lo_b:.2f} – {hi_b:.2f} nm)")
    ax2.axvline(DENDRIMERO_TECHO_MEDIDO_nm, color=VERDE, ls="-", lw=2.0,
                label=f"G10 medido, Prosa 2001: {DENDRIMERO_TECHO_MEDIDO_nm:.2f} nm")
    ax2.axvline(techo, color=VERDE, ls="--", lw=2.0,
                label=f"G10 calculado, Maiti 2004: {techo:.2f} nm  (techo usado)")
    ax2.axvline(u["envolvimiento"], color=ROJO, ls="-", lw=2.2,
                label=f"envolvimiento exige: {u['envolvimiento']:.3f} nm")

    ax2.annotate("", xy=(techo, 0.62), xytext=(u["envolvimiento"], 0.62),
                 arrowprops=dict(arrowstyle="<->", color="black", lw=1.3))
    ax2.text((techo + u["envolvimiento"]) / 2, 0.66,
             f"{u['envolvimiento'] - techo:.2f} nm", ha="center", va="bottom",
             fontsize=10, fontweight="bold")
    # El aviso va centrado en el EJE, no entre las dos líneas: ahí cruzaba la
    # línea roja y se leía partido.
    ax2.text(0.5, 0.16, "EMPATE TÉCNICO — la banda de incertidumbre de la medida "
             f"llega a {hi_b:.2f} nm\ny cruza el umbral de envolvimiento. "
             "Firme en valores centrales, marginal en el extremo.",
             transform=ax2.transAxes, ha="center", va="center", fontsize=9,
             color=ROJO,
             bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                       edgecolor=ROJO, alpha=0.92, linewidth=0.8))

    ax2.set_xlim(13.0, 15.1)
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("diámetro (nm)")
    ax2.set_title("Zoom del techo: por qué el resultado NO es un caso cerrado")
    ax2.grid(axis="x", alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def figura_dendrimero_generaciones(nombre="dendrimero_generaciones.png"):
    """Diámetro por generación: lo medido frente a lo calculado, y los umbrales."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    u = _umbrales_modelo()
    # Maiti et al. 2004, Tabla 4, columna R_SAV (Å) -> diámetro en nm.
    maiti = {1: 1.79, 2: 2.30, 3: 2.78, 4: 3.55, 5: 4.63, 6: 5.77,
             7: 7.29, 8: 9.18, 9: 11.49, 10: 14.16, 11: 17.62}
    VERDE, ROJO, AZUL, GRIS = "#2e7d32", "#c62828", "#1565c0", "#9e9e9e"

    fig, ax = plt.subplots(figsize=(11, 6.4))
    gp = sorted(DENDRIMERO_GENERACIONES_nm)
    dp = [DENDRIMERO_GENERACIONES_nm[g] for g in gp]
    ax.errorbar(gp, dp, yerr=[d * DENDRIMERO_PRECISION for d in dp],
                fmt="o-", color=VERDE, lw=2, ms=7, capsize=4,
                label="MEDIDO · Prosa 2001, SAXS en metanol (±5 %)")
    gm = sorted(maiti)
    ax.plot(gm, [maiti[g] for g in gm], "s--", color="#6a1b9a", lw=1.6, ms=5,
            alpha=0.85, label="CALCULADO · Maiti 2004, MD en fase gas")

    ax.axhspan(0, u["glicocalix"], color=AZUL, alpha=0.10)
    ax.axhline(u["glicocalix"], color=AZUL, lw=1.8)
    ax.text(1.1, u["glicocalix"] - 0.35, f"tamiz del glicocálix: {u['glicocalix']:.2f} nm",
            color=AZUL, fontsize=8.5, va="top")
    ax.axhline(u["envolvimiento"], color=ROJO, lw=1.8)
    ax.text(1.1, u["envolvimiento"] + 0.25,
            f"envolvimiento: {u['envolvimiento']:.2f} nm", color=ROJO,
            fontsize=8.5, va="bottom")

    # facecolor, no color: con 'color' matplotlib avisa de que pisa edgecolor.
    ax.axvspan(10.5, 11.5, facecolor=GRIS, alpha=0.28, hatch="\\\\",
               edgecolor="#bdbdbd", linewidth=0.4)
    ax.text(11, 1.0, "G11\nno sintetizable\n(tensión estérica)", ha="center",
            fontsize=8, color="#555555")
    ax.axvline(10.5, color="black", ls=":", lw=1.4)
    ax.text(10.42, 17.0, "techo: G10 es la última generación completable",
            rotation=90, ha="right", va="top", fontsize=8.5)

    ax.set_xticks(range(1, 12))
    ax.set_xticklabels([f"G{g}" for g in range(1, 12)])
    ax.set_xlabel("generación")
    ax.set_ylabel("diámetro (nm)")
    ax.set_ylim(0, 19)
    ax.set_title("Diámetro del PAMAM por generación\n"
                 "las dos fuentes divergen abajo y coinciden en G10, que es "
                 "donde se decide el resultado")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def figura_dendrimero_carga(nombre="dendrimero_carga.png"):
    """Constante de complejación por generación. El dato de la carga útil."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as _np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
    cols = {4: "#8d6e63", 7: "#2e7d32", 10: "#1565c0"}

    for ax, sup, titulo_ in ((ax1, "amina", "superficie de aminas (–NH₂)"),
                             (ax2, "ester", "superficie de éster (–COOCH₃)")):
        gens = sorted({g for ph in DEVARAKONDA_K[sup].values() for g in ph})
        x = _np.arange(len(gens))
        anchos = 0.26
        for k, ph in enumerate((4, 7, 10)):
            v = [DEVARAKONDA_K[sup][ph].get(g, (0, 0))[0] for g in gens]
            e = [DEVARAKONDA_K[sup][ph].get(g, (0, 0))[1] for g in gens]
            ax.bar(x + (k - 1) * anchos, v, anchos, yerr=e, capsize=3,
                   color=cols[ph], alpha=0.85, label=f"pH {ph}")
        ax.set_xticks(x)
        ax.set_xticklabels([f"G{g:g}" for g in gens])
        ax.set_title(titulo_, fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8)

    ax1.set_ylabel("K(1:1)  (M⁻¹)")
    ax1.text(0.02, 0.96, "G0 no aparece:\nno formó complejo medible",
             transform=ax1.transAxes, fontsize=8, va="top", color="#c62828")
    fig.suptitle("Carga útil del dendrímero — Devarakonda 2004, Tabla 2 "
                 "(nifedipino, no fingolimod)\n"
                 "Todos los perfiles son A_L con pendiente < 1: la estequiometría "
                 "es 1:1, UNA molécula de fármaco por dendrímero", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def figura_dendrimero_vs_liposoma(nombre="dendrimero_vs_liposoma.png"):
    """Las dos clases contra las dos ventanas. Ninguna cabe en las dos."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    u = _umbrales_modelo()
    d_min_lipo = G.diametro_liposoma_minimo_nm(4.0, 4.0)
    VERDE, ROJO, AZUL, MORADO = "#2e7d32", "#c62828", "#1565c0", "#6a1b9a"

    fig, ax = plt.subplots(figsize=(12.5, 5.0))
    ax.axvspan(1, u["glicocalix"], color=AZUL, alpha=0.13)
    ax.axvspan(u["envolvimiento"], 2000, color=ROJO, alpha=0.10)
    ax.text(2.4, 2.62, f"ventana del glicocálix\nd ≤ {u['glicocalix']:.2f} nm",
            color=AZUL, fontsize=9, ha="center")
    ax.text(120, 2.62, f"ventana del envolvimiento\nd ≥ {u['envolvimiento']:.2f} nm",
            color=ROJO, fontsize=9, ha="center")

    ax.barh(1, DENDRIMERO_TECHO_nm - DENDRIMERO_SUELO_nm, left=DENDRIMERO_SUELO_nm,
            height=0.42, color=MORADO, alpha=0.75, edgecolor="black", linewidth=0.6)
    ax.text(DENDRIMERO_SUELO_nm * 0.93, 1, f"{DENDRIMERO_SUELO_nm:.1f}", fontsize=8,
            ha="right", va="center")
    ax.text(DENDRIMERO_TECHO_nm * 1.07, 1, f"{DENDRIMERO_TECHO_nm:.1f}", fontsize=8,
            va="center")
    ax.barh(0, 2000 - d_min_lipo, left=d_min_lipo, height=0.42, color=VERDE,
            alpha=0.7, edgecolor="black", linewidth=0.6)
    ax.text(d_min_lipo * 0.93, 0, f"{d_min_lipo:.0f}", fontsize=8, ha="right",
            va="center")

    ax.set_xscale("log")
    ax.set_xlim(1, 2000)
    ax.set_ylim(-0.8, 3.0)
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["dendrímero PAMAM\n(G3 – G10)",
                        "liposoma\n(bicapa cerrada)"], fontsize=9.5)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("diámetro (nm)")
    ax.set_title("Por qué ninguna clase cabe en las dos ventanas\n"
                 "el dendrímero solo alcanza la de la izquierda · "
                 "el liposoma solo la de la derecha")
    ax.grid(axis="x", alpha=0.3, which="both")
    ax.legend(handles=[
        Patch(facecolor=MORADO, alpha=0.75, label="tamaños que la arquitectura permite"),
        Patch(facecolor=AZUL, alpha=0.13, label="exige el glicocálix"),
        Patch(facecolor=ROJO, alpha=0.10, label="exige el envolvimiento"),
    ], loc="upper right", bbox_to_anchor=(1.0, 0.78), fontsize=8)
    fig.tight_layout()
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def informe_polimero():
    """Resultado de G.1b en texto, sin prosa."""
    u = _umbrales_modelo()
    print("=" * 79)
    print(" POLÍMERO MACIZO · límite geométrico (tarea G.1b)")
    print("=" * 79)
    print("  POLÍMERO MACIZO — suelo del glóbulo de cadena colapsada")
    print("    d = (6M / pi rho N_A)^(1/3)   ·   sin techo arquitectónico")
    print()
    print("    POLÍMERO       Mw (kDa)   rho (g/cm3)   suelo (nm)   ¿pasa glicocálix?")
    print("    " + "-" * 68)
    for etq, mw in (("PLA15", 15), ("PLA24", 24), ("PLGA 85:15", 53), ("PLA60", 60)):
        rho = POLIMERO_DENSIDAD_g_cm3[etq]
        d = diametro_globulo_colapsado_nm(mw * 1000.0, rho)
        print(f"    {etq:<13} {mw:6.0f}     {rho:7.2f}      {d:7.2f}      "
              f"{'sí' if d <= u['glicocalix'] else 'no':^12s}")
    print()
    a = diametro_globulo_colapsado_nm(53000.0, 1.19)
    for n in (1, 10):
        b = diametro_globulo_colapsado_nm(53000.0 + n * 307.48, 1.19)
        print(f"    sensibilidad al fármaco: {n:2d} molécula(s) de fingolimod "
              f"-> {b:.4f} nm ({100*(b/a-1):+.2f} %)")
    print()
    print("  UMBRALES DEL MODELO (leídos del código)")
    print(f"    tamiz del glicocálix        d <= {u['glicocalix']:7.3f} nm")
    print(f"    envolvimiento de membrana   d >= {u['envolvimiento']:7.3f} nm")
    print()
    print("  RESULTADO  el polímero macizo es la única clase sin techo y con suelo")
    print("             por debajo del glicocálix.")
    print()
    print("  Salvedades y fuentes: verificacion/verificacion_polimero_micela_tarea_G_1b.md")
    print("=" * 79)


def figura_polimero(nombre="polimero_suelo.png"):
    """Suelo del glóbulo colapsado frente a la masa molar, con las 4 densidades."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    u = _umbrales_modelo()
    VERDE, ROJO, AZUL = "#2e7d32", "#c62828", "#1565c0"
    mws = np.linspace(1.0, 500.0, 400)

    fig, ax = plt.subplots(figsize=(11, 6.2))
    for etq, rho in POLIMERO_DENSIDAD_g_cm3.items():
        ax.plot(mws, [diametro_globulo_colapsado_nm(m * 1000.0, rho) for m in mws],
                lw=1.8, label=f"{etq}  (ρ = {rho:.2f} g/cm³)")

    ax.axhspan(0, u["glicocalix"], color=AZUL, alpha=0.10)
    ax.axhline(u["glicocalix"], color=AZUL, lw=1.8)
    ax.text(1.5, u["glicocalix"] - 0.4, f"tamiz del glicocálix: {u['glicocalix']:.2f} nm",
            color=AZUL, fontsize=8.5, va="top")
    ax.axhline(u["envolvimiento"], color=ROJO, lw=1.8)
    # A la izquierda lo tapaba la leyenda. Se pone a la derecha, donde no hay
    # ni curvas ni puntos a esa altura.
    ax.text(480, u["envolvimiento"] + 0.3,
            f"envolvimiento: {u['envolvimiento']:.2f} nm", color=ROJO,
            fontsize=8.5, va="bottom", ha="right")

    # PLGA 53 kDa y PLA60 caen casi encima (5.21 y 5.38 nm) y sus etiquetas se
    # pisaban hasta quedar ilegibles. Se separan a mano, una arriba y otra
    # abajo, con línea de guía.
    for m, rho, txt, dx, dy in ((15, 1.14, "PLA15", 8, -22),
                                (53, 1.19, "PLGA 53 kDa", -66, -30),
                                (60, 1.22, "PLA60", 14, 16)):
        d = diametro_globulo_colapsado_nm(m * 1000.0, rho)
        ax.plot([m], [d], "o", color="black", ms=6, zorder=5)
        ax.annotate(f"{txt}\n{d:.2f} nm", (m, d), textcoords="offset points",
                    xytext=(dx, dy), fontsize=8, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="#bbbbbb", lw=0.6, alpha=0.9),
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="#888888"))

    ax.set_xscale("log")
    ax.set_xlim(1, 500)
    ax.set_ylim(0, 18)
    ax.set_xlabel("masa molar de la cadena (kDa)")
    ax.set_ylabel("diámetro mínimo (nm)")
    ax.set_title("Polímero macizo: suelo de una sola cadena colapsada\n"
                 "d = (6M/πρN_A)^(1/3) · densidades medidas por Parker 2010, Tabla 2")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def figura_clases(nombre="clases_ventanas.png"):
    """Las cuatro clases contra las dos ventanas, en un solo eje."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    u = _umbrales_modelo()
    AZUL, ROJO, GRIS = "#1565c0", "#c62828", "#9e9e9e"
    d_lipo = G.diametro_liposoma_minimo_nm(4.0, 4.0)
    d_poli = diametro_globulo_colapsado_nm(53000.0, POLIMERO_DENSIDAD_POR_DEFECTO)

    filas = [
        ("liposoma\n(bicapa cerrada)", d_lipo, 2000.0, "#2e7d32"),
        ("dendrímero PAMAM\n(G3 – G10)", DENDRIMERO_SUELO_nm, DENDRIMERO_TECHO_nm, "#6a1b9a"),
        ("polímero macizo\n(PLGA 53 kDa)", d_poli, 2000.0, "#ef6c00"),
    ]

    fig, ax = plt.subplots(figsize=(12.5, 5.4))
    ax.axvspan(1, u["glicocalix"], color=AZUL, alpha=0.13)
    ax.axvspan(u["envolvimiento"], 2000, color=ROJO, alpha=0.10)
    ax.text(2.6, 3.62, f"ventana del glicocálix\nd ≤ {u['glicocalix']:.2f} nm",
            color=AZUL, fontsize=9, ha="center")
    ax.text(140, 3.62, f"ventana del envolvimiento\nd ≥ {u['envolvimiento']:.2f} nm",
            color=ROJO, fontsize=9, ha="center")
    ax.axvspan(u["glicocalix"], u["envolvimiento"], facecolor=GRIS, alpha=0.30,
               hatch="\\\\", edgecolor="#bdbdbd", linewidth=0.4)

    for i, (etq, lo, hi, col) in enumerate(filas):
        y = len(filas) - 1 - i
        ax.barh(y, hi - lo, left=lo, height=0.44, color=col, alpha=0.75,
                edgecolor="black", linewidth=0.6)
        ax.text(lo * 0.92, y, f"{lo:.1f}", fontsize=8, ha="right", va="center")
        if hi < 1000:
            ax.text(hi * 1.08, y, f"{hi:.1f}", fontsize=8, va="center")

    ax.set_xscale("log")
    ax.set_xlim(1, 2000)
    ax.set_ylim(-0.8, 4.0)
    ax.set_yticks(range(len(filas) - 1, -1, -1))
    ax.set_yticklabels([e for e, _, _, _ in filas], fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("diámetro (nm)")
    # OJO con el título: decir "ninguna clase cruza el hueco" sería FALSO — el
    # liposoma y el polímero macizo lo cruzan de sobra. Lo que no existe es un
    # DIÁMETRO que esté en las dos ventanas, porque el hueco las separa.
    ax.set_title("Las tres clases contra las dos ventanas\n"
                 "barra = tamaños que la arquitectura permite · la franja rayada "
                 "no pertenece a ninguna ventana")
    ax.grid(axis="x", alpha=0.3, which="both")
    ax.legend(handles=[
        Patch(facecolor=AZUL, alpha=0.13, label="exige el glicocálix"),
        Patch(facecolor=GRIS, alpha=0.30, hatch="\\\\", label="hueco entre ventanas"),
        Patch(facecolor=ROJO, alpha=0.10, label="exige el envolvimiento"),
    ], loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(nombre, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return nombre


def figuras_dendrimero():
    """Las seis figuras del dendrímero. No toca ninguna de las de liposoma."""
    import contextlib as _c
    import io as _io

    figura_dendrimero()
    figura_dendrimero_generaciones()
    figura_dendrimero_carga()
    figura_dendrimero_vs_liposoma()
    # El set común de tres, con el catálogo de dendrímeros. Se silencia su print
    # para que la lista de abajo salga una sola vez y completa.
    with _c.redirect_stdout(_io.StringIO()):
        figuras(prefijo="dend_reales", catalogo=catalogo_dendrimero())

    hechas = [
        ("dend_reales_ventanas.png", "ventanas de tamaño de cada compuerta"),
        ("dend_reales_matriz.png", "qué generación puede usar qué ruta"),
        ("dend_reales_recorrido.png", "dónde se cae cada generación, paso a paso"),
        ("dendrimero_ventana.png", "ventana geométrica + zoom del techo"),
        ("dendrimero_generaciones.png", "diámetro por generación, medido vs calculado"),
        ("dendrimero_vs_liposoma.png", "las dos clases contra las dos ventanas"),
        ("dendrimero_carga.png", "carga útil: K(1:1) por generación y pH"),
    ]
    print("\n  Figuras guardadas:")
    for f, q in hechas:
        print(f"    {f:32s} {q}")
    print(f"\n  en {__import__('os').getcwd()}")
    return [f for f, _ in hechas]


def figuras_liposoma_separadas():
    """Las tres figuras de rutas, una tanda por cada mitad del catálogo.

    Las de prefijo `rutas_` siguen existiendo con el catálogo COMPLETO (no se
    tocan: las usan el informe y quien ya las cite). Estas dos tandas son para
    la web, que separa liposomas reales de liposomas teóricos y no puede
    enseñar una figura que mezcle los dos.
    """
    import contextlib as _c
    import io as _io

    hechas = []
    for prefijo, catalogo in (("lip_reales", CATALOGO_REAL),
                              ("lip_teoricos", CATALOGO_TEORICO)):
        with _c.redirect_stdout(_io.StringIO()):
            figuras(prefijo=prefijo, catalogo=catalogo, incluir_ventanas=True)
        hechas += [f"{prefijo}_ventanas.png", f"{prefijo}_matriz.png",
                   f"{prefijo}_recorrido.png"]

    print("\n  Figuras del liposoma separadas (reales / teóricos):")
    for f in hechas:
        print(f"    {f}")
    return hechas


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]

    # --callado omite lo que ya se ha impreso antes. Lo usa correr.sh, que
    # encadena varios pasos y si no repetiría las pruebas y la tabla.
    callado = "--callado" in args

    if "--dendrimero" in args:
        informe_dendrimero()
        figuras_dendrimero()
    elif "--teoricos" in args:
        informe_dendrimeros_teoricos()
        figuras_dendrimeros_teoricos()
        informe_polimeros_teoricos()
        figuras_polimeros_teoricos()
    elif "--polimero" in args:
        informe_polimero()
        print(f"\n  Figuras guardadas:")
        print(f"    {figura_polimero():32s} suelo frente a masa molar")
        print(f"    {figura_clases():32s} las cuatro clases contra las ventanas")
    elif "--notas" in args:
        notas()
    elif "--detalle" in args:
        if not callado:
            validar_contra_experimentos()
            print()
        informe(detalle=True, tabla=not callado)
    elif "--figuras" in args:
        if not callado:
            resumen()
        figuras()
        figuras_liposoma_separadas()
    elif "--tests" in args:
        validar_contra_experimentos()
    else:
        resumen(pistas=not callado)
