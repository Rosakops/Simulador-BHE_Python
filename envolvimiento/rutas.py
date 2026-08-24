#!/usr/bin/env python3
# =============================================================================
#  SIMULADOR DE RUTAS COMPLETAS · proyecto BHE / SENACYT
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

from dataclasses import dataclass
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

def g_transportador_fabricable(d: Diseno):
    """¿Es geométricamente posible un liposoma de ese tamaño?

    El suelo del liposoma sale de que es una BICAPA CERRADA
    (d_ext = d_núcleo + 2·t_bicapa, con t medido en 3.6–4.6 nm), y se exige
    un núcleo acuoso mínimo de 4 nm para que encapsule algo.

    El simulador es exclusivo para liposomas (decisión de Jhovan,
    2026-08-24): cualquier otra clase devuelve DESCONOCIDA en vez de
    aplicarle el límite del liposoma.
    """
    if d.clase == "liposoma":
        minimo = G.diametro_liposoma_minimo_nm(4.0, 4.0)
        return _cmp("Transportador fabricable", d.diametro_nm, minimo, "nm",
                    "Pan et al. 2008, PRL 100:198103, Fig. 3c", True)

    return Resultado("Transportador fabricable", DESCONOCIDA,
                     d.diametro_nm, None, "nm", None, "",
                     f"clase '{d.clase}' desconocida para esta compuerta")

def g_glicocalix_tamiz(d: Diseno):
    """¿Atraviesa el tamiz de la matriz de fibras del glicocálix?

    CAMBIADO el 2026-08-18: por encima del poro geométrico (Weinbaum et al.
    2003) esta compuerta YA NO devuelve FALLA, devuelve DESCONOCIDA. Motivo:
    dos fuentes primarias con BHE real (no cultivo) muestran cruce medible,
    dependiente de carga, de partículas mucho más grandes que el poro de
    9nm:
      - Lockman et al. 2004, J Drug Target 12(9-10):635-41, DOI
        10.1080/10611860400015936 (perfusión cerebral in situ en rata, BHE
        nativa): nanopartículas de 74.7-127.1nm SÍ cruzan, con permeabilidad
        dependiente de ζ (aniónica 7.93e-3 > neutra 4.10e-3 > catiónica
        2.73e-3 ml/s/g).
      - Gromnicova et al. 2016, PLoS ONE 11(8):e0161610, PMC4999129
        (hCMEC/D3): remover el glicocálix con enzimas NO cambió la captación
        de nanopartículas en endotelio cerebral.
    La exclusión binaria por tamaño de poro no predice el desenlace real a
    esta escala, probablemente porque el cruce real es vía vesicular, que
    no requiere pasar por el poro. El dato geométrico (Weinbaum 2003) sigue
    siendo real y se reporta en `valor`/`umbral`; deja de ser el veredicto
    final por sí solo. Por DEBAJO del poro sigue dando PASA sin cambios.
    """
    dmax = 2.0 * G.radio_exclusion_nm()
    r = _cmp("Tamiz del glicocálix", d.diametro_nm, dmax, "nm",
             "Weinbaum et al. 2003, PNAS 100:7988, Transport Model", False)
    if r.estado == FALLA:
        r.estado = DESCONOCIDA
        r.motivo = (f"{d.diametro_nm:.1f}nm supera el poro geométrico "
                    f"({dmax:.2f}nm, Weinbaum 2003), pero Lockman 2004 (in "
                    "situ, BHE nativa) y Gromnicova 2016 (hCMEC/D3) muestran "
                    "cruce real medible por encima de ese tamaño; el poro "
                    "no es excluyente por sí solo a esta escala")
    return r


_F_KABEDEV = ("Kabedev & Lobaskin 2022, Nanomedicine 17:979, Fig.4B/Fig.6 "
              "(digitalizado) · ajuste kT_hombro() en glicocalix.py")

#  UMBRALES DE VEREDICTO: decisión METODOLÓGICA del proyecto, NO un dato de
#  Kabedev & Lobaskin 2022 (esa fuente no publica ningún corte pasa/no-pasa).
#  Se toma prestada la convención ESTÁNDAR de ciencia de coloides sobre
#  cuándo una barrera de energía es cinéticamente insuperable, aplicada por
#  analogía (mismo marco de Boltzmann/DLVO, distinto sistema: partícula-
#  partícula en la fuente original, partícula-malla aquí):
#
#    Tadros T. (2007), "General Principles of Colloid Stability and the Role
#    of Surface Forces", cap. 1 de Colloid Stability: The Role of Surface
#    Forces, Part I (Colloids and Interface Science Series vol. 1),
#    Wiley-VCH, ISBN 978-3-527-31462-1. Cita a su vez a Verwey EJW & Overbeek
#    JTG (1948), Theory of Stability of Lyophobic Colloids, Elsevier.
#    "The condition for colloid stability is to have an energy maximum
#    (barrier) that is much larger than the thermal energy... In general,
#    one requires Gmax > 25kT." El mismo capítulo describe la floculación
#    débil y reversible con un mínimo secundario de solo "a few kT units"
#    (1-5 kT), insuficiente para bloquear nada.
#
#  De ahí las dos bandas: <5kT barrera insignificante frente a la agitación
#  térmica (PASA); >25kT convención estándar de barrera insuperable (FALLA);
#  entre medio, zona de transición sin criterio limpio (DESCONOCIDA, mismo
#  patrón que la "zona gris" ya usada en g_difusion_ecs()).
UMBRAL_PMF_BAJO_kT = 5.0     # por debajo: barrera insignificante frente a kT
UMBRAL_PMF_ALTO_kT = 25.0    # por encima: convención DLVO de barrera insuperable
_F_DLVO_UMBRAL = ("Tadros 2007, cap.1 de Colloid Stability: Role of Surface "
                   "Forces Part I (Wiley-VCH) · Verwey & Overbeek 1948: "
                   "convención Gmax>25kT, aplicada por analogía al glicocálix, "
                   "NO es un umbral publicado para este sistema")


def g_glicocalix_pmf(d: Diseno):
    """Barrera de energía (kT) del glicocálix en la meseta/hombro, vía el
    ajuste kT_hombro() (Kabedev & Lobaskin 2022). Complementa (NO sustituye)
    a g_glicocalix_tamiz(): esa compuerta tiene un umbral MEDIDO (el hueco de
    9nm de Weinbaum, verificado a 3 decimales) y sigue siendo el criterio
    geométrico independiente.

    El veredicto PASA/FALLA/DESCONOCIDA de ESTA compuerta usa un umbral que
    NO viene de Kabedev; es la convención Gmax>25kT de ciencia de coloides
    (Tadros 2007 / Verwey-Overbeek 1948), aplicada por analogía. Declarado
    así en `fuente` y en `motivo` para que nunca se lea como si fuera un
    dato de la fuente primaria del PMF.
    """
    z = d.zeta_mV
    if z > 0.0:
        return Resultado("Barrera del glicocálix (kT)", DESCONOCIDA,
                         None, None, "kT", None, _F_KABEDEV,
                         motivo="ζ positivo: la carga positiva NO está modelada "
                                "en kT_hombro() (no es un múltiplo del neutro; "
                                "en los datos digitalizados a veces casi cancela "
                                "la barrera y a veces forma un pozo atractivo). "
                                "Ver PMF_HOMBRO_kT_A/_B en glicocalix.py para "
                                "los números brutos por radio y panel.")

    carga = "neutro" if z == 0.0 else "negativo"
    r = G.kT_hombro(d.radio_nm, carga)
    kt = r["kT"]
    extra = (", EXTRAPOLADO fuera del rango medido (3.5-10nm), menos "
             "confiable" if r["extrapolado"] else ", dentro del rango medido")
    fuente = _F_KABEDEV + " · " + _F_DLVO_UMBRAL
    base_motivo = f"carga {carga}, R={d.radio_nm:.1f}nm{extra}."

    if kt < UMBRAL_PMF_BAJO_kT:
        return Resultado("Barrera del glicocálix (kT)", PASA,
                         kt, UMBRAL_PMF_BAJO_kT, "kT", UMBRAL_PMF_BAJO_kT - kt,
                         fuente,
                         advertencia=f"{base_motivo} Umbral de PASA (<5kT) es "
                                     "convención propia, no dato de Kabedev; "
                                     "ver g_glicocalix_tamiz() para el criterio "
                                     "geométrico medido.")
    if kt > UMBRAL_PMF_ALTO_kT:
        return Resultado("Barrera del glicocálix (kT)", FALLA,
                         kt, UMBRAL_PMF_ALTO_kT, "kT", UMBRAL_PMF_ALTO_kT - kt,
                         fuente,
                         motivo=f"{base_motivo} Supera la convención DLVO de "
                                "barrera cinéticamente insuperable (Gmax>25kT, "
                                "Tadros 2007/Verwey-Overbeek 1948); no es un "
                                "umbral publicado para el glicocálix.")
    return Resultado("Barrera del glicocálix (kT)", DESCONOCIDA,
                     kt, None, "kT", None, fuente,
                     motivo=f"{base_motivo} Cae en la zona de transición "
                            f"({UMBRAL_PMF_BAJO_kT:.0f}-{UMBRAL_PMF_ALTO_kT:.0f}"
                            "kT) entre la convención de barrera insignificante "
                            "y la de barrera insuperable; ninguna de las dos "
                            "aplica con confianza aquí.")


def g_envolvimiento(d: Diseno, kappa_kT=25.0, sigma_mNm=0.03, hamaker_J=4.5e-21):
    """¿Es lo bastante grande para que la membrana lo envuelva?"""
    w = E.w_adhesion(d.radio_nm, d.zeta_mV, d.peg_nm, hamaker_J)
    d_min = 2.0 * E.radio_critico_nm(w, kappa_kT)
    return _cmp("Envolvimiento de membrana", d.diametro_nm, d_min, "nm",
                "Deserno 2004, PRE 69:031903, Sec. III C", True)


def g_caveola(d: Diseno):
    """¿Cabe en una caveola?

    CAMBIADO el 2026-08-18: por encima del diámetro de caveola esta
    compuerta YA NO devuelve FALLA, devuelve DESCONOCIDA. Motivo: Wang et
    al. 2026, J Nanobiotechnology 24:164, DOI 10.1186/s12951-026-04023-y,
    PMC12903337 (hCMEC/D3 Transwell) muestra que NP-120 (119.5nm, muy por
    encima de los 60-80nm de caveola) cruza igual, por una vía confirmada
    farmacológicamente como mediada por clatrina/dinamina (clorpromazina y
    dinasor reducen el transporte sin afectar TEER ni la permeabilidad de
    trazador), NO caveolar. El simulador no modela la vía de clatrina;
    decir FALLA aquí implicaría que ninguna vía vesicular es posible, lo
    cual la propia fuente contradice. Por DEBAJO del tamaño de caveola sigue
    dando PASA sin cambios.
    """
    r = _cmp("Compuerta de caveola", d.diametro_nm, E.DIAM_CAVEOLA_MAX_nm, "nm",
             "Bastiani & Parton 2010, J Cell Sci 123:3831", False)
    if r.estado == FALLA:
        r.estado = DESCONOCIDA
        r.motivo = (f"{d.diametro_nm:.1f}nm no cabe en una caveola "
                    f"({E.DIAM_CAVEOLA_MAX_nm:.1f}nm), pero Wang 2026 "
                    "muestra cruce real vía clatrina/dinamina a este tamaño; "
                    "el simulador no modela esa vía, no se puede excluir "
                    "por esta compuerta sola")
    return r


# -----------------------------------------------------------------------------
#  DIFUSIÓN EN EL ESPACIO EXTRACELULAR: compuerta de DOS VARIABLES
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
#  que la compuerta pregunta (infusión IV de monocitos cargados con
#  nanopartícula y recuento en cerebro inflamado) y da un PICO a las 48 h, con
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

    Tiempo de tránsito · Tong et al. 2016, S1 Text, textual: "Brain tissues
    were collected at 1, 2, 3 and 7 days following MDM transfers into LPS ICI
    treated mice. The number of recruited donor-derived cells peaked at 48h and
    decreased afterwards". Las células ya se detectan en el primer punto de
    muestreo (24 h), así que el tránsito real cae entre 24 y 48 h.

    Tiempo de descarga · Mao et al. 2014, textual: "The release over a time
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
    FITC tres tamaños, 150/550/700 nm; ese es el dato que define el umbral.
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
#   (1) CARGA DEL LIPOSOMA: RESUELTO. Mouzoura et al. 2025 (Int J Nanomedicine
#       20:239, doi 10.2147/IJN.S494512) da relación molar fármaco:lípido 1:8
#       con eficiencia de carga 94-97.2 %, con FTY720 directo. No es la
#       eficiencia de encapsulación de Mao 2014, que era lo que había antes y
#       no servía.
#   (2) UMBRAL EN PARÉNQUIMA: RESUELTO. Foster et al. 2007 (JPET 323(2):469,
#       Tabla 3) mide 398 ± 186 ng/g de FTY720-P en cerebro a una dosis que FUE
#       terapéuticamente eficaz en EAE (0.3 mg/kg, días 11-33). Eso da la
#       concentración diana contra la que comparar.
#   (3) CARGA DEL DENDRÍMERO CON FINGOLIMOD: sigue sin existir. El 1:1 de
#       Devarakonda 2004 es con nifedipino y es estequiometría de complejo
#       inferida de una pendiente, no carga útil medida. Congelado: el
#       dendrímero quedó fuera del foco de simulación el 2026-08-17.
#
# LO QUE FALTA AHORA es UN solo eslabón, no tres: qué FRACCIÓN de la dosis
# inyectada llega al parénquima. Y ese dato NO es un hueco independiente, es
# exactamente el mismo que bloquea la transcitosis. Decirlo importa: el modelo
# tiene menos agujeros independientes de los que aparentaba, y cerrar la
# transcitosis cerraría dos compuertas, no una.
#
# Sigue DESCONOCIDA porque sin la fracción entregada no se puede decidir
# "suficiente". No se fuerza a PASA teniendo dos de tres números: eso sería
# fabricar el eslabón que falta.
g_carga_util = _sin_dato(
    "Carga útil suficiente",
    "carga del liposoma SÍ conocida (1:8 molar, 94-97 % de eficiencia, "
    "Mouzoura 2025) y umbral en parénquima SÍ conocido (398 ± 186 ng/g de "
    "FTY720-P a dosis eficaz en EAE, Foster 2007 Tabla 3). Falta UN eslabón: "
    "la fracción de dosis que llega al parénquima, que es el MISMO hueco de la "
    "transcitosis, no uno independiente. Para el dendrímero falta además la "
    "carga con fingolimod (el 1:1 de Devarakonda es con nifedipino)",
    "G.2")


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
        g_carga_util,
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
    "D · mediada por receptor": [
        g_transportador_fabricable,
        g_acceso_receptor,
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

    # -- CAMBIADO 2026-08-18: por encima del poro, DESCONOCIDA, no FALLA.
    #    Lockman 2004 (in situ, BHE nativa) y Gromnicova 2016 (hCMEC/D3)
    #    muestran cruce real por encima del poro de 9nm de Weinbaum.
    r_tamiz_grande = g_glicocalix_tamiz(Diseno("liposoma 150nm", 150.0, 0.0))
    chequeo("C5f tamiz del glicocálix >9nm da DESCONOCIDA, no FALLA",
            r_tamiz_grande.estado == DESCONOCIDA,
            f"da {r_tamiz_grande.estado} (Lockman 2004 / Gromnicova 2016)")

    # -- CAMBIADO 2026-08-18: por encima de caveola, DESCONOCIDA, no FALLA.
    #    Wang 2026 muestra cruce vía clatrina a 119.5nm, fuera de caveola.
    r_caveola_pequena = g_caveola(Diseno("furtivo 31nm", 31.0, 0.0))
    chequeo("C5g caveola <=80nm sigue dando PASA sin cambios",
            r_caveola_pequena.estado == PASA,
            f"da {r_caveola_pequena.estado}")

    r_caveola_grande = g_caveola(Diseno("liposoma 150nm", 150.0, 0.0))
    chequeo("C5h caveola >80nm da DESCONOCIDA, no FALLA",
            r_caveola_grande.estado == DESCONOCIDA,
            f"da {r_caveola_grande.estado} (Wang 2026, vía clatrina)")

    # -- g_glicocalix_pmf (Kabedev, 2026-08-18): SIEMPRE DESCONOCIDA, con o sin
    #    dato: es información continua sin umbral pasa/no-pasa medido, no un
    #    reemplazo de g_glicocalix_tamiz(). No se conecta a evaluar_ruta()
    #    (forzaría NO EVALUABLE en cualquier ruta que la incluyera) hasta que
    #    Jhovan decida cómo combinarla con el veredicto binario.
    r_pmf_pasa = g_glicocalix_pmf(Diseno("convencional D=20 (R=10)", 20.0, 0.0))
    chequeo("C5b g_glicocalix_pmf: R=10nm neutro (2.55kT) PASA el umbral bajo",
            r_pmf_pasa.estado == PASA and r_pmf_pasa.valor is not None,
            f"kT={r_pmf_pasa.valor:.2f} (umbral PASA <{UMBRAL_PMF_BAJO_kT:.0f}kT)")

    r_pmf_pos = g_glicocalix_pmf(Diseno("catiónico", 17.5, +6.7))
    chequeo("C5c g_glicocalix_pmf con ζ positivo da DESCONOCIDA SIN valor "
            "(carga positiva no modelada)",
            r_pmf_pos.estado == DESCONOCIDA and r_pmf_pos.valor is None)

    r_pmf_zonagris = g_glicocalix_pmf(Diseno("furtivo D=31 (R=15.5)", 31.0, -8.0))
    chequeo("C5d g_glicocalix_pmf: R=15.5nm negativo (22.4kT) cae en zona de "
            "transición y marca extrapolación",
            r_pmf_zonagris.estado == DESCONOCIDA
            and "EXTRAPOLADO" in r_pmf_zonagris.motivo,
            f"kT={r_pmf_zonagris.valor:.2f}, R=15.5nm")

    r_pmf_falla = g_glicocalix_pmf(Diseno("grande D=40 (R=20)", 40.0, 0.0))
    chequeo("C5e g_glicocalix_pmf: R=20nm neutro (31.0kT) FALLA la convención "
            "DLVO de barrera insuperable",
            r_pmf_falla.estado == FALLA,
            f"kT={r_pmf_falla.valor:.2f} (umbral FALLA >{UMBRAL_PMF_ALTO_kT:.0f}kT, "
            "Tadros 2007/Verwey-Overbeek 1948)")

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
    v, res_c = evaluar_ruta(Diseno("sigiloso", 50.0, -2.0),
                            RUTAS["C · adsortiva (carga positiva)"])
    chequeo("F1 una ruta con datos faltantes NO sale aprobada",
            v == "NO EVALUABLE" and not any(r.estado == FALLA for r in res_c),
            f"da '{v}'")

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

    # -- carga útil, añadida el 2026-08-13. Es DESCONOCIDA permanente mientras
    #    falten los tres números, igual que la transcitosis.
    chequeo("F26 la carga útil está en las CUATRO rutas",
            all(any(c(CATALOGO[0]).compuerta == "Carga útil suficiente"
                    for c in comps) for comps in RUTAS.values()))
    chequeo("F27 la carga útil sale DESCONOCIDA para todo el catálogo",
            all(g_carga_util(d).estado == DESCONOCIDA for d in CATALOGO),
            "(faltan carga con fingolimod, moléculas por liposoma y dosis)")

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
    print(" 4. CAMBIADO 2026-08-18: el tamiz del glicocálix y la compuerta de caveola")
    print("    ya NO excluyen por sí solos por encima de su umbral geométrico (Lockman")
    print("    2004, Gromnicova 2016 y Wang 2026 muestran cruce real por encima de esos")
    print("    tamaños). Los tres diseños teóricos ya no salen EXCLUIDOS de la ruta A por")
    print("    esa vía; quedan NO EVALUABLE por las incógnitas que sí siguen abiertas")
    print("    (transcitosis fuera de alcance, carga útil desconocida). Siguen EXCLUIDOS")
    print("    de la ruta B por captación fagocítica insuficiente (compuerta sin cambios).")
    print(" 5. La ruta B es además la única donde lo que tiene que difundir por el")
    print("    espacio extracelular es el FÁRMACO (1.68 nm) y no el transportador.")
    print("    Por eso Muselman 2026, de 700 nm, sale EXCLUIDO de C y D pero no de B:")
    print("    una partícula de 700 nm no atraviesa el parénquima, su carga sí.")
    print(" 6. LA CARGA ÚTIL es incógnita en las CUATRO rutas, y se añadió el")
    print("    2026-08-13. Hasta entonces el modelo preguntaba si el transportador")
    print("    LLEGA y nunca si entrega FÁRMACO SUFICIENTE: un diseño que llegara con")
    print("    una sola molécula habría salido igual de bien que uno con miles. Era un")
    print("    hueco INVISIBLE. Faltan tres números: carga de un dendrímero de")
    print("    generación alta CON FINGOLIMOD (el 1:1 de Devarakonda 2004 es con")
    print("    nifedipino y es un complejo de solubilización, no una carga útil),")
    print("    moléculas de fingolimod por liposoma (Mao 2014 da eficiencia de")
    print("    encapsulación, que no es lo mismo) y la dosis necesaria en parénquima.")
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

    El simulador es exclusivo para liposomas (decisión de Jhovan,
    2026-08-24): solo hay ventana verificada para clase="liposoma"; cualquier
    otra clase devuelve una ventana indefinida (barra rayada, TOPE 2000 nm).
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
    # 2026-08-19: con catálogos de un solo diseño (tarjetas de detalle por
    # liposoma) el eje sale muy bajo y el hueco de -0.155 no basta: la nota
    # de abajo queda pisada por la leyenda. -0.24 deja margen incluso con
    # un solo diseño; con catálogos más grandes el eje es más alto y el
    # mismo valor solo deja más aire, no rompe nada.
    ax.text(0.5, -0.24, "Cada '!' de un SÍ = una compuerta que pasa pero con "
            "salvedad declarada. SÍ!!! se apoya en tres.",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
            color="#555555")
    fig.tight_layout()
    fig.savefig(f"{prefijo}_matriz.png", dpi=160, bbox_inches="tight")
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


def _umbrales_modelo():
    """Los tres umbrales de tamaño, leídos del propio código, no a mano."""
    return dict(
        glicocalix=2.0 * G.radio_exclusion_nm(),
        envolvimiento=g_envolvimiento(Diseno("_", 20.0, 0.0)).umbral,
        fagocitosis=550.0,
    )
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

    if "--notas" in args:
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
