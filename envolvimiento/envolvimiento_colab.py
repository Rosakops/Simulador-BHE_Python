#@title Envolvimiento de membrana · BHE/SENACYT (Fase 2: tareas 2.4, 2.5, 2.6)
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

# =============================================================================
#  MÓDULO DE ENVOLVIMIENTO DE MEMBRANA: Fase 2 del proyecto BHE / SENACYT
#  Tareas 2.4, 2.5 y 2.6 del cronograma v3.
# =============================================================================
#
#  QUÉ HACE
#  --------
#  Toma la energía libre de interacción DLVO entre el nanotransportador y la
#  membrana endotelial, y responde: ¿la membrana puede envolverlo?
#
#  Encadena dos modelos:
#    (1) DLVO extendida  ->  G(D), energía de interacción esfera-plano.
#    (2) Deserno 2004    ->  criterio de envolvimiento a partir de la energía
#                            de adhesión por unidad de área, w.
#
#  El puente entre ambos es la aproximación de Derjaguin (ver más abajo, es el
#  punto delicado de todo el módulo).
#
#  QUÉ NO HACE
#  -----------
#  No predice permeación. Predice si se cumplen las condiciones ENERGÉTICAS
#  NECESARIAS para que se inicie la transcitosis por envolvimiento de membrana.
#  Un diseño que no las cumple queda EXCLUIDO por esa vía; uno que las cumple
#  es un candidato, no una predicción de éxito.
#
#  FUENTE PRINCIPAL
#  ----------------
#  Deserno M. "Elastic deformation of a fluid membrane upon colloid binding".
#  Phys. Rev. E 69, 031903 (2004). Preprint: arXiv:cond-mat/0303656.
#  Verificado íntegro. Ver verificacion/verificacion_wrapping_tarea_2_1.md
#
#  DOS DIFERENCIAS DELIBERADAS RESPECTO A cribado_bhe_tres_liposomas.py
#  --------------------------------------------------------------------
#  (a) Aquí NO se usa la regularización `D_sat` del van der Waals. En el cribado,
#      `_derivar()` calibra D_sat para que el pozo de vdW en el contacto valga
#      siempre 15 kT, y lo hace PARA CUALQUIER RADIO. Eso deja el vdW sin
#      dependencia de R (a R=100 nm el valor correcto es 7 veces más profundo),
#      lo cual es inaceptable aquí porque el radio es el eje principal de este
#      módulo. Se usa el corte estándar D0 = 0.165 nm (Israelachvili), que ya
#      impide la divergencia por sí solo.
#  (b) La salida es energía en kT y radios críticos, no fracción adherida. La
#      fracción satura al 100% y no discrimina (tarea 3.2 del cronograma).
#
# =============================================================================

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz

# =============================================================================
#  BLOQUE 1: PARÁMETROS
#  Todo lo ajustable vive AQUÍ. Nada de esto está enterrado en el código.
#  Cambiar un rango debe ser editar una línea y volver a ejecutar.
# =============================================================================

# ---- constantes físicas (SI) ----
EPS0 = 8.8541878128e-12   # F/m
KB   = 1.380649e-23       # J/K
QE   = 1.602176634e-19    # C
NM   = 1e-9               # m por nm

# ---- condiciones ----
TEMP_K = 310.15           # 37 °C
KT_J   = KB * TEMP_K

# ---- diana: barrera hematoencefálica (parámetros CON fuente verificada) ----
ZETA_BHE_mV = -11.4       # endotelio cerebral humano hCMEC/D3 a 37 °C.
                          # Santa-Maria et al. 2019, BBA Biomembranes 1861:1579, Fig. 4A.
                          # Corroborado: -12.7 mV en Kincses et al. 2020, Lab Chip 20:3792, Fig. 6F.
EPS_R    = 74.0           # permitividad relativa del agua a 37 °C
DEBYE_nm = 0.8            # longitud de Debye a ~150 mM y 37 °C
D0_nm    = 0.165          # separación de contacto (convención DLVO, Israelachvili)
D_MAX_nm = 20.0           # frontera lejana de integración

# ---- términos repulsivos de corto alcance (calibrados, NO medidos) ----
B_HID_kT = 1.5            # amplitud de la repulsión de hidratación
W_HID_nm = 0.3            # alcance de la hidratación
B_PEG_kT = 8.0            # amplitud de la repulsión estérica del cepillo de PEG

# ---- BARRIDOS DE IGNORANCIA -------------------------------------------------
# Estos NO son variables de diseño. Son magnitudes de valor fijo pero NO medido
# para endotelio cerebral. Van en la sección de LIMITACIONES del artículo, no en
# resultados. La conclusión solo es válida si sobrevive en los tres niveles.

KAPPA_kT = [15.0, 25.0, 50.0, 100.0]
# Rigidez de curvatura de la membrana ENDOTELIAL (no del liposoma: en este
# modelo el liposoma es el coloide rígido y la membrana es la que se deforma).
# No existe medida para endotelio cerebral. Anclajes:
#   15 -> extremo blando de bicapas PC fluidas medidas (Nagle 2017, Tabla 1)
#   25 -> K_C tilt-dependent de POPC (25.7) y SOPC (24.6), los mono-insaturados
#         más parecidos a una membrana plasmática. Nagle 2017, Chem Phys Lipids
#         205:18, Tabla 1. Coincide con los 20 kT que usa Deserno.
#   50 -> cubre diC22:1PC (46.2), el más rígido de los diez lípidos de Nagle.
#  100 -> AÑADIDO el 2026-08-13 (tareas C5 y C6). Cubre el efecto del colesterol
#         sobre lípidos de cadenas SATURADAS, que es el caso que el nivel de 50
#         decía cubrir y NO cubría. Ver abajo.
#
# DOS ERRORES CORREGIDOS EL 2026-08-13, los dos verificados sobre el PDF:
#
# (1) El comentario decía que el nivel de 50 cubre «el efecto del colesterol
#     sobre lípidos con cadenas saturadas (Pan 2008)». Es FALSO. Pan 2008
#     (PRL 100:198103) dice textualmente que K_C de DMPC saturada «increases
#     more than four-fold» con 30 % de colesterol. Con el K_Ctd de DMPC de
#     Nagle 2017 (24.6 kT) eso da MÁS DE 98 kT, no 50. El nivel de 50 solo
#     cubre diC22:1PC. De ahí el cuarto nivel.
#
# (2) El barrido MEZCLABA las dos columnas de la Tabla 1 de Nagle. El 25 y el
#     50 salen de K_Ctd (tilt-DEPENDIENTE): POPC 25.7, SOPC 24.6, diC22:1PC
#     46.2. El 15 sale de K_Cti (tilt-INDEPENDIENTE): DLPC 14.3, DMPC 15.6,
#     DOPC 16.3; en K_Ctd el mínimo es diPhyPC con 17.4. Las dos columnas
#     difieren por un factor de 1.17 a 1.63 según el lípido. Mezclarlas hace el
#     barrido MÁS ANCHO, o sea conservador, así que no invalida nada, pero
#     estaba mal justificado y ahora queda dicho.
#
# LO QUE ESTE BARRIDO NO PUEDE DECIDIR: Pan 2008 solo rigidiza con cadenas
# saturadas. Campbell et al. 2014 (Mol Pharm 11:3541, Tabla 1) da la
# composición del endotelio cerebral humano por clase de cabeza polar
# (esfingomielina 33.4 %, colesterol 20.8 %) pero NO el perfil de ácidos
# grasos, que es justo la variable que Pan necesita. Tarea C5, sigue abierta.
#
# Dispersión entre métodos para un mismo lípido: factor 2.3 (Nagle 2017, Tabla 3).

SIGMA_mNm = [0.003, 0.03, 0.3]
# Tensión lateral de la membrana. No existe medida para endotelio cerebral.
# Rango canónico de membranas plasmáticas: 0.003-0.3 mN/m (Morris & Homann 2001,
# J Membr Biol 179:79, citado a través de Shi & Baumgart 2015, Nat Commun 6:5974).
# AVISO: los métodos de tether miden tensión APARENTE = tensión de bicapa +
# adhesión membrana-citoesqueleto. Son COTAS SUPERIORES de la sigma de Helfrich.
# ESTE es el parámetro más consecuente del módulo: barre un factor 100.

HAMAKER_J = [3.0e-21, 4.5e-21, 6.5e-21]
# Sin valor primario verificado para esta composición lipídica. Se barre.

# ---- EJES DE DISEÑO ---------------------------------------------------------
# Estos SÍ son variables que el diseñador elige. El barrido ES el resultado y
# va en la sección de resultados, no en limitaciones.
#
# Rangos fijados en la tarea 0.2 a partir de formulaciones de fingolimod
# publicadas y verificadas (ver verificacion/verificacion_formulaciones_tarea_0_2.md):
#   · Gong et al. 2022, Nanophotonics 11:5133  -> 145 nm, zeta -28.33 mV, PEGilado
#   · Chow et al. 2025, Drug Deliv Transl Res 15:2022 -> 134 nm, zeta -0.24 mV
#     (rango del diseño de experimentos: 107 a 187 nm)
# El eje de radio llega hasta 125 nm (diámetro 250 nm) para cubrirlas todas con
# margen. El eje de zeta conserva la parte positiva porque es espacio de diseño
# legítimo, aunque no haya formulación publicada ahí.

EJE_R_nm    = np.linspace(5.0, 125.0, 121)   # radio externo (diámetro 10 a 250 nm)
EJE_ZETA_mV = np.linspace(-30.0, 30.0, 61)   # potencial zeta del nanotransportador
EJE_PEG_nm  = [0.0, 2.0, 5.0, 10.0]          # espesor del cepillo de PEG

# ---- formulaciones reales, para situarlas sobre el mapa (tarea 3.9) ----
# Marcarlas de forma DISTINTA a los tres diseños teóricos: estas sí son medidas.
FORMULACIONES_REALES = [
    dict(nombre="Mao 2014 (LP-FTY720)", R_nm=78.75, zeta_mV=+3.99, peg_nm=5.0,
         fuente="Nanomedicine 10:393, Tabla 1", via="intravenosa, leucemia (no BHE)"),
    dict(nombre="Gong 2022 (PSL-FTY720/AB)", R_nm=72.5, zeta_mV=-28.33, peg_nm=5.0,
         fuente="Nanophotonics 11:5133, Sec. 3.1", via="hemorragia intracerebral"),
    dict(nombre="Chow 2025 (nanosuspensión)", R_nm=67.0, zeta_mV=-0.24, peg_nm=0.0,
         fuente="Drug Deliv Transl Res 15:2022", via="nasal directa al cerebro"),
]
# NOTA (Mao 2014, Tabla 1 y Resultados): el liposoma VACÍO tiene zeta -4.10 +- 0.34 mV
# y al cargar el fingolimod sube a +3.99 +- 1.67 mV. El fármaco es catiónico a pH 7.4
# (amino protonado, pKa ~8.9), así que desplaza la carga de la partícula hacia positivo.
# Es la medida real que respalda el supuesto de zeta ligeramente positivo del estudio.

# ---- compuerta superior independiente ----
DIAM_CAVEOLA_MAX_nm = 80.0   # Bastiani & Parton 2010, "Caveolae at a glance".
                             # Límite geométrico superior, física distinta a la
                             # compuerta de envolvimiento (que es inferior).

# ---- constante de escalado de Deserno ----
A_DESERNO = 5.650   # Deserno 2004, Sec. V. Ajuste asintótico a alta tensión;
                    # el mismo valor sale de cuatro variables distintas con
                    # dispersión relativa de 6e-4.


# =============================================================================
#  BLOQUE 2 · DLVO: energía de interacción esfera-plano G(D)
# =============================================================================

def energia_libre_J(D_nm, R_nm, zeta_mV, peg_nm, hamaker_J,
                    zeta_bhe_mV=ZETA_BHE_mV):
    """G(D) en JULIOS. D = separación superficie-a-superficie en nm.

    Cuatro términos:
      1. van der Waals, esfera-plano: -A*R/(6D). ATRAE siempre. Escala con R.
      2. Doble capa eléctrica (superposición lineal). El signo sale solo:
         cargas opuestas -> atractivo; mismo signo -> repulsivo.
      3. Hidratación: repulsiva, corto alcance.
      4. Cepillo de PEG: repulsivo, solo si peg_nm > 0.

    A DIFERENCIA del cribado, aquí NO hay regularización D_sat: el corte
    D >= D0 = 0.165 nm ya impide la divergencia, y añadir D_sat destruiría
    la dependencia con R (ver cabecera del módulo).
    """
    D = np.maximum(np.asarray(D_nm, dtype=float), D0_nm)
    D_m = D * NM
    R_m = R_nm * NM

    # 1) van der Waals
    g_vdw = -hamaker_J * R_m / (6.0 * D_m)

    # 2) doble capa eléctrica
    kT_e = KT_J / QE
    g1 = np.tanh((zeta_mV * 1e-3) / (4.0 * kT_e))
    g2 = np.tanh((zeta_bhe_mV * 1e-3) / (4.0 * kT_e))
    pref = 64.0 * np.pi * EPS0 * EPS_R * R_m * (kT_e ** 2)
    g_edl = pref * g1 * g2 * np.exp(-D_m / (DEBYE_nm * NM))

    # 3) hidratación   4) PEG
    g_hid = B_HID_kT * KT_J * np.exp(-(D - D0_nm) / W_HID_nm)
    g_peg = (B_PEG_kT * KT_J * np.exp(-(D - D0_nm) / peg_nm)) if peg_nm > 0 else 0.0

    return g_vdw + g_edl + g_hid + g_peg


# =============================================================================
#  BLOQUE 3 · EL PUENTE: de G(D) esfera-plano a w (energía de adhesión por área)
#
#  ESTE ES EL PUNTO DELICADO DEL MÓDULO. Léelo antes de tocar nada.
#
#  El cronograma v3 proponía  w = G_min / (2*pi*R).  Esa fórmula está MAL:
#  G_min es energía [J] y 2*pi*R es longitud [m], así que J/m = N, una FUERZA.
#  Pero w es energía por unidad de área [J/m^2 = N/m]. No cuadran las unidades.
#
#  La aproximación de Derjaguin correcta para esfera-plano relaciona la FUERZA
#  con la energía por área de dos placas planas:
#
#      F_esfera-plano(D) = 2*pi*R * W_plano(D)
#
#  y como F = -dG/dD, invirtiendo:
#
#      W_plano(D) = -(1 / (2*pi*R)) * dG/dD
#
#  La energía de adhesión de Deserno es la profundidad de ese pozo plano-plano:
#
#      w = -min_D [ W_plano(D) ]        (positiva por convención)
#
#  COMPROBACIÓN ANALÍTICA (la hace test_limites()): para van der Waals puro,
#  G = -A*R/(6D), esta inversión debe dar W_plano = -A/(12*pi*D^2), que es
#  independiente de R. Coincide en cuatro decimales.
#
#  LIMITACIÓN QUE HAY QUE DECLARAR EN EL ARTÍCULO: Derjaguin supone superficies
#  RÍGIDAS e indeformables, mientras que el modelo de Deserno describe
#  precisamente la DEFORMACIÓN de la membrana. Encadenarlos es un híbrido: se
#  usa Derjaguin sobre la geometría NO deformada para estimar w, y esa w
#  constante se le entrega al modelo elástico. Es lo que se hace en la
#  literatura, pero no es gratis y no debe darse por supuesto.
# =============================================================================

def w_adhesion(R_nm, zeta_mV, peg_nm, hamaker_J, n_puntos=200000,
               devolver_perfil=False):
    """Energía de adhesión por unidad de área w, en J/m^2. Positiva = adhiere.

    Devuelve 0.0 si el potencial plano-plano no tiene pozo atractivo.
    """
    D = np.linspace(D0_nm, D_MAX_nm, n_puntos)
    D_m = D * NM
    G = energia_libre_J(D, R_nm, zeta_mV, peg_nm, hamaker_J)
    W_plano = -np.gradient(G, D_m) / (2.0 * np.pi * R_nm * NM)
    w = -float(W_plano.min())
    w = max(w, 0.0)
    if devolver_perfil:
        return w, D, G, W_plano
    return w


def pozo_y_barrera_kT(R_nm, zeta_mV, peg_nm, hamaker_J, n_puntos=40000):
    """Descriptores de G(D) en kT. Sustituyen a la fracción adherida, que
    satura al 100% y no discrimina (tarea 3.2 del cronograma).

    Devuelve un dict con:
      pozo_kT        profundidad del mínimo (NEGATIVO si atrae)
      D_pozo_nm      posición del mínimo
      barrera_entrada_kT
                     altura del máximo local de G que hay que cruzar VINIENDO
                     DESDE LEJOS para llegar al pozo, medida respecto al campo
                     lejano (G -> 0). Es 0 si G decrece monótonamente hacia el
                     contacto, es decir, si NO hay barrera. Esta es la barrera
                     DLVO clásica.
      escape_kT      profundidad desde el pozo hasta el campo lejano, = |pozo|.
                     Es lo que costaría DESPEGARSE. No confundir con la anterior.

    OJO: son magnitudes distintas. Una versión anterior de este módulo las
    confundía y reportaba `escape` llamándolo barrera.
    """
    D = np.linspace(D0_nm, D_MAX_nm, n_puntos)
    G = energia_libre_J(D, R_nm, zeta_mV, peg_nm, hamaker_J) / KT_J
    i_min = int(np.argmin(G))
    pozo = float(G[i_min])
    G_lejos = float(G[-1])

    # barrera de entrada: máximo de G por FUERA del pozo, respecto al campo lejano
    if i_min < len(G) - 1:
        pico = float(G[i_min:].max())
        barrera_entrada = max(0.0, pico - G_lejos)
    else:
        barrera_entrada = 0.0

    return dict(pozo_kT=pozo,
                D_pozo_nm=float(D[i_min]),
                barrera_entrada_kT=barrera_entrada,
                escape_kT=float(G_lejos - pozo))


# =============================================================================
#  BLOQUE 4: CRITERIOS DE ENVOLVIMIENTO (Deserno 2004)
# =============================================================================

def adimensionales(R_nm, w, kappa_kT, sigma_mNm):
    """(w_tilde, sigma_tilde) de Deserno, Ec. (3).

        w_tilde     = 2*w*a^2 / kappa
        sigma_tilde = sigma*a^2 / kappa       (= (a/lambda)^2)
    """
    a = R_nm * NM
    kappa = kappa_kT * KT_J
    sigma = sigma_mNm * 1e-3            # mN/m -> N/m
    return 2.0 * w * a * a / kappa, sigma * a * a / kappa


def radio_critico_nm(w, kappa_kT):
    """R_min = sqrt(2*kappa/w). Deserno 2004, Sec. III C (textual).

    Por debajo de este radio la membrana NO envuelve, por mucha afinidad que
    haya: la adhesión no puede pagar el coste de curvatura.

    OJO: es el umbral del INICIO del envolvimiento (parcial), no del
    envolvimiento COMPLETO. Para eso, ver frontera_envolvimiento_completo().
    """
    if w <= 0:
        return np.inf
    kappa = kappa_kT * KT_J
    return float(np.sqrt(2.0 * kappa / w) / NM)


def _razon_frontera(sigma_tilde):
    """Cociente (w_tilde - 4)/sigma_tilde en la frontera de envolvimiento completo.

    Deserno 2004 da dos asíntotas y demuestra que el cruce entre ellas es
    LENTÍSIMO (Fig. 6):
      · tensión baja  (sigma_tilde <~ 1):   cociente -> 2     [w_tilde = 4 + 2*sigma_tilde]
      · tensión alta  (sigma_tilde >> 1):   cociente -> 4     [límite de Young-Dupré]
        con aproximación a la asíntota, Ec. (35):
            (w_tilde-4)/sigma_tilde ~= 4 - 3*A^(2/3)*sigma_tilde^(-1/3)

    En medio no vale ninguna fórmula cerrada. Aquí se INTERPOLA de forma monótona
    en log10(sigma_tilde) entre 2 y la rama alta. Es una interpolación que
    reproduce la Fig. 6 cualitativamente, NO una derivación. Declarar como tal.
    """
    st = np.maximum(np.asarray(sigma_tilde, dtype=float), 1e-12)
    alta = 4.0 - 3.0 * (A_DESERNO ** (2.0 / 3.0)) * st ** (-1.0 / 3.0)
    alta = np.clip(alta, 2.0, 4.0)
    # peso logístico en log10(sigma_tilde), centrado en el cruce sigma_tilde ~ 10
    peso = 1.0 / (1.0 + 10.0 ** (-(np.log10(st) - 1.0)))
    return (1.0 - peso) * 2.0 + peso * alta


def frontera_envolvimiento_completo(sigma_tilde):
    """w_tilde mínimo para que el envolvimiento sea COMPLETO."""
    return 4.0 + _razon_frontera(sigma_tilde) * np.asarray(sigma_tilde, dtype=float)


def barrera_envolvimiento_kT(sigma_tilde, kappa_kT):
    """Altura de la barrera energética de la transición de envolvimiento, en kT.

    La transición de parcialmente envuelto a completamente envuelto es
    DISCONTINUA y tiene una barrera (tensión almacenada en la membrana libre).
    Deserno 2004, Sec. III C y Fig. 9:
      · alta tensión:  E_tilde ~ (3/4)*(2*sqrt(3)-3)*A^(4/3)*sigma_tilde^(1/3)
      · baja tensión:  ley de potencias empírica con exponente 0.86
      · las dos ramas se cruzan en sigma_tilde ~= 4.72
      · la barrera SE ANULA en el límite sigma_tilde -> 0
    con E_barrier = E_tilde * pi * kappa.

    !! PRECISIÓN LIMITADA. Esto es una interpolación de escalados, no una
    fórmula exacta: reproduce el ejemplo numérico de Deserno (22 kT a
    sigma_tilde=0.22 con kappa=20 kT) por construcción, pero se desvía hasta un
    ~50% del segundo ejemplo (66 kT a sigma_tilde=1). USAR COMO ORDEN DE
    MAGNITUD, y decirlo así en el artículo.

    Es una compuerta CINÉTICA, no termodinámica: un diseño puede ser
    termodinámicamente favorable y estar bloqueado por decenas de kT.
    """
    st = np.maximum(np.asarray(sigma_tilde, dtype=float), 0.0)
    st_cruce = 4.72
    # rama alta
    coef_alta = 0.75 * (2.0 * np.sqrt(3.0) - 3.0) * A_DESERNO ** (4.0 / 3.0)
    # rama baja, prefactor calibrado al ejemplo numérico de Deserno
    #   sigma_tilde = 0.22, kappa = 20 kT  ->  E_barrier = 22 kT  ->  E_tilde = 0.3501
    coef_baja = 0.3501 / (0.22 ** 0.86)
    with np.errstate(divide="ignore", invalid="ignore"):
        e_baja = coef_baja * np.where(st > 0, st, 0.0) ** 0.86
        e_alta = coef_alta * np.where(st > 0, st, 0.0) ** (1.0 / 3.0)
    e_tilde = np.where(st < st_cruce, e_baja, e_alta)
    return np.asarray(e_tilde * np.pi * kappa_kT, dtype=float)


def clasificar(R_nm, zeta_mV, peg_nm, hamaker_J, kappa_kT, sigma_mNm):
    """Evalúa TODAS las compuertas para un diseño. Devuelve un dict.

    Compuertas:
      G1  radio critico     R > R_min = sqrt(2*kappa/w)        (inferior, energética)
      G2  cociente w/sigma  w/sigma > 1.37                     (independiente del tamaño)
      G3  envolvimiento completo   w_tilde > 4 + razon*sigma_tilde
      G4  caveola           diametro <= 80 nm                  (superior, geométrica)
    """
    w = w_adhesion(R_nm, zeta_mV, peg_nm, hamaker_J)
    wt, st = adimensionales(R_nm, w, kappa_kT, sigma_mNm)
    r_min = radio_critico_nm(w, kappa_kT)
    sigma_SI = sigma_mNm * 1e-3
    perfil = pozo_y_barrera_kT(R_nm, zeta_mV, peg_nm, hamaker_J)

    g1 = R_nm > r_min
    g2 = (w / sigma_SI) > 1.37 if sigma_SI > 0 else True
    g3 = wt > float(frontera_envolvimiento_completo(st))
    g4 = (2.0 * R_nm) <= DIAM_CAVEOLA_MAX_nm

    return dict(
        R_nm=R_nm, zeta_mV=zeta_mV, peg_nm=peg_nm,
        hamaker_J=hamaker_J, kappa_kT=kappa_kT, sigma_mNm=sigma_mNm,
        w_uNm=w * 1e6, w_tilde=wt, sigma_tilde=st,
        R_min_nm=r_min, w_sobre_sigma=(w / sigma_SI if sigma_SI > 0 else np.inf),
        pozo_dlvo_kT=perfil["pozo_kT"],
        D_pozo_nm=perfil["D_pozo_nm"],
        barrera_entrada_kT=perfil["barrera_entrada_kT"],
        escape_kT=perfil["escape_kT"],
        barrera_envolv_kT=float(barrera_envolvimiento_kT(st, kappa_kT)),
        G1_radio_critico=bool(g1), G2_w_sobre_sigma=bool(g2),
        G3_envolv_completo=bool(g3), G4_caveola=bool(g4),
        NO_EXCLUIDO=bool(g1 and g2 and g3 and g4),
    )


# =============================================================================
#  BLOQUE 5: VALIDACIÓN (tarea 2.6) Y CONTROL DE FALSABILIDAD
#  Un modelo que solo sabe decir "sí" no prueba nada. Este sabe decir "no".
# =============================================================================

def test_limites(verbose=True):
    """Casos límite con respuesta conocida. Devuelve True si pasan todos."""
    ok = []

    def chequeo(nombre, condicion, detalle=""):
        ok.append(bool(condicion))
        if verbose:
            print(f"  [{'OK ' if condicion else 'FALLA'}] {nombre}{'  ' + detalle if detalle else ''}")

    if verbose:
        print("=" * 78)
        print(" VALIDACIÓN DEL MÓDULO (tarea 2.6)")
        print("=" * 78)

    # T1: Derjaguin contra el resultado analítico del vdW puro.
    A = 4.5e-21
    for R in (20.0, 50.0, 100.0):
        D = np.linspace(D0_nm, 3.0, 200000)
        D_m = D * NM
        G = -A * (R * NM) / (6.0 * D_m)
        W_num = -np.gradient(G, D_m) / (2.0 * np.pi * R * NM)
        W_teo = -A / (12.0 * np.pi * D_m ** 2)
        err = abs(W_num[0] - W_teo[0]) / abs(W_teo[0])
        chequeo(f"T1 Derjaguin vs analítico (R={R:.0f} nm)", err < 1e-3,
                f"error relativo {err:.2e}")

    # T2: w independiente de R para vdW puro (lo exige Derjaguin).
    ws = [w_adhesion(R, 0.0, 0.0, A) for R in (20.0, 50.0, 100.0)]
    disp = (max(ws) - min(ws)) / np.mean(ws)
    chequeo("T2 w casi independiente de R (vdW dominante)", disp < 0.05,
            f"dispersión {disp:.2%}")

    # T3: adhesión nula => nunca envuelve.
    chequeo("T3 w=0 -> R_min infinito", np.isinf(radio_critico_nm(0.0, 25.0)))

    # T4: membrana infinitamente blanda => siempre envuelve.
    r_blanda = radio_critico_nm(4.4e-3, 1e-6)
    chequeo("T4 kappa->0 -> R_min->0", r_blanda < 1e-2,
            f"R_min = {r_blanda:.2e} nm")

    # T5: sigma=0 => barrera nula (Deserno: catenoide, transición sin barrera).
    chequeo("T5 sigma_tilde=0 -> barrera de envolvimiento = 0",
            abs(float(barrera_envolvimiento_kT(0.0, 25.0))) < 1e-9)

    # T6: reproduce el ejemplo numérico de Deserno (22 kT).
    b = float(barrera_envolvimiento_kT(0.22, 20.0))
    chequeo("T6 barrera a sigma_tilde=0.22, kappa=20 kT", abs(b - 22.0) < 1.0,
            f"da {b:.1f} kT (Deserno: ~22 kT)")

    # T7 · FALSABILIDAD: aniónico fuerte debe adherir mucho peor que catiónico.
    w_cat = w_adhesion(20.0, +20.0, 0.0, A)
    w_ani = w_adhesion(20.0, -30.0, 0.0, A)
    chequeo("T7 falsabilidad: w(aniónico) < w(catiónico)", w_ani < w_cat,
            f"{w_ani*1e6:.0f} vs {w_cat*1e6:.0f} uN/m")

    # T8 · FALSABILIDAD: el PEG debe reducir la adhesión.
    w_sin = w_adhesion(20.0, +5.0, 0.0, A)
    w_con = w_adhesion(20.0, +5.0, 10.0, A)
    chequeo("T8 falsabilidad: el PEG reduce w", w_con < w_sin,
            f"{w_con*1e6:.0f} vs {w_sin*1e6:.0f} uN/m")

    # T9 · monotonía: más kappa => radio crítico mayor.
    r15 = radio_critico_nm(4.4e-3, 15.0)
    r50 = radio_critico_nm(4.4e-3, 50.0)
    chequeo("T9 monotonía: kappa mayor -> R_min mayor", r50 > r15,
            f"{r15:.1f} -> {r50:.1f} nm")

    # T10: sin NaN ni infinitos en un barrido representativo.
    vals = [w_adhesion(R, z, p, h)
            for R in (10.0, 40.0, 75.0)
            for z in (-30.0, 0.0, 30.0)
            for p in (0.0, 5.0)
            for h in HAMAKER_J]
    chequeo("T10 sin NaN ni infinitos en el barrido", np.all(np.isfinite(vals)))

    # T11: el barrido de kappa tiene que cubrir el caso saturado + colesterol.
    #       Pan 2008: K_C de DMPC sube MÁS DE 4 veces con 30 % de colesterol.
    #       Con el K_Ctd de DMPC de Nagle 2017 (24.6 kT) eso pasa de 98 kT, así
    #       que el techo del barrido no puede quedarse en 50 como estuvo hasta
    #       el 2026-08-13.
    chequeo("T11 el barrido de kappa cubre saturado + colesterol (>98 kT)",
            max(KAPPA_kT) >= 4.0 * 24.6,
            f"techo {max(KAPPA_kT):.0f} kT frente a {4.0*24.6:.1f} kT")

    # T12: MARGEN REAL del diseño más pequeño frente a la rigidez. Es el número
    #       que el barrido no decía: R_min = sqrt(2k/w) = R  ->  k = w·R²/2.
    #       Con R = 15.5 nm (el furtivo, el más pequeño del catálogo) el modelo
    #       aguanta hasta ~116 kT. TRIPWIRE: si entra un diseño más pequeño o
    #       cambia w, esta prueba debe fallar para obligar a recalcularlo.
    _R_min_catalogo_nm = 15.5
    _w = w_adhesion(_R_min_catalogo_nm, 2.0, 5.0, HAMAKER_J[1])
    _k_critico = _w * (_R_min_catalogo_nm * NM) ** 2 / 2.0 / KT_J
    chequeo("T12 el diseño más pequeño aguanta el techo del barrido",
            _k_critico > max(KAPPA_kT),
            f"rompe a {_k_critico:.1f} kT, techo del barrido {max(KAPPA_kT):.0f} kT "
            f"(margen {100*(_k_critico/max(KAPPA_kT)-1):.0f} %)")

    if verbose:
        print("-" * 78)
        print(f" RESULTADO: {sum(ok)}/{len(ok)} pruebas superadas")
        print("=" * 78)
    return all(ok)


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
    print(f" Endotelio: zeta = {ZETA_BHE_mV} mV (Santa-Maria 2019, Fig. 4A). T = 37 °C.")
    print(f" Barridos de IGNORANCIA  kappa {KAPPA_kT} kT | sigma {SIGMA_mNm} mN/m"
          f" | Hamaker {[f'{h:.1e}' for h in HAMAKER_J]} J")
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
        r = clasificar(lip["R_nm"], lip["zeta_mV"], lip["peg_nm"],
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
    for s in SIGMA_mNm:
        for k in KAPPA_kT:
            r = clasificar(lip["R_nm"], lip["zeta_mV"], lip["peg_nm"],
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
    D = np.linspace(D0_nm, 6.0, 4000)
    for lip in LIPOSOMAS:
        G = energia_libre_J(D, lip["R_nm"], lip["zeta_mV"],
                              lip["peg_nm"], HAMAKER_REF) / KT_J
        ax.plot(D, G, label=lip["nombre"])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("separación D (nm)"); ax.set_ylabel("G(D)  [kT]")
    ax.set_title("Energía libre de interacción DLVO")
    ax.set_ylim(-30, 15); ax.legend(); fig.tight_layout()
    fig.savefig(f"{prefijo}_G_de_D.png", dpi=160); plt.close(fig)

    # Fig 2: radio crítico frente a w, para los tres kappa
    fig, ax = plt.subplots(figsize=(7, 4.5))
    w_grid = np.logspace(-5, -1, 400)
    for k in KAPPA_kT:
        ax.plot(w_grid * 1e6, [radio_critico_nm(w, k) for w in w_grid],
                label=f"kappa = {k:.0f} kT")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("w  (uN/m)"); ax.set_ylabel("R_min  (nm)")
    ax.set_title("Radio crítico de envolvimiento  R_min = sqrt(2k/w)")
    ax.legend(); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(f"{prefijo}_radio_critico.png", dpi=160); plt.close(fig)

    # Fig 3: barrera de envolvimiento frente a sigma_tilde
    fig, ax = plt.subplots(figsize=(7, 4.5))
    st = np.logspace(-3, 3, 600)
    for k in KAPPA_kT:
        ax.plot(st, barrera_envolvimiento_kT(st, k), label=f"kappa = {k:.0f} kT")
    ax.axhline(3.0, color="gray", ls="--", lw=0.8)
    ax.text(1e-3, 3.4, "~3 kT: cruzable térmicamente", fontsize=8, color="gray")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("tensión reducida  sigma_tilde"); ax.set_ylabel("barrera  (kT)")
    ax.set_title("Barrera de la transición de envolvimiento (orden de magnitud)")
    ax.legend(); ax.grid(alpha=0.3, which="both"); fig.tight_layout()
    fig.savefig(f"{prefijo}_barrera.png", dpi=160); plt.close(fig)

    print(f"\n Figuras guardadas: {prefijo}_G_de_D.png, "
          f"{prefijo}_radio_critico.png, {prefijo}_barrera.png")



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
