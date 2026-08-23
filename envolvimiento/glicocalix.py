#!/usr/bin/env python3
# =============================================================================
#  GLICOCÁLIX: tarea 3.1 del cronograma v3
#  Proyecto BHE / SENACYT
#
#  VERSIÓN 2 (2026-08-08): reescrita tras conseguir la fuente PRIMARIA.
#  La versión 1 usaba tamaños de malla tomados de resúmenes de búsqueda, sin
#  verificar. Al leer el artículo original resultó que uno de esos números
#  estaba mal interpretado. Ver "CORRECCIÓN" más abajo.
# =============================================================================
#
#  FUENTE PRIMARIA (leída íntegra, acceso libre)
#  ---------------------------------------------
#  Weinbaum S, Zhang X, Han Y, Vink H, Cowin SC.
#  "Mechanotransduction and flow across the endothelial glycocalyx".
#  PNAS 100(13):7988-7995 (2003).  doi:10.1073/pnas.1332808100
#  Sección "Transport Model". Datos estructurales de:
#  Squire JM, Chew M, Nneji G, Neal C, Barry J, Michel C.
#  J. Struct. Biol. 136:239-255 (2001)  [NO leído: citado a través de Weinbaum].
#
#  CORRECCIÓN respecto a la versión 1 de este módulo
#  --------------------------------------------------
#  La v1 tomaba 20 nm como el TAMAÑO DE HUECO de la malla. Es incorrecto.
#  Weinbaum, Sección "Transport Model", es explícito:
#     · 20 nm  = espaciado CENTRO A CENTRO de los centros de dispersión
#     · 10-12 nm = diámetro de esos centros (radio de fibra rf = 6 nm)
#     · 8 nm   = HUECO real entre fibras (Delta), que es 20 - 2*6
#  El hueco por el que puede pasar algo es 8 nm, no 20.
#
#  CONSECUENCIA: la "coincidencia" que la v1 señalaba (diámetro máximo de
#  53.5 nm, cerca del '< 50 nm' de la hipótesis) era un ARTEFACTO de ese error.
#  Queda RETIRADA. Con los valores correctos el límite es de ~9 nm.
#
#  MODELO PRINCIPAL: teoría de matriz de fibras
#  --------------------------------------------
#  Weinbaum, Sección "Transport Model", ecuaciones citadas de su ref. 28:
#
#      Vf   = 2 pi rf^2 / ( sqrt(3) (2 rf + Delta)^2 )      (red hexagonal)
#      phi  = 1 - Vf (1 + a/rf)^2                           (coef. de partición)
#      sigma = (1 - phi)^2                                  (coef. de reflexión)
#
#  con a = radio del soluto, rf = radio de fibra, Delta = hueco entre fibras.
#  phi es la fracción del espacio de poro accesible al soluto: phi = 0 significa
#  EXCLUSIÓN ESTÉRICA TOTAL.
#
#  VERIFICACIÓN: estas ecuaciones reproducen los dos valores que el propio
#  artículo publica, a tres decimales (ver test_glicocalix, T1 y T2):
#      · modelo nuevo (rf=6 nm, Delta=8 nm), albúmina a=3.5 nm  -> sigma = 0.670
#        (el artículo dice 0.67)
#      · modelo antiguo (rf=0.6 nm, Delta=8 nm), misma albúmina -> sigma = 0.519
#        (el artículo dice 0.52)
#
#  POR QUÉ NO SE USA ALEXANDER–DE GENNES, COMO DECÍA EL CRONOGRAMA
#  ---------------------------------------------------------------
#  AdG describe la COMPRESIÓN de dos superficies planas recubiertas de cepillo,
#  y exige que la partícula sea grande frente al cepillo. Aquí el grosor del
#  glicocálix es 150-400 nm y el radio de los liposomas 15-20 nm: R/L ~ 0.04-0.13.
#  La partícula no comprime nada, tiene que atravesar una red mayor que ella.
#  AdG se conserva SOLO en regimen_adg(), como comprobación de que no aplica.
# =============================================================================

import numpy as np

# =============================================================================
#  PARÁMETROS: todos con fuente
# =============================================================================

# ---- geometría de la matriz de fibras (Weinbaum 2003, "Transport Model") ----
RF_nm = 6.0        # radio de fibra efectivo (centros de dispersión de 10-12 nm
                   # de diámetro). Weinbaum lo usa como valor de trabajo.
DELTA_nm = 8.0     # hueco entre fibras. Weinbaum: "rf = 6 nm, Delta = 8 nm".
ESPACIADO_nm = 20.0  # centro a centro = 2*rf + Delta. Solo informativo.

# Modelo ANTIGUO de cadenas de GAG extendidas, para comparar (mismo artículo):
RF_GAG_nm = 0.6
DELTA_GAG_nm = 8.0

# ---- grosor del glicocálix ----
# PRIMARIO (Weinbaum 2003): 150 nm en capilar de mesenterio de rana (su ref. 15)
#                           400 nm en cremáster de hámster (su ref. 2)
# Weinbaum discute explícitamente esa discrepancia y le da explicación funcional.
L_GLICO_nm = [150.0, 400.0]
# SEGUNDA MANO (Walter 2021, pág. e1904773-3, su ref. 20): 0.2-5 um para
# endotelio en general. Rango mucho más ancho. Se conserva solo como contraste.
L_GLICO_WALTER_nm = [200.0, 5000.0]

# ---- grosor en ENDOTELIO CEREBRAL. Añadido: 2026-08-13 ----------------
# Hasta hoy este archivo decía "NO hay medida para endotelio cerebral" y usaba
# rana y hámster. Es FALSO desde 2025: hay dos medidas directas en capilar
# cerebral de ratón, las dos por microscopía electrónica con tinción de nitrato
# de lantano y sobre animal perfundido (o sea, con flujo, no cultivo).
#
#   Shi et al. 2025, Nature 639:985, Fig. 1c   -> 540 ± 86 nm  (ratón 3 meses)
#                                                 232 ± 92 nm  (ratón 21 meses)
#   Larsen et al. 2025, bioRxiv 2025.04.07.645297 -> 726 ± 148 nm total
#                                                    412 ±  98 nm capa densa
#
# Los dos SUPERAN el rango de rana/hámster que usa el proyecto: el glicocálix
# cerebral es el más grueso que se ha medido, y Larsen lo compara de frente con
# músculo, corazón, riñón e hígado.
#
# NO gobiernan todavía ninguna compuerta. Tres motivos para la prudencia:
#   1. Larsen es PREPRINT, sin revisión por pares (comprobado el 2026-08-13:
#      sigue sin publicar y con 0 citas en Europe PMC).
#   2. Shi (540) y Larsen (726) no coinciden pese a usar la misma técnica en el
#      mismo animal; puede que no midan lo mismo con el mismo nombre.
#   3. Ninguno de los dos da radio de fibra ni espaciado, que es de lo que
#      depende el tamiz de 9 nm. El tamiz NO cambia con esto.
L_GLICO_CEREBRAL_nm = [540.0, 726.0]     # Shi 2025 joven · Larsen 2025 total
L_GLICO_CEREBRAL_ENVEJECIDO_nm = 232.0   # Shi 2025, ratón de 21 meses
COBERTURA_GLICO_CEREBRAL = 0.93          # Larsen 2025, 93 ± 2.5 %

# ---- barrido de incertidumbre sobre el hueco ----
# NO hay medida para endotelio cerebral. Walter 2021 señala que el glicocálix
# de microvasos cerebrales es MÁS DENSO que el de pulmón o hígado (su ref. 23),
# lo que empujaría hacia huecos MENORES, pero no da número.
DELTA_BARRIDO_nm = [4.0, 8.0, 16.0]


# =============================================================================
#  MODELO PRINCIPAL: teoría de matriz de fibras
# =============================================================================

def fraccion_volumen_fibra(rf_nm=RF_nm, delta_nm=DELTA_nm):
    """Vf = 2 pi rf^2 / (sqrt(3) (2 rf + Delta)^2). Red hexagonal."""
    return 2.0 * np.pi * rf_nm ** 2 / (np.sqrt(3.0) * (2.0 * rf_nm + delta_nm) ** 2)


def coef_particion(a_nm, rf_nm=RF_nm, delta_nm=DELTA_nm):
    """phi = 1 - Vf (1 + a/rf)^2, acotado a [0, 1].

    phi es la fracción del espacio de poro accesible a un soluto de radio a.
    phi = 0  ->  exclusión estérica total.
    """
    Vf = fraccion_volumen_fibra(rf_nm, delta_nm)
    phi = 1.0 - Vf * (1.0 + np.asarray(a_nm, dtype=float) / rf_nm) ** 2
    return np.clip(phi, 0.0, 1.0)


def coef_reflexion(a_nm, rf_nm=RF_nm, delta_nm=DELTA_nm):
    """sigma = (1 - phi)^2. sigma = 1 -> el soluto se refleja por completo."""
    return (1.0 - coef_particion(a_nm, rf_nm, delta_nm)) ** 2


def radio_exclusion_nm(rf_nm=RF_nm, delta_nm=DELTA_nm):
    """Radio del soluto al que phi llega a 0: exclusión estérica total.

    De phi = 0:   Vf (1 + a/rf)^2 = 1   ->   a = rf (1/sqrt(Vf) - 1)
    """
    Vf = fraccion_volumen_fibra(rf_nm, delta_nm)
    return rf_nm * (1.0 / np.sqrt(Vf) - 1.0)


def compuerta_glicocalix(R_nm, rf_nm=RF_nm, delta_nm=DELTA_nm):
    """Compuerta G5: ¿es la partícula estéricamente admisible en la matriz?"""
    phi = float(coef_particion(R_nm, rf_nm, delta_nm))
    a_max = radio_exclusion_nm(rf_nm, delta_nm)
    return dict(R_nm=R_nm, rf_nm=rf_nm, delta_nm=delta_nm,
                phi=phi, sigma=float(coef_reflexion(R_nm, rf_nm, delta_nm)),
                a_max_nm=a_max, d_max_nm=2.0 * a_max,
                G5_admisible=bool(phi > 0.0))


# =============================================================================
#  PMF DEL GLICOCÁLIX: Kabedev & Lobaskin 2022 (datos DIGITALIZADOS)
#  Añadido: 2026-08-18. Punto 4 del roadmap de viabilidad predictiva (18b/18c).
# =============================================================================
#
#  FUENTE (leída íntegra, open access + Supplementary Material Text 1, 5 pág.)
#  ---------------------------------------------------------------------------
#  Kabedev A, Lobaskin V (2022). "Endothelial glycocalyx permeability for
#  nanoscale solutes". Nanomedicine (Lond) 17(13):979-996.
#  DOI: 10.2217/nnm-2021-0367. PMID: 35815713.
#
#  El artículo da el PMF(z) del glicocálix SOLO EN FIGURAS (Fig. 4B y Fig. 6),
#  sin fórmula cerrada ni tabla numérica, confirmado tras leer el Supplementary
#  Material completo (solo trae teoría de vdW descartada + detalle del modelo
#  hidrodinámico + heatmaps de distribución en el plano xy; nada de PMF
#  tabulado). Los números de aquí son DIGITALIZACIÓN de esas gráficas, medida a
#  pixel sobre el PDF a 600 dpi con calibración de ejes contra las marcas de
#  los ticks (no a ojo). NO son datos tabulados por los autores, declarado
#  así por la regla dura del proyecto de no inventar números.
#
#  QUÉ MODELO ES: Core Protein Model (CPM), sigma = 0.0019 nm^-2, el mismo
#  modelo de matriz de fibras que usa compuerta_glicocalix() arriba (aunque la
#  parametrización geométrica de Kabedev NO es la misma que la de Weinbaum
#  2003; son dos modelos independientes del mismo tipo de estructura).
#
#  DOS DEFINICIONES DE "BARRERA" QUE NO SON INTERCAMBIABLES
#  ---------------------------------------------------------------------------
#  1. PICO DE CONTACTO (z ~ 0-10 nm): el valor extremo justo al tocar la
#     primera fibra. Es lo único que Fig. 4B muestra con claridad (un pico
#     agudo por R, sin solape entre curvas).
#  2. HOMBRO/MESETA (z ~ 20-450 nm): el valor sostenido mientras la partícula
#     recorre el grosor del glicocálix. Es lo único que Fig. 6 permite separar
#     con confianza entre neutra/negativa/positiva (el pico de contacto queda
#     fuera del rango de esos ejes para casi todas las curvas).
#  Decisión de Jhovan (2026-08-18): usar HOMBRO/MESETA como definición de
#  barrera para G5, por ser la resistencia sostenida a lo largo del grosor de
#  la capa (coherente con "atravesar la malla"), no un pico puntual de
#  contacto, y por ser la única definición que Fig. 6 deja separar por carga.
#
#  PMF_HOMBRO_NEUTRO_kT: Fig. 4B, meseta tras el pico de contacto (NO el pico)
#  ---------------------------------------------------------------------------
#  Para R = 3.5, 5, 7.5, 10 nm el pico de contacto SÍ se resolvió con
#  confianza alta (línea aislada, sin solape). Para R = 15 y 20 nm, la propia
#  figura avisa: "For all the curves, only the last noninfinite values are
#  shown (for the chosen sampling bin)", el muestreo se cortó ANTES de
#  alcanzar el pico real, así que el valor mostrado es un PISO (cota inferior),
#  no el pico verdadero. Se conserva la distinción explícita en el dict de
#  abajo con "tipo": "pico" (dato directo) o "piso" (cota inferior, dato real
#  desconocido y probablemente mayor).
PMF_CONTACTO_NEUTRO_kT = {
    # R_nm: (valor_kT, z_del_punto_nm, tipo)
    3.5:  (0.70,  5.0,   "pico"),
    5.0:  (1.76,  4.6,   "pico"),
    7.5:  (8.36,  5.0,   "pico"),
    10.0: (12.5,  5.5,   "pico"),
    15.0: (6.93,  43.0,  "piso"),   # muestreo cortado antes del pico real
    20.0: (9.36,  325.0, "piso"),   # muestreo cortado antes del pico real
}

#  PMF_HOMBRO_kT_A / _B: Fig. 6, región de meseta (z ~ 20-450 nm), por carga.
#  ---------------------------------------------------------------------------
#  _A = panel A, ρq-EG = 8.3 mEq/l. _B = panel B, ρq-EG = 25 mEq/l (misma
#  figura, digitalizado en sesión separada del mismo día). El texto del
#  artículo dice que el efecto de triplicar la carga del glicocálix es
#  "relativamente pequeño" para partículas negativas, con _A y _B ya no hace
#  falta creerle esa frase, se puede comparar directamente.
#
#  "negativo" (panel A) viene de Fig. 6A (línea punteada larga); R=10nm
#  negativo viene del INSET de Fig. 6A (el eje principal lo corta fuera de
#  escala), su valor es TAMBIÉN un piso: el inset empieza a mostrar la curva
#  en z=111 nm con 7.0 kT, pero para z<111 nm el valor sigue fuera del rango
#  del inset (0-8 kT), el pico real es mayor y desconocido.
#  "positivo" (panel A, R=10nm) viene de Fig. 6A (línea de cruces '+'); no es
#  un valor único sino DECRECIENTE con z desde ~2.2 kT en z=50nm, se guarda
#  el primer punto resuelto, no un pico ni una meseta plana.
#  R=5nm bajo carga en panel A: NO resuelto en la primera pasada (línea negra
#  confundida con bordes de otras curvas); SÍ resuelto después con muestreo
#  denso (30 rebanadas en z, promediado) para "neutro" y "negativo", ver _B,
#  mismo método aplicado retroactivamente a A donde fue posible.
#
#  PANEL B · "positivo" de R=10nm y "negativo" de R=10nm NO resueltos: la
#  curva "+" en este panel forma un POZO NO MONÓTONO (baja, toca un mínimo
#  hacia z~400, vuelve a subir) que se CRUZA en pantalla con la línea
#  punteada larga (negativa) alrededor de z=250-350, separarlas requiere
#  seguir cada curva por continuidad a través de z, no solo agrupar por
#  altura en una rebanada vertical (lo que sí funcionó para las demás curvas,
#  que no se cruzan entre sí en esa región). Forzar un número aquí sería
#  mezclar dos curvas físicamente distintas. Se deja "no resuelto" a
#  propósito. El resto de R=10nm negativo/positivo de panel B queda igual:
#  sin digitalizar por la misma razón.
PMF_HOMBRO_kT_A = {
    # R_nm: {"neutro": ..., "negativo": ..., "positivo": ...} en kT, o None si
    # no se resolvió con confianza.
    3.5: {"neutro": (0.47, "meseta"),
          "negativo": (0.95, "meseta"),
          "positivo": (0.03, "meseta, casi plana")},
    5.0: {"neutro": None, "negativo": None, "positivo": None},  # no resuelto
    10.0: {"neutro": (2.5, "meseta, consistente con Fig.4B post-pico"),
           "negativo": (7.0, "piso en z=111nm (inset); pico real mayor y desconocido"),
           "positivo": (2.2, "en z=50nm, DECRECIENTE con z, no es meseta plana")},
}

PMF_HOMBRO_kT_B = {
    3.5: {"neutro": (0.46, "meseta, promedio z=40-400nm, n=37 rebanadas"),
          "negativo": (1.03, "meseta, promedio z=40-400nm, n=27 rebanadas"),
          "positivo": None},  # pozo no monótono, cruza con la curva negativa
    5.0: {"neutro": (0.81, "meseta, promedio z=40-400nm, n=36 rebanadas"),
          "negativo": (1.49, "meseta, promedio z=40-400nm, n=13 rebanadas; menos denso, confianza media"),
          "positivo": None},  # ídem: pozo, llega a ~-5.7kT hacia z=420 (atractivo, no barrera)
    10.0: {"neutro": (2.48, "meseta, promedio z=40-400nm, n=37 rebanadas; consistente con panel A"),
           "negativo": None,  # sube a ~10.5kT hacia z=280-300 (lectura visual, no separada por
                               # pixel de la curva '+' que cruza ahí, no se declara como número)
           "positivo": None},
}

#  QUÉ NO HACE ESTE MÓDULO TODAVÍA
#  ---------------------------------------------------------------------------
#  No hay función kT(R, zeta) continua. Con huecos reales en R=15 y 20 nm
#  neutro (piso, no pico real) y en R=10nm cargado de panel B (curvas '+' y
#  punteada larga se cruzan, no separadas), interpolar ahora sería rellenar
#  esos huecos con una suposición, contra la regla dura del proyecto de no
#  inventar números. compuerta_glicocalix() (el booleano G5_admisible de
#  arriba) sigue siendo la única compuerta activa; estos diccionarios son
#  datos de referencia para cuando haya más puntos o un modelo teórico que
#  los conecte, no un reemplazo todavía.
#
#  PENDIENTE, en orden: (a) separar por trazado de curva (no por rebanada) el
#  negativo/positivo de R=10nm en panel B, que se cruzan alrededor de
#  z=250-350; (b) decidir con Jhovan si vale la pena un modelo teórico (p.ej.
#  Debye-Hückel de la propia Fig.6, Ec. 3-5 del artículo) que conecte los
#  puntos en vez de seguir digitalizando figura por figura, con _A y _B ya
#  cubiertos, el rendimiento marginal de seguir leyendo píxeles es bajo.


# =============================================================================
#  kT(R, carga): AJUSTE continuo, NO dato digitalizado
#  Añadido: 2026-08-18. Cierra los huecos de PMF_CONTACTO_NEUTRO_kT /
#  PMF_HOMBRO_kT_A / PMF_HOMBRO_kT_B con un modelo, no con más lectura de
#  píxeles. Decisión de Jhovan (2026-08-18): con solo 3 anclajes limpios y
#  huecos reales en R=15/20 y en la carga, el rendimiento marginal de seguir
#  digitalizando es bajo, cerrar con un ajuste declarado como tal.
#
#  QUÉ SE AJUSTÓ Y CON QUÉ
#  ---------------------------------------------------------------------------
#  Se usan ÚNICAMENTE los 3 anclajes de MESETA/HOMBRO que son limpios y
#  consistentes entre panel A (8.3 mEq/l) y panel B (25 mEq/l):
#      R = 3.5, 5, 10 nm  ->  kT_neutro = 0.465, 0.81, 2.49  (promedio A/B)
#  Se ajustó kT_neutro(R) = A * exp(k*R) (mínimos cuadrados sobre ln(kT), 2
#  parámetros, 3 puntos). NO se mezclaron con los picos de contacto de
#  Fig. 4B: son una magnitud física distinta (pico puntual vs. resistencia
#  sostenida); mezclarlos habría sesgado el ajuste, no lo habría hecho más
#  completo. Resultado: A = 0.2099, k = 0.2497 (1/nm).
#
#  A = 0.2099, k = 0.2497:
#      R=3.5 -> 0.50 kT   (dato: 0.465)
#      R=5   -> 0.73 kT   (dato: 0.81)
#      R=10  -> 2.55 kT   (dato: 2.49)
#      R=15  -> 8.9  kT   (EXTRAPOLADO, sin dato limpio; el piso de Fig.4B es
#                           >=6.93 kT, el ajuste no lo contradice)
#      R=20  -> 31.0 kT   (EXTRAPOLADO, sin dato limpio; el piso de Fig.4B es
#                           >=9.36 kT, el ajuste no lo contradice)
#
#  FACTOR DE CARGA NEGATIVA: cociente negativo/neutro medido en los 4 casos
#  limpios: 2.02 (R=3.5, panel A), 2.24 (R=3.5, panel B), 1.84 (R=5, panel B),
#  2.80 (R=10, panel A, pero ese numerador es un PISO, el cociente real en
#  R=10 es >= 2.80, no exactamente 2.80). Promedio de los 4 = 2.22. Se usa un
#  factor CONSTANTE en R por falta de suficientes puntos para ajustar cómo
#  cambia con el radio, es la simplificación más honesta disponible, no una
#  afirmación de que el cociente sea realmente constante.
#
#  POR QUÉ LA CARGA POSITIVA NO ENTRA EN ESTE ESQUEMA
#  ---------------------------------------------------------------------------
#  Un factor multiplicativo sobre kT_neutro asume que la carga siempre AUMENTA
#  la barrera. Para carga positiva eso es falso: R=3.5 positivo mide ~0.03 kT
#  (casi cancela la barrera estérica) y en panel B las curvas positivas de
#  R=3.5/5/10 forman POZOS atractivos (PMF negativo en parte del recorrido).
#  Multiplicar por un factor <1 constante tampoco serviría, el efecto no es
#  una barrera reducida, es una física distinta (atracción neta en parte del
#  trayecto). No se modela aquí. Quien necesite carga positiva debe usar
#  PMF_HOMBRO_kT_A / _B directamente y tratar la ausencia de barrera como
#  información en sí misma, no rellenarla con este ajuste.
FACTOR_CARGA_NEGATIVA = 2.22   # promedio de 4 cocientes medidos, ver arriba
_A_AJUSTE_kT = 0.2099
_K_AJUSTE_kT = 0.2497          # 1/nm


def kT_hombro_neutro_kT(R_nm):
    """kT_neutro(R) = A*exp(k*R). Ajuste, NO dato; ver cabecera de esta
    sección. Reproduce los 3 anclajes medidos (R=3.5,5,10) a +-0.08 kT;
    fuera de ese rango es EXTRAPOLACIÓN, declarada como tal en el docstring
    de kT_hombro().
    """
    return _A_AJUSTE_kT * np.exp(_K_AJUSTE_kT * np.asarray(R_nm, dtype=float))


def kT_hombro(R_nm, carga="neutro"):
    """kT(R, carga) en la meseta/hombro del glicocálix. AJUSTE, no dato.

    carga: "neutro" (ajuste exponencial), "negativo" (ajuste x factor
    constante 2.22, ver cabecera), o "positivo" (NO soportado: lanza
    ValueError a propósito, ver "POR QUÉ LA CARGA POSITIVA NO ENTRA EN ESTE
    ESQUEMA" arriba; usar PMF_HOMBRO_kT_A/_B directamente para ese caso).

    Devuelve un dict con el valor, si R está dentro del rango medido
    (3.5-10 nm, interpolación) o fuera (extrapolación, menos confiable).
    """
    R_nm = float(R_nm)
    if carga == "positivo":
        raise ValueError(
            "kT_hombro() no soporta carga positiva: no es un múltiplo del "
            "neutro (a veces cancela la barrera, a veces es un pozo "
            "atractivo). Usar PMF_HOMBRO_kT_A / PMF_HOMBRO_kT_B directamente."
        )
    base = float(kT_hombro_neutro_kT(R_nm))
    if carga == "negativo":
        valor = base * FACTOR_CARGA_NEGATIVA
    elif carga == "neutro":
        valor = base
    else:
        raise ValueError(f"carga debe ser 'neutro' o 'negativo', recibido {carga!r}")
    extrapolado = not (3.5 <= R_nm <= 10.0)
    return dict(R_nm=R_nm, carga=carga, kT=valor, extrapolado=extrapolado,
                nota=("fuera del rango medido (3.5-10nm), extrapolación del "
                      "ajuste, no dato" if extrapolado else
                      "dentro del rango medido, ajuste interpolado"))


# =============================================================================
#  MODELO SECUNDARIO: coste osmótico de inserción (estimación de escalado)
# =============================================================================

def energia_insercion_kT(R_nm, xi_nm):
    """dG ~ Pi V = kT (4/3) pi R^3 / xi^3, con Pi ~ kT/xi^3.

    Es una estimación de ESCALADO, no una teoría cerrada. Se conserva porque da
    una idea del ORDEN de la penalización cuando phi ya es 0 y la teoría de
    matriz de fibras solo dice "excluido" sin cuantificar cuánto.
    """
    return (4.0 / 3.0) * np.pi * (np.asarray(R_nm, dtype=float) / float(xi_nm)) ** 3


# =============================================================================
#  DEGRADACIÓN DEL GLICOCÁLIX
#
#  Pregunta: ¿cuánto tendría que degradarse la matriz para que pase un
#  nanotransportador de diámetro d?
#
#  Modelo: si se elimina una fracción f de las fibras, la densidad areal cae a
#  (1 - f), el espaciado centro a centro crece como s = s0 / sqrt(1 - f), y el
#  hueco efectivo pasa a ser Delta = s - 2 rf.
#
#  AVISO IMPORTANTE SOBRE LOS DATOS EXPERIMENTALES. Lo que la literatura mide
#  casi siempre es el GROSOR de la capa, no la densidad areal de fibras. Son
#  magnitudes distintas: unas fibras más cortas con el mismo espaciado dejan el
#  tamiz IGUAL. No hay conversión directa entre "50% menos de grosor" y "50% de
#  fibras eliminadas". Aquí se usa la lectura más GENEROSA posible, que es
#  suponer que un x% menos de grosor equivale a un x% de fibras eliminadas.
#  Es un límite superior del efecto, no una estimación.
# =============================================================================

def hueco_tras_degradacion_nm(f, rf_nm=RF_nm, espaciado_nm=ESPACIADO_nm):
    """Hueco entre fibras tras eliminar una fracción f de ellas."""
    f = float(np.clip(f, 0.0, 0.999999))
    return espaciado_nm / np.sqrt(1.0 - f) - 2.0 * rf_nm


def diametro_admisible_tras_degradacion_nm(f, rf_nm=RF_nm, espaciado_nm=ESPACIADO_nm):
    """Diámetro máximo admisible tras eliminar una fracción f de fibras."""
    D = hueco_tras_degradacion_nm(f, rf_nm, espaciado_nm)
    return 0.0 if D <= 0 else 2.0 * radio_exclusion_nm(rf_nm, D)


def degradacion_necesaria(d_objetivo_nm, rf_nm=RF_nm, espaciado_nm=ESPACIADO_nm):
    """Fracción de fibras que hay que eliminar para admitir un diámetro dado.

    Búsqueda por bisección; la función es monótona creciente en f.
    """
    lo, hi = 0.0, 0.9999999
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if diametro_admisible_tras_degradacion_nm(m, rf_nm, espaciado_nm) < d_objetivo_nm:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


# Degradación MEDIDA en la literatura, para contraste. Todas son reducciones de
# GROSOR, no de densidad de fibras, y ninguna es de cerebro ni de esclerosis
# múltiple. Convergen en torno al 50% en sepsis, que es de los estados
# inflamatorios más severos que existen.
DEGRADACION_MEDIDA = {
    "sepsis, ratón, AFM (266 -> 137 nm)": 0.49,
    "trombina / LPS / TNF-alfa in vitro": 0.50,
    "endotoxina en humano, sublingual (0.60 -> 0.30 um)": 0.50,
}


# =============================================================================
#  LÍMITE GEOMÉTRICO INFERIOR DE UN LIPOSOMA
#
#  Responde a la objeción evidente: "si 40 nm no pasa, hagámoslo más pequeño".
#
#  Un liposoma es una bicapa cerrada sobre sí misma, así que su diámetro externo
#  no puede ser menor que dos espesores de bicapa:
#
#      d_externo = d_núcleo + 2 t_bicapa
#
#  El espesor de bicapa está MEDIDO. Pan et al. 2008, Phys. Rev. Lett. 100,
#  198103, Fig. 3(c): la distancia cabeza-cabeza D_HH va de 36 a 46 Angstrom,
#  es decir 3.6 a 4.6 nm, según el lípido y el contenido de colesterol. El
#  espesor total incluyendo las cabezas polares es algo mayor.
#
#  Consecuencia: incluso un liposoma hipotético con núcleo acuoso NULO mediría
#  entre 7.2 y 9.2 nm de diámetro, que ya está en el límite del glicocálix, y no
#  podría encapsular nada en su interior acuoso. Cualquier liposoma real, con un
#  núcleo utilizable, supera con holgura ese límite.
# =============================================================================

T_BICAPA_nm = [3.6, 4.0, 4.6]   # D_HH medido (Pan et al. 2008, PRL 100:198103, Fig. 3c)


def diametro_liposoma_minimo_nm(d_nucleo_nm, t_bicapa_nm=4.0):
    """Diámetro externo mínimo de un liposoma con un núcleo acuoso dado."""
    return d_nucleo_nm + 2.0 * t_bicapa_nm


def puede_existir_liposoma_que_pase(rf_nm=RF_nm, delta_nm=DELTA_nm):
    """¿Existe algún liposoma que pase el tamiz del glicocálix?

    Compara el límite superior del glicocálix con el límite inferior geométrico
    de un liposoma. Devuelve un dict con el veredicto y el margen.
    """
    d_max_glico = 2.0 * radio_exclusion_nm(rf_nm, delta_nm)
    d_min_lipo = min(diametro_liposoma_minimo_nm(0.0, t) for t in T_BICAPA_nm)
    return dict(d_max_glicocalix_nm=d_max_glico,
                d_min_liposoma_nm=d_min_lipo,
                margen_nm=d_max_glico - d_min_lipo,
                existe=bool(d_min_lipo < d_max_glico),
                nota=("solo con núcleo acuoso nulo, que no encapsularía nada"
                      if d_min_lipo < d_max_glico else "ni siquiera con núcleo nulo"))


# =============================================================================
#  ALEXANDER–DE GENNES: solo comprobación de régimen
# =============================================================================

def regimen_adg(R_nm, L_nm):
    """AdG describe compresión y exige R >> L. Devuelve el cociente y el fallo."""
    ratio = R_nm / L_nm
    if ratio > 1.0:
        v = "compresión: AdG aplica"
    elif ratio > 0.1:
        v = "intermedio: AdG dudoso"
    else:
        v = "PENETRACIÓN: AdG NO aplica"
    return dict(R_nm=R_nm, L_nm=L_nm, R_sobre_L=ratio, veredicto=v)


# =============================================================================
#  VALIDACIÓN
# =============================================================================

def test_glicocalix(verbose=True):
    ok = []

    def chequeo(nombre, cond, detalle=""):
        ok.append(bool(cond))
        if verbose:
            print(f"  [{'OK ' if cond else 'FALLA'}] {nombre}{'  ' + detalle if detalle else ''}")

    if verbose:
        print("=" * 78)
        print(" VALIDACIÓN DEL MÓDULO DE GLICOCÁLIX (v2, contra fuente primaria)")
        print("=" * 78)

    # T1: reproduce sigma del modelo NUEVO publicado por Weinbaum (0.67).
    s = float(coef_reflexion(3.5, RF_nm, DELTA_nm))
    chequeo("T1 sigma(albúmina) con rf=6, Delta=8", abs(s - 0.67) < 0.005,
            f"da {s:.3f} (Weinbaum: 0.67)")

    # T2: reproduce sigma del modelo ANTIGUO publicado por Weinbaum (0.52).
    s = float(coef_reflexion(3.5, RF_GAG_nm, DELTA_GAG_nm))
    chequeo("T2 sigma(albúmina) con rf=0.6, Delta=8", abs(s - 0.52) < 0.005,
            f"da {s:.3f} (Weinbaum: 0.52)")

    # T3: reproduce la fracción de volumen de fibra publicada (c = 0.326).
    c = fraccion_volumen_fibra(RF_nm, DELTA_nm)
    chequeo("T3 fracción de volumen de fibra", abs(c - 0.326) < 0.002,
            f"da {c:.4f} (Weinbaum: 0.326)")

    # T4 · coherencia: en el radio de exclusión, phi debe ser exactamente 0.
    a_max = radio_exclusion_nm()
    chequeo("T4 phi = 0 justo en el radio de exclusión",
            abs(float(coef_particion(a_max))) < 1e-9, f"a_max = {a_max:.2f} nm")

    # T5 · monotonía: phi debe decrecer al crecer el soluto.
    phis = [float(coef_particion(a)) for a in (0.5, 1.0, 2.0, 3.0, 4.0)]
    chequeo("T5 monotonía: phi decrece con el tamaño del soluto",
            all(phis[i] > phis[i + 1] for i in range(len(phis) - 1)))

    # T6 · FALSABILIDAD: un hueco mucho mayor debe dejar pasar lo que el real no.
    g_estrecho = compuerta_glicocalix(15.5, DELTA_nm)
    g_ancho = compuerta_glicocalix(15.5, RF_nm, 200.0)
    chequeo("T6 falsabilidad: con hueco de 200 nm sí pasa, con 8 nm no",
            (not g_estrecho["G5_admisible"]) and g_ancho["G5_admisible"])

    # T7: el régimen del proyecto no es el de AdG (caso menos favorable).
    #      El caso MENOS favorable es el glicocálix más FINO, porque R/L sube.
    #      Se evalúa sobre rana/hámster Y sobre el cerebral, incluido el
    #      envejecido de Shi 2025, que es el más fino de todos (232 nm).
    _todos_L = L_GLICO_nm + L_GLICO_CEREBRAL_nm + [L_GLICO_CEREBRAL_ENVEJECIDO_nm]
    peor = max(regimen_adg(R, L)["R_sobre_L"] for R in (15.5, 20.0) for L in _todos_L)
    chequeo("T7 el régimen del proyecto NO es el de AdG", peor <= 0.15,
            f"R/L máximo = {peor:.4f}")
    # T10: el glicocálix CEREBRAL es más grueso que el de rana/hámster, así que
    #       el trayecto hasta la membrana es MÁS largo de lo que asumía D.3.
    #       OJO: el nombre T9 ya estaba cogido más abajo (liposoma con núcleo de
    #       5 nm). Se detectó al ver la salida, no en las pruebas: dos pruebas
    #       con el mismo nombre pasan igual y solo se nota mirando el informe.
    chequeo("T10 el glicocálix cerebral supera el rango de rana/hámster",
            min(L_GLICO_CEREBRAL_nm) > max(L_GLICO_nm),
            f"{min(L_GLICO_CEREBRAL_nm):.0f} nm frente a {max(L_GLICO_nm):.0f} nm")

    # T8: coherencia del límite geométrico del liposoma.
    chequeo("T8 liposoma de núcleo nulo = 2 espesores de bicapa",
            abs(diametro_liposoma_minimo_nm(0.0, 4.0) - 8.0) < 1e-9)

    # T9: ningún liposoma con núcleo utilizable pasa el tamiz.
    v = puede_existir_liposoma_que_pase()
    d_util = diametro_liposoma_minimo_nm(5.0, 4.0)   # núcleo de 5 nm, mínimo utilizable
    chequeo("T9 un liposoma con núcleo de 5 nm NO pasa el glicocálix",
            d_util > v["d_max_glicocalix_nm"],
            f"{d_util:.1f} nm frente a un límite de {v['d_max_glicocalix_nm']:.1f} nm")

    # T11: el ajuste kT_hombro_neutro_kT reproduce los 3 anclajes medidos.
    anclajes = {3.5: 0.465, 5.0: 0.81, 10.0: 2.49}
    errores = {r: abs(float(kT_hombro_neutro_kT(r)) - v) for r, v in anclajes.items()}
    chequeo("T11 el ajuste exponencial reproduce los anclajes medidos (+-0.1 kT)",
            all(e < 0.1 for e in errores.values()),
            f"errores = {[round(e, 3) for e in errores.values()]}")

    # T12: el ajuste es monótono creciente en R (barrera estérica no baja).
    vals = [float(kT_hombro_neutro_kT(r)) for r in (3.5, 5, 7.5, 10, 15, 20)]
    chequeo("T12 monotonía: kT_hombro_neutro_kT crece con R",
            all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)))

    # T13: la extrapolación a R=15/20 no contradice los pisos de Fig.4B
    #       (PMF_CONTACTO_NEUTRO_kT es una magnitud distinta: pico, no
    #       hombro; pero un hombro nunca debería caer por debajo de su
    #       propio piso medido en otra figura si el piso es de verdad un
    #       límite inferior de la física real).
    piso15 = PMF_CONTACTO_NEUTRO_kT[15.0][0]
    piso20 = PMF_CONTACTO_NEUTRO_kT[20.0][0]
    chequeo("T13 extrapolación R=15/20 no contradice los pisos de Fig.4B",
            float(kT_hombro_neutro_kT(15.0)) >= piso15 and
            float(kT_hombro_neutro_kT(20.0)) >= piso20,
            f"R=15: {float(kT_hombro_neutro_kT(15.0)):.1f} vs piso {piso15}; "
            f"R=20: {float(kT_hombro_neutro_kT(20.0)):.1f} vs piso {piso20}")

    # T14: kT_hombro() marca correctamente interpolación vs extrapolación.
    r_interp = kT_hombro(5.0, "neutro")
    r_extrap = kT_hombro(15.0, "neutro")
    chequeo("T14 kT_hombro() distingue interpolación de extrapolación",
            (not r_interp["extrapolado"]) and r_extrap["extrapolado"])

    # T15: el factor de carga negativa se aplica multiplicativamente.
    neu = kT_hombro(5.0, "neutro")["kT"]
    neg = kT_hombro(5.0, "negativo")["kT"]
    chequeo("T15 factor de carga negativa aplicado correctamente",
            abs(neg / neu - FACTOR_CARGA_NEGATIVA) < 1e-9)

    # T16 · FALSABILIDAD: la carga positiva debe rechazarse explícitamente,
    #       no devolver un número silencioso (no es un múltiplo del neutro).
    rechazo_positivo = False
    try:
        kT_hombro(5.0, "positivo")
    except ValueError:
        rechazo_positivo = True
    chequeo("T16 carga positiva rechazada explícitamente (no modelada)",
            rechazo_positivo)

    if verbose:
        print("-" * 78)
        print(f" RESULTADO: {sum(ok)}/{len(ok)} pruebas superadas")
        print("=" * 78)
    return all(ok)


# =============================================================================
#  INFORME
# =============================================================================

LIPOSOMAS = [("Convencional", 20.0), ("Furtivo/PEG", 15.5), ("Catiónico", 17.5)]


def informe():
    print("=" * 79)
    print(" GLICOCÁLIX: ¿puede el nanotransportador atravesar la matriz?")
    print("=" * 79)
    print(" Fuente primaria: Weinbaum et al., PNAS 100:7988 (2003), 'Transport Model'.")
    print(f" Geometría: radio de fibra {RF_nm} nm, hueco {DELTA_nm} nm,")
    print(f"            espaciado centro a centro {ESPACIADO_nm} nm.")
    print(f" Grosor: {L_GLICO_nm} nm (rana / hámster).")
    print("=" * 79)

    a_max = radio_exclusion_nm()
    print(f"\n### 1. Límite estérico de la matriz intacta\n")
    print(f"  Fracción de volumen de fibra Vf = {fraccion_volumen_fibra():.4f}")
    print(f"  Exclusión total (phi = 0) a partir de radio {a_max:.2f} nm")
    print(f"  ->  DIÁMETRO MÁXIMO ADMISIBLE = {2*a_max:.1f} nm")

    print("\n### 2. Coeficiente de partición por tamaño\n")
    print(f"{'diámetro (nm)':>14} {'radio (nm)':>11} {'phi':>9} {'sigma':>8}   estado")
    print("-" * 62)
    for d in (2, 5, 7, 9, 10, 20, 31, 35, 40, 50, 100):
        a = d / 2.0
        p = float(coef_particion(a))
        s = float(coef_reflexion(a))
        est = "pasa" if p > 0 else "EXCLUIDO"
        marca = "  <- albúmina" if d == 7 else ""
        print(f"{d:14d} {a:11.1f} {p:9.4f} {s:8.3f}   {est}{marca}")
    print("-" * 62)

    print("\n### 3. Los tres diseños del estudio\n")
    print(f"{'diseño':14s} {'Ø (nm)':>8} {'phi':>8}   veredicto")
    print("-" * 55)
    for nombre, R in LIPOSOMAS:
        g = compuerta_glicocalix(R)
        print(f"{nombre:14s} {2*R:8.1f} {g['phi']:8.4f}   "
              f"{'admisible' if g['G5_admisible'] else 'EXCLUIDO'}")
    print("-" * 55)

    print("\n### 4. Sensibilidad al hueco entre fibras (no medido en cerebro)\n")
    print(f"{'hueco (nm)':>11} {'Vf':>9} {'Ø máx (nm)':>12}   ¿pasan los diseños?")
    print("-" * 62)
    for D in DELTA_BARRIDO_nm:
        dmax = 2 * radio_exclusion_nm(RF_nm, D)
        pasan = "sí" if dmax >= 40.0 else ("alguno" if dmax >= 31.0 else "NINGUNO")
        print(f"{D:11.1f} {fraccion_volumen_fibra(RF_nm, D):9.4f} {dmax:12.1f}   {pasan}")
    print("-" * 62)
    print(" Haría falta un hueco de ~40 nm para que los diseños pasen. El valor")
    print(" medido es 8 nm, y en cerebro el glicocálix es MÁS denso, no menos.")

    print("\n### 5. ¿Y si hacemos el liposoma más pequeño?\n")
    print(f"  Un liposoma es una bicapa cerrada: d_externo = d_núcleo + 2 t_bicapa.")
    print(f"  Espesor de bicapa medido (Pan 2008, Fig. 3c): {T_BICAPA_nm} nm\n")
    print(f"{'núcleo acuoso (nm)':>20} | " + " | ".join(f"t = {t} nm".rjust(12)
                                                       for t in T_BICAPA_nm))
    print("-" * 66)
    for c in (0, 2, 5, 10, 20):
        fila = [f"{diametro_liposoma_minimo_nm(c, t):9.1f} nm" for t in T_BICAPA_nm]
        print(f"{c:20d} | " + " | ".join(s.rjust(12) for s in fila))
    print("-" * 66)
    v = puede_existir_liposoma_que_pase()
    print(f"  Límite superior del glicocálix ....... {v['d_max_glicocalix_nm']:.1f} nm")
    print(f"  Liposoma más pequeño concebible ...... {v['d_min_liposoma_nm']:.1f} nm "
          f"(núcleo nulo)")
    print(f"  Margen ............................... {v['margen_nm']:+.1f} nm")
    print(f"  -> {v['nota']}")

    print("\n### 6. Por qué no se usa Alexander–de Gennes\n")
    print(f"{'R (nm)':>8} {'L (nm)':>9} {'R/L':>9}   veredicto")
    print("-" * 58)
    for R in (15.5, 20.0):
        for L in L_GLICO_nm:
            r = regimen_adg(R, L)
            print(f"{R:8.1f} {L:9.0f} {r['R_sobre_L']:9.4f}   {r['veredicto']}")
    print("-" * 58)

    print("\n" + "=" * 79)
    print(" LECTURA: el resultado más fuerte y más incómodo del proyecto")
    print("=" * 79)
    print(" 1. Con la geometría medida, el glicocálix intacto excluye estéricamente")
    print("    todo lo que pase de ~9 nm de diámetro. Tus liposomas miden 31-40 nm.")
    print("    Quedan EXCLUIDOS por un factor de 4.")
    print()
    print(" 2. No es sorprendente visto de frente: la función fisiológica del")
    print("    glicocálix es ser un tamiz que retiene la albúmina (7 nm). Un objeto")
    print("    cuatro veces mayor que la albúmina no lo atraviesa por difusión.")
    print()
    print(" 3. LO QUE ESTO NO DEMUESTRA. El modelo es una red ORDENADA, RÍGIDA y")
    print("    ESTÁTICA, ajustada a capilar de mesenterio de RANA. El glicocálix real")
    print("    es heterogéneo, dinámico, se degrada en patología y tiene defectos.")
    print("    Y es un hecho experimental que hay nanopartículas que llegan al")
    print("    endotelio cerebral. Así que o hay aberturas transitorias mayores, o")
    print("    entran por zonas adelgazadas o dañadas, o comprimen la capa.")
    print()
    print(" 4. LA PREGUNTA DEL PROYECTO CAMBIA. Deja de ser '¿qué diseño adhiere")
    print("    mejor?' y pasa a ser '¿por qué vía atraviesa algo de 30-40 nm una capa")
    print("    que filtra a 9 nm?'. Esa pregunta es más interesante y es exactamente")
    print("    la tercera pregunta que hay que hacerle al doctor.")
    print()
    print(" 5. FALTA la electrostática. El glicocálix es fuertemente ANIÓNICO. Una")
    print("    partícula catiónica no solo es excluida estéricamente: puede quedar")
    print("    ATRAPADA. Eso empeoraría, no mejoraría, la cláusula de carga de la")
    print("    hipótesis.")
    print("=" * 79)


if __name__ == "__main__":
    test_glicocalix()
    print()
    informe()
