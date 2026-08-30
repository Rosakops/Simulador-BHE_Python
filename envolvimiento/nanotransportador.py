#!/usr/bin/env python3
# =============================================================================
#  MODELO DE CLASES DE NANOTRANSPORTADOR · NeuroCross
#  Tarea de ARQUITECTURA (no de datos ni de física nueva).
# =============================================================================
#
#  QUÉ ES ESTO
#  -----------
#  El simulador evaluaba candidatos con compuertas físicas hardcodeadas para
#  liposoma, heredadas del repo viejo de SENACYT. Este módulo introduce la
#  estructura que permite soportar MÚLTIPLES clases de nanotransportador
#  (orgánico, inorgánico, híbrido, con subtipos dentro de cada uno) sin
#  reescribir el sistema cada vez que llegue un subtipo nuevo.
#
#  Modelo de TRES NIVELES:
#      categoría  (orgánico | inorgánico | híbrido)
#          -> subtipo       (ej. liposoma, dendrímero — lista ABIERTA)
#              -> parámetros característicos de ese subtipo
#
#  HOY solo existe UN subtipo soportado: "liposoma", dentro de "organico".
#  La lista real de subtipos la están investigando Jhovan y Kiel aparte y
#  llega por partes vía el protocolo Bridge (Bridge/pending.md). Este módulo
#  NO adivina esa lista: agregar un subtipo nuevo es responsabilidad de una
#  tarea futura, validada por Jhovan, no de este archivo.
#
#  REGLA DURA DEL PROYECTO (irrenunciable)
#  ----------------------------------------------------------------------------
#  Cada compuerta física declara, como metadato explícito en el código (el
#  decorador `compuerta_para` de aquí abajo), para qué subtipo(s) fue
#  calibrada/validada. Si se evalúa un candidato de un subtipo para el que
#  una compuerta NO está validada, esa compuerta devuelve DESCONOCIDA
#  automáticamente -- nunca PASA ni FALLA, sin excepción. Verificado por
#  test_nanotransportador() (T5) y por rutas.validar_contra_experimentos()
#  (F28, prueba de extremo a extremo sobre TODAS las rutas del simulador).
# =============================================================================

import functools
from dataclasses import dataclass
from typing import Optional

# =============================================================================
#  VOCABULARIO COMPARTIDO DE LAS COMPUERTAS
#  Vive aquí (y no en rutas.py) para que este módulo no dependa de rutas.py:
#  rutas.py importa de aquí, nunca al revés. Evita el ciclo de importación
#  que tendría el decorador compuerta_para() si necesitara construir un
#  Resultado importándolo desde rutas.py.
# =============================================================================

PASA, FALLA, DESCONOCIDA = "PASA", "FALLA", "DESCONOCIDA"


@dataclass
class Resultado:
    """Resultado de una compuerta. Ver docstring extendido en rutas.py."""
    compuerta: str
    estado: str
    valor: Optional[float] = None
    umbral: Optional[float] = None
    unidad: str = ""
    margen: Optional[float] = None
    fuente: str = ""
    motivo: str = ""
    advertencia: str = ""


# =============================================================================
#  REGISTRO: categoría -> {subtipo: metadata}
#  ÚNICA fuente de verdad de qué categorías/subtipos existen HOY. Lista
#  ABIERTA y extensible por diseño (agregar una fila es agregar un subtipo),
#  pero esta tarea NO agrega ninguno más allá de lo que ya existía: liposoma.
# =============================================================================

REGISTRO = {
    "organico": {
        "liposoma": {
            # Documentación de qué campos de Diseno (rutas.py) son
            # característicos de este subtipo. No aplica validación de tipos
            # por sí sola: es lo que un futuro subtipo debe declarar también.
            "parametros": [
                "diametro_nm", "zeta_mV", "peg_nm", "farmaco_diametro_nm",
            ],
            "nota": ("único subtipo con compuertas calibradas hoy (decisión "
                     "de Jhovan, 2026-08-24: simulador exclusivo de "
                     "liposomas). Ver rutas.g_transportador_fabricable()."),
        },
    },
    "inorganico": {},
    "hibrido": {},
}


class SubtipoNoSoportado(ValueError):
    """El subtipo (o la categoría) pedidos no están soportados todavía.

    Nunca se debe capturar esta excepción para "seguir de todos modos" con
    parámetros libres: el punto de esta excepción es impedir justamente eso.
    """


def categorias_disponibles():
    """Lista de categorías reconocidas (aunque no tengan subtipos todavía)."""
    return sorted(REGISTRO)


def subtipos_disponibles(categoria):
    """Subtipos SOPORTADOS hoy dentro de una categoría.

    Lanza SubtipoNoSoportado si la categoría en sí no se reconoce.
    """
    if categoria not in REGISTRO:
        raise SubtipoNoSoportado(
            f"categoría no reconocida: '{categoria}'. Disponibles: "
            f"{', '.join(categorias_disponibles())}")
    return sorted(REGISTRO[categoria])


def validar_categoria_subtipo(categoria, subtipo):
    """Punto ÚNICO de validación del flujo de entrada del usuario.

    Primero se pregunta la CATEGORÍA, luego se valida el SUBTIPO dentro de
    esa categoría. Si el subtipo pedido no está soportado, lo dice de forma
    explícita en vez de dejar que el candidato se construya de todos modos
    con parámetros libres.

    Devuelve la metadata del subtipo si es válido; lanza SubtipoNoSoportado
    si no.
    """
    disponibles = subtipos_disponibles(categoria)   # valida la categoría
    if subtipo not in disponibles:
        if disponibles:
            detalle = f"soportados hoy en '{categoria}': {', '.join(disponibles)}"
        else:
            detalle = f"'{categoria}' no tiene ningún subtipo soportado todavía"
        raise SubtipoNoSoportado(
            f"subtipo no soportado todavía: '{subtipo}' ({detalle}). No se "
            "evalúa con parámetros libres ni con compuertas de otro subtipo.")
    return REGISTRO[categoria][subtipo]


# =============================================================================
#  DECORADOR: metadato de calibración + cumplimiento automático de la regla dura
# =============================================================================

def compuerta_para(*subtipos_validados, nombre=None):
    """Declara, como metadato explícito en el código, para qué subtipo(s)
    fue calibrada/validada una compuerta física.

    Uso:
        @compuerta_para("liposoma")
        def g_mi_compuerta(d: Diseno):
            ...

    Si `d.clase` (el subtipo del diseño evaluado, ver Diseno.subtipo en
    rutas.py) NO está en `subtipos_validados`, la función decorada NUNCA
    ejecuta su lógica interna: el decorador corta antes y devuelve
    DESCONOCIDA directamente. Esto es lo que hace la regla dura IRROMPIBLE
    por construcción, en vez de depender de que cada compuerta se acuerde de
    chequear la clase a mano (que es justo lo que fallaba antes de esta
    tarea: la mayoría de las compuertas ni miraban `d.clase`).

    `nombre` es solo para la etiqueta del Resultado devuelto cuando el
    subtipo NO calza; si se omite se usa el nombre de la función. No afecta
    al nombre que la compuerta real le pone a su propio Resultado cuando sí
    se ejecuta.
    """
    subtipos_validados = frozenset(subtipos_validados)
    if not subtipos_validados:
        raise ValueError("compuerta_para() necesita declarar al menos un subtipo")

    def decorador(fn):
        etiqueta = nombre or fn.__name__

        @functools.wraps(fn)
        def envoltura(d, *args, **kwargs):
            if d.clase not in subtipos_validados:
                return Resultado(
                    etiqueta, DESCONOCIDA, motivo=(
                        f"compuerta no calibrada ni validada para el subtipo "
                        f"'{d.clase}'; validada solo para: "
                        f"{', '.join(sorted(subtipos_validados))}"))
            return fn(d, *args, **kwargs)

        envoltura.subtipos_validados = subtipos_validados
        return envoltura

    return decorador


# =============================================================================
#  VALIDACIÓN
# =============================================================================

def test_nanotransportador(verbose=True):
    """Prueba el registro, el flujo de entrada y la regla dura del decorador."""
    ok = []

    def chequeo(nombre, cond, detalle=""):
        ok.append(bool(cond))
        if verbose:
            print(f"  [{'OK ' if cond else 'FALLA'}] {nombre}{'  ' + detalle if detalle else ''}")

    if verbose:
        print("=" * 78)
        print(" VALIDACIÓN DEL MODELO DE CLASES DE NANOTRANSPORTADOR")
        print("=" * 78)

    # T1: el registro hoy solo tiene liposoma soportado, dentro de orgánico.
    chequeo("T1 'liposoma' está soportado en la categoría 'organico'",
            "liposoma" in subtipos_disponibles("organico"))
    chequeo("T2 ningún subtipo más está soportado hoy (no se hardcodeó nada extra)",
            sum(len(v) for v in REGISTRO.values()) == 1,
            f"{sum(len(v) for v in REGISTRO.values())} subtipo(s) en el registro")

    # T3: pedir un subtipo no soportado (pero de categoría válida) se rechaza
    #     de forma explícita, no se deja pasar con parámetros libres.
    lanzo = False
    try:
        validar_categoria_subtipo("organico", "dendrimero")
    except SubtipoNoSoportado:
        lanzo = True
    chequeo("T3 subtipo no soportado se rechaza explícitamente (dendrímero)", lanzo)

    # T4: pedir una categoría que no existe también se rechaza explícitamente.
    lanzo = False
    try:
        validar_categoria_subtipo("marciano", "x")
    except SubtipoNoSoportado:
        lanzo = True
    chequeo("T4 categoría no reconocida se rechaza explícitamente", lanzo)

    # T5 · REGLA DURA: una compuerta decorada nunca da PASA/FALLA para un
    #     subtipo fuera de su lista declarada, sin importar lo que su lógica
    #     interna habría decidido con esos mismos números. Se prueba con una
    #     compuerta de mentira que SIEMPRE da PASA si llega a ejecutarse, para
    #     demostrar que el corte es del decorador, no de la lógica interna.
    class _DisenoFalso:
        def __init__(self, clase):
            self.clase = clase

    @compuerta_para("liposoma", nombre="Compuerta de prueba")
    def _compuerta_de_prueba(d):
        return Resultado("Compuerta de prueba", PASA)

    r_valido = _compuerta_de_prueba(_DisenoFalso("liposoma"))
    r_ajeno = _compuerta_de_prueba(_DisenoFalso("dendrimero"))
    chequeo("T5a subtipo validado SÍ ejecuta la lógica interna",
            r_valido.estado == PASA)
    chequeo("T5b subtipo NO validado da DESCONOCIDA, no PASA -- sin excepción",
            r_ajeno.estado == DESCONOCIDA, f"da {r_ajeno.estado}")

    if verbose:
        print("-" * 78)
        print(f" RESULTADO: {sum(ok)}/{len(ok)} pruebas superadas")
        print("=" * 78)
    return all(ok)


if __name__ == "__main__":
    test_nanotransportador()
