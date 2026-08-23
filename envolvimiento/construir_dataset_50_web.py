#!/usr/bin/env python3
# =============================================================================
#  construir_dataset_50_web.py · última día, 2026-08-18
# =============================================================================
#  Página estática independiente para los 50 liposomas sintéticos generados en
#  dataset_50_liposomas.py. Sigue la misma regla de oro que construir_web.py:
#  esta página NO calcula nada, solo ejecuta rutas.py (evaluar, figuras) y
#  escribe los números y las imágenes que acaba de obtener.
#
#  Tabla: nombre + parámetros físico-químicos por fila.
#  Al desplegar una fila (<details>, sin JS): veredicto y desglose de
#  compuertas por las 4 rutas + las 3 figuras propias de ESE liposoma
#  (ventanas/matriz/recorrido con catalogo=[d], igual que hace
#  figuras_liposoma_separadas() para reales/teóricos).
#
#  SALE: web/dataset_50.html, web/img/dataset_50/liposoma_NN_*.png
# =============================================================================

import html
from pathlib import Path

import rutas as R
from dataset_50_liposomas import generar

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
SALIDA = RAIZ / "web"
IMG_DIR = SALIDA / "img" / "dataset_50"


def _fmt(v):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.3g}"
    return str(v)


def _color(estado):
    return {"PASA": "#2e7d32", "FALLA": "#c62828",
            "DESCONOCIDA": "#9e9e9e"}.get(estado, "#333")


def _color_ruta(v):
    return {"NO EXCLUIDA": "#2e7d32", "EXCLUIDA": "#c62828",
            "NO EVALUABLE": "#9e9e9e"}.get(v, "#333")


def _tabla_compuertas(resultados):
    filas = []
    for r in resultados:
        filas.append(
            f'<tr><td>{html.escape(r.compuerta)}</td>'
            f'<td style="color:{_color(r.estado)};font-weight:600">{r.estado}</td>'
            f'<td>{_fmt(r.valor)}</td><td>{_fmt(r.umbral)}</td>'
            f'<td>{html.escape(r.unidad or "")}</td>'
            f'<td class="fuente">{html.escape(r.fuente or "")}</td></tr>'
        )
    return ("<table class='compuertas'><thead><tr><th>Compuerta</th><th>Estado</th>"
            "<th>Valor</th><th>Umbral</th><th>Unidad</th><th>Fuente</th></tr></thead>"
            f"<tbody>{''.join(filas)}</tbody></table>")


def construir():
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    disenos = generar()

    filas_tabla = []
    detalles = []

    for i, d in enumerate(disenos, start=1):
        veredictos = R.evaluar(d)

        prefijo = str(IMG_DIR / f"liposoma_{i:02d}")
        R.figuras(prefijo=prefijo, catalogo=[d], incluir_ventanas=True)
        img_rel = f"img/dataset_50/liposoma_{i:02d}"

        celdas_rutas = "".join(
            f'<td style="color:{_color_ruta(v)};font-weight:600">{v}</td>'
            for v, _ in veredictos.values()
        )
        filas_tabla.append(
            f'<tr><td>{html.escape(d.nombre)}</td>'
            f'<td>{d.diametro_nm:.1f}</td><td>{d.zeta_mV:+.2f}</td>'
            f'<td>{d.peg_nm:.2f}</td>{celdas_rutas}</tr>'
        )

        bloques_rutas = []
        for nombre_ruta, (v, resultados) in veredictos.items():
            bloques_rutas.append(
                f'<h4>{html.escape(nombre_ruta)} · '
                f'<span style="color:{_color_ruta(v)}">{v}</span></h4>'
                f'{_tabla_compuertas(resultados)}'
            )

        detalles.append(f'''
<details class="liposoma">
<summary>{html.escape(d.nombre)} · Ø {d.diametro_nm:.1f} nm · ζ {d.zeta_mV:+.2f} mV ·
PEG {d.peg_nm:.2f} nm</summary>
<div class="detalle">
  <div class="figuras3">
    <img src="{img_rel}_ventanas.png" alt="ventanas">
    <img src="{img_rel}_matriz.png" alt="matriz">
    <img src="{img_rel}_recorrido.png" alt="recorrido">
  </div>
  {''.join(bloques_rutas)}
</div>
</details>''')

    nombres_rutas = list(R.RUTAS.keys())
    cab_rutas = "".join(f"<th>{html.escape(n)}</th>" for n in nombres_rutas)

    html_out = f'''<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>Dataset de 50 liposomas sintéticos · JC155</title>
<link rel="stylesheet" href="estilo.css">
<style>
table.principal {{ border-collapse: collapse; width: 100%; margin-bottom: 2rem; }}
table.principal th, table.principal td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: center; font-size: 0.9rem; }}
table.principal th {{ background: #f4f4f4; }}
table.compuertas {{ border-collapse: collapse; width: 100%; margin: 0.5rem 0 1rem; }}
table.compuertas th, table.compuertas td {{ border: 1px solid #eee; padding: 4px 8px; font-size: 0.82rem; text-align: left; }}
table.compuertas td.fuente {{ font-size: 0.75rem; color: #666; }}
details.liposoma {{ border: 1px solid #ddd; border-radius: 6px; margin-bottom: 8px; padding: 6px 10px; }}
details.liposoma summary {{ cursor: pointer; font-weight: 600; }}
.figuras3 {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0; }}
.figuras3 img {{ max-width: 32%; border: 1px solid #ccc; }}
.aviso {{ background:#fff3cd; border:1px solid #ffe08a; padding:8px 12px; border-radius:6px; margin-bottom:1rem; font-size:0.9rem; }}
</style>
</head>
<body>
<div class="envoltorio">
<h1>Dataset de 50 liposomas sintéticos</h1>
<p class="aviso"><strong>SINTÉTICOS.</strong> Los 50 diseños son combinaciones
aleatorias (semilla 42) de diámetro, ζ y PEG dentro del rango que cubren los
diseños reales ya validados del proyecto (Mao 2014, Gong 2022, Chow 2025,
Muselman 2026). No son mediciones ni predicciones citables, igual que
CATALOGO_TEORICO. La carga útil (G.2) queda DESCONOCIDA en todos, por decisión
vigente del 2026-08-13.</p>

<table class="principal">
<thead><tr><th>Liposoma</th><th>Ø nm</th><th>ζ mV</th><th>PEG nm</th>{cab_rutas}</tr></thead>
<tbody>
{''.join(filas_tabla)}
</tbody>
</table>

<h2>Detalle por liposoma</h2>
<p>Clic sobre cada fila para desplegar sus 3 figuras y el desglose de compuertas por ruta.</p>
{''.join(detalles)}

</div>
</body>
</html>'''

    (SALIDA / "dataset_50.html").write_text(html_out, encoding="utf-8")
    print(f"  {SALIDA / 'dataset_50.html'}")
    print(f"  {len(disenos) * 3} figuras en {IMG_DIR}")


if __name__ == "__main__":
    construir()
