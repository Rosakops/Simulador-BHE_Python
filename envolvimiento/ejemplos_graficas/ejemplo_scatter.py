import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_PATH = "dataset_50_liposomas.csv"

with open(CSV_PATH, encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

RUTAS_COLS = {
    "A": "A · pasiva (adhesión y envolvimiento)",
    "B": "B · celular (macrófago de Troya)",
    "C": "C · adsortiva (carga positiva)",
    "D": "D · mediada por receptor",
}
TITULOS = {
    "A": "A · pasiva",
    "B": "B · celular",
    "C": "C · adsortiva",
    "D": "D · mediada por receptor",
}

diametro = [float(r["diametro_nm"]) for r in rows]
zeta = [float(r["zeta_mV"]) for r in rows]

xmin, xmax = 20, 720
ymin, ymax = -31, 9

fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
fig.patch.set_facecolor("white")

handles = None
for ax, letra in zip(axes.flat, ["A", "B", "C", "D"]):
    col = RUTAS_COLS[letra]
    x_excl, y_excl, x_noeval, y_noeval = [], [], [], []
    for r, x, y in zip(rows, diametro, zeta):
        if r[col] == "EXCLUIDA":
            x_excl.append(x); y_excl.append(y)
        else:
            x_noeval.append(x); y_noeval.append(y)

    if letra == "A":
        ax.axvspan(60, 80, color="#f0a500", alpha=0.12, zorder=0)
        ax.text(70, ymax - 0.8, "caveola\n60–80 nm", ha="center", va="top",
                fontsize=6.5, color="#a67c00")

    ax.axhline(0, color="#888888", lw=0.8, zorder=1)
    h1 = ax.scatter(x_excl, y_excl, color="#e02424", s=32,
                     label="EXCLUIDA", zorder=3)
    h2 = ax.scatter(x_noeval, y_noeval, color="#2255cc", s=32,
                     label="NO EVALUABLE", zorder=3)
    handles = [h1, h2]

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(True, alpha=0.3)
    ax.set_title(f"{TITULOS[letra]}  (excl. {len(x_excl)} · no eval. {len(x_noeval)})",
                 fontsize=10.5)

# Ejes compartidos: solo se rotulan los de la fila/columna de borde,
# igual que "Diámetro" solo aparece una vez abajo y "Potencial zeta" una
# vez a la izquierda.
fig.text(0.5, 0.035, "Diámetro (nm)", ha="center", fontsize=12, color="#c0392b")
fig.text(0.02, 0.5, "Potencial zeta (mV)", va="center", rotation="vertical",
         fontsize=12, color="#c0392b")

fig.suptitle("Diámetro vs Potencial zeta por ruta de entrada a la BHE\n"
             "(50 liposomas sintéticos)", fontsize=14, fontweight="bold", y=0.995)

# Leyenda ÚNICA, compartida por las 4 rutas (mismo código de color en
# todas), colocada fuera de los ejes, debajo de toda la figura -- igual
# que los rótulos de los ejes.
fig.legend(handles=handles, labels=["EXCLUIDA", "NO EVALUABLE"],
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02),
           frameon=False, fontsize=10)

fig.tight_layout(rect=[0.035, 0.06, 1, 0.96])
fig.savefig("ejemplo_scatter_4rutas.png", dpi=160, bbox_inches="tight")
print("OK")
