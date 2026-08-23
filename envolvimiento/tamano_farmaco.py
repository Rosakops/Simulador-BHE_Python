#!/usr/bin/env python3
"""Tamaño molecular del fingolimod y del FTY720-fosfato (tarea C9).

POR QUÉ EXISTE ESTE ARCHIVO
---------------------------
`rutas.py` llevaba `farmaco_diametro_nm = 1.0` con el comentario
«fingolimod, 307 Da: ~1 nm  [sin fuente, tarea C9]». Era una estimación a ojo
dentro de una cadena de cálculo, que es justo lo que la regla de método del
proyecto prohíbe.

Este script NO mide nada: DERIVA el tamaño de la estructura química, con un
método declarado y reproducible. El resultado es un cálculo, no un dato
experimental, y así hay que citarlo. Se ejecuta con:

    python3 tamano_farmaco.py

MÉTODO
------
1. Estructura de partida en SMILES, escrita desde la fórmula del expediente:
   fingolimod = 2-amino-2-[2-(4-octilfenil)etil]propano-1,3-diol.
   Comprobación cruzada: RDKit devuelve C19H33NO2 y 307.48 Da, que son
   exactamente la fórmula y la masa del expediente JC155. Si algún día no
   coincide, el SMILES está mal.
2. 50 confórmeros con ETKDGv3 (semilla fija, resultado reproducible) y
   optimización MMFF94.
3. Volumen de van der Waals por rejilla de 0.15 Å, promediado sobre los 50.
4. Dos tamaños, porque la molécula NO es esférica:
     · diámetro de la esfera de igual volumen  -> D = 2·(3V/4π)^(1/3)
     · dimensión máxima extremo a extremo, promediada sobre los confórmeros
5. Se hace lo mismo con el FTY720-fosfato, que es la especie que exporta SPNS2
   y la que de verdad tiene que difundir (ver F.1).

QUÉ NÚMERO USA EL MODELO Y POR QUÉ
----------------------------------
El modelo usa la DIMENSIÓN MÁXIMA, no el diámetro esférico. Es la elección
conservadora: si la molécula pasa con su dimensión más larga, pasa seguro. Y es
coherente con la lección de Cheng 2019 anotada en la ficha de C.2, con una
nanovarilla, los propios autores atribuyeron el bloqueo al EJE LARGO y no al
diámetro. Cuál de las dos dimensiones gobierna el paso por un poro no está
resuelto para una molécula flexible, así que se toma la peor.

DE PASO RESUELVE UNA DISCREPANCIA
---------------------------------
La crítica al simulador de acople (2026-08-12) señalaba que el fingolimod medía
«~1.5-2 nm» en las fichas de Jhovan y «~1 nm» en el proyecto, sin fuente ninguno.
No eran dos valores en conflicto: son las DOS dimensiones del mismo objeto.
El ~1 nm es el diámetro esférico equivalente (0.86) y el ~1.5-2 nm es la
dimensión máxima (1.68). Los dos estaban bien y a la vez mal explicados.

SALVEDADES
----------
· Es geometría en el vacío. NO es un radio hidrodinámico medido, que sería mayor
  por la capa de solvatación. Para cerrar C9 del todo haría falta un radio
  hidrodinámico experimental o una estructura cristalina.
· El volumen de van der Waals depende del juego de radios atómicos de RDKit.
· El fosfato se modela NEUTRO. A pH fisiológico está desprotonado y cargado, lo
  que cambia la solvatación pero no apreciablemente el volumen de van der Waals.
"""

import math

SMILES = {
    "fingolimod neutro (FTY720)": "CCCCCCCCc1ccc(CCC(N)(CO)CO)cc1",
    "FTY720-fosfato":             "CCCCCCCCc1ccc(CCC(N)(CO)COP(=O)(O)O)cc1",
}

N_CONFORMEROS = 50
SEMILLA = 0xf00d
REJILLA_A = 0.15

# --- resultados de la última ejecución, para poder comprobarlos sin RDKit ----
# Fingolimod neutro:  C19H33NO2, 307.48 Da, V = 331.0 A^3
#                     D esfera = 0.858 nm · dimensión máxima = 1.683 nm
# FTY720-fosfato:     C19H34NO5P, 387.46 Da, V = 377.9 A^3
#                     D esfera = 0.897 nm · dimensión máxima = 1.777 nm
FINGOLIMOD_D_ESFERA_nm = 0.858
FINGOLIMOD_D_MAXIMO_nm = 1.683
FOSFATO_D_ESFERA_nm = 0.897
FOSFATO_D_MAXIMO_nm = 1.777


def calcular(smiles, n=N_CONFORMEROS):
    """Devuelve (fórmula, masa, volumen A^3, D esfera nm, D máximo nm)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    import numpy as np

    plano = Chem.MolFromSmiles(smiles)
    m = Chem.AddHs(plano)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEMILLA
    cids = AllChem.EmbedMultipleConfs(m, numConfs=n, params=ps)
    AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)

    vols = [AllChem.ComputeMolVolume(m, confId=c, gridSpacing=REJILLA_A)
            for c in cids]
    v = sum(vols) / len(vols)
    d_esfera = 2.0 * (3.0 * v / (4.0 * math.pi)) ** (1 / 3) / 10.0

    maximos = []
    for c in cids:
        pos = m.GetConformer(c).GetPositions()
        dist = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1))
        maximos.append(dist.max())
    d_max = sum(maximos) / len(maximos) / 10.0

    return (rdMolDescriptors.CalcMolFormula(plano), Descriptors.MolWt(plano),
            v, d_esfera, d_max, min(maximos) / 10.0, max(maximos) / 10.0)


def informe():
    print()
    print("=" * 79)
    print(" TAMAÑO MOLECULAR DEL FÁRMACO (tarea C9) · CÁLCULO, no medida")
    print("=" * 79)
    print(f" {N_CONFORMEROS} confórmeros ETKDGv3 (semilla fija) + MMFF94 · "
          f"volumen vdW por rejilla de {REJILLA_A} Å")
    print()
    for nombre, smi in SMILES.items():
        f, mw, v, de, dm, dmin, dmax = calcular(smi)
        print(f" {nombre}")
        print(f"   fórmula                {f}")
        print(f"   masa molar             {mw:.2f} Da")
        print(f"   volumen de van der Waals {v:.1f} Å³")
        print(f"   D de esfera equivalente  {de:.3f} nm")
        print(f"   dimensión máxima         {dm:.3f} nm  "
              f"(entre {dmin:.3f} y {dmax:.3f} según el confórmero)")
        print()
    print(" El modelo usa la DIMENSIÓN MÁXIMA, que es la elección conservadora.")
    print(" Los dos valores quedan muy por debajo de los 38 nm de Thorne, así que")
    print(" el veredicto de la compuerta del fármaco no depende de cuál se elija.")
    print("=" * 79)
    print()


if __name__ == "__main__":
    informe()
