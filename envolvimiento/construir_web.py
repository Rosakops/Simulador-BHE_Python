#!/usr/bin/env python3
# =============================================================================
#  construir_web.py  —  proyecto BHE / SENACYT
#
#  Genera una web estática con los resultados del simulador.
#
#  REGLA DE ORO DE ESTE ARCHIVO: la web NO calcula nada. Este script ejecuta el
#  simulador de verdad (rutas.py, glicocalix.py, envolvimiento_core.py) y
#  escribe el HTML con los números que acaba de obtener. Por eso la página no
#  puede desviarse del modelo: si el modelo cambia, se regenera y ya.
#
#  Este archivo NO modifica ningún módulo del simulador. Solo importa y lee.
#
#  Lo único que va escrito a mano aquí son los bloques marcados con
#  "MANTENIDO A MANO": bibliografía, recuento de la validación y pendientes,
#  porque viven en los .md y no en el código. Están señalados en la propia web.
#
#  USO:  python3 construir_web.py        (o: sh correr.sh web)
#  SALE: web/index.html, web/estilo.css, web/fichas/*.html y los PNG copiados.
# =============================================================================

import html
import os
import re
import shutil
import datetime
from pathlib import Path

import rutas as R
import glicocalix as G
from dataset_50_liposomas import generar as _generar_dataset_50

AQUI = Path(__file__).resolve().parent
_BUILD_TS = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
RAIZ = AQUI.parent
SALIDA = RAIZ / "web"


# =============================================================================
#  MANTENIDO A MANO  —  lo que no vive en el código
# =============================================================================

PROYECTO = dict(
    expediente="JC155 · SENACYT",
    titulo=("Modelado computacional de parámetros físico-químicos de "
            "nanopartículas para la permeabilidad de la barrera hematoencefálica"),
    subtitulo="Un estudio in silico como aproximación a terapias para desmielinización",
    investigadores="Jhovan Watson · Astrid Torres",
    farmaco="Fingolimod (FTY720)",
    enfermedad="Esclerosis múltiple",
)

VALIDACION = dict(
    fuente="verificacion/resultados_validacion.md",
    fecha="2026-08-10",
    conteo=[("C1 · compatible con el modelo", 5),
            ("C2 · no atraviesa una BHE intacta", 1),
            ("C3 · barrera comprometida", 1),
            ("C4 · contraejemplo confirmado", 0),
            ("C4 · contraejemplo sin resolver", 0)],
    veredicto=("Se aplica el primer supuesto de la sección 7 del protocolo: "
               "la predicción P1 SOBREVIVE. No queda ningún punto donde el "
               "modelo pueda caerse por contraejemplo."),
)

PENDIENTES = [
    ("G.6", "baja", "Densidad de PLGA por picnometría: la de Parker es una "
     "derivación acústica que sus propios autores llaman estimación."),
    ("C8", "media", "La compuerta de difusión del fármaco no discrimina: aprueba "
     "cualquier molécula por debajo de 38 nm."),
    ("C9", "CERRADA", "El «~1 nm» del fingolimod era una estimación sin fuente. "
     "Sustituido el 2026-08-13 por un valor DERIVADO y reproducible: 1.683 nm de "
     "dimensión máxima (0.858 nm de esfera equivalente), 50 confórmeros "
     "ETKDGv3 + MMFF94, en <code>tamano_farmaco.py</code>. Es un cálculo, no una "
     "medida: falta un radio hidrodinámico experimental."),
    ("G.2", "alta", "Carga útil. La compuerta EXISTE desde el 2026-08-13 en las "
     "cuatro rutas y es DESCONOCIDA permanente. Faltan tres números: carga de un "
     "dendrímero de generación alta con fingolimod, moléculas de fingolimod por "
     "liposoma y dosis necesaria en parénquima."),
    ("G.3", "media", "Tamaño del PAMAM en agua a pH fisiológico. Prosa mide en "
     "metanol y el techo queda a 0.40 nm del umbral de envolvimiento."),
    ("V-7", "media", "Tamaño de malla del glicocálix cerebral envejecido, no "
     "solo el espesor."),
    ("Cuantitativa", "alta", "Sustituir el booleano de cada compuerta por una "
     "probabilidad con banda de incertidumbre. Es lo que elimina el "
     "«no evaluable» sin inventarse ningún dato."),
]

FUENTES = [
    ("Weinbaum S, Zhang X, Han Y, Vink H, Cowin SC (2003)",
     "Mechanotransduction and flow across the endothelial glycocalyx",
     "PNAS 100(13):7988-7995", "https://doi.org/10.1073/pnas.1332808100",
     "Tamiz del glicocálix"),
    ("Kabedev A, Lobaskin V (2022)",
     "Potential of mean force of a nanoparticle in a fibrous glycocalyx layer",
     "Nanomedicine 17:979-996", "https://doi.org/10.2217/nnm-2021-0387",
     "kT_hombro() / g_glicocalix_pmf"),
    ("Tadros T (2007)",
     "General Principles of Colloid Stability and the Role of Surface Forces, "
     "cap.1 de Colloid Stability: The Role of Surface Forces, Part I",
     "Wiley-VCH, ISBN 978-3-527-31462-1", "",
     "Umbral Gmax>25kT (por analogía, no dato de Kabedev)"),
    ("Verwey EJW, Overbeek JTG (1948)",
     "Theory of Stability of Lyophobic Colloids",
     "Elsevier", "",
     "Origen del umbral Gmax>25kT citado por Tadros 2007"),
    ("Lockman PR, Koziara JM, Mumper RJ, Allen DD (2004)",
     "Nanoparticle surface charges alter blood-brain barrier integrity and permeability",
     "J Drug Target 12(9-10):635-641", "https://doi.org/10.1080/10611860400015936",
     "g_glicocalix_tamiz: contraste real, BHE nativa in situ"),
    ("Gromnicova R, Kaya M, Romero IA, Williams P, Satchell S, Sharrack B, Male D (2016)",
     "Transport of Gold Nanoparticles by Vascular Endothelium from Different Human Tissues",
     "PLoS ONE 11(8):e0161610", "https://doi.org/10.1371/journal.pone.0161610",
     "g_glicocalix_tamiz: glicocálix cerebral no gatilla captación en hCMEC/D3"),
    ("Wang X, Shen B, Yang W, Wang X, Li C, Wu H (2026)",
     "A physics-informed neural network framework for quantitative analysis of "
     "transcytosis and physical diffusion in an in vitro BBB",
     "J Nanobiotechnology 24:164", "https://doi.org/10.1186/s12951-026-04023-y",
     "g_caveola: cruce vía clatrina/dinamina fuera del rango de caveola"),
    ("Deserno M (2004)",
     "Elastic deformation of a fluid membrane upon colloid binding",
     "Phys Rev E 69:031903", "https://doi.org/10.1103/PhysRevE.69.031903",
     "Envolvimiento de membrana"),
    ("Pan J, Tristram-Nagle S, Kučerka N, Nagle JF (2008)",
     "Temperature dependence of structure, bending rigidity, and bilayer "
     "interactions of DOPC bilayers",
     "Phys Rev Lett 100:198103", "https://doi.org/10.1103/PhysRevLett.100.198103",
     "Suelo geométrico del liposoma"),
    ("Nance EA et al. (2012)",
     "A dense poly(ethylene glycol) coating improves penetration of large "
     "polymeric nanoparticles within brain tissue",
     "Sci Transl Med 4(149):149ra119",
     "https://doi.org/10.1126/scitranslmed.3003594",
     "Difusión en el espacio extracelular"),
    ("Maiti PK, Çağin T, Wang G, Goddard WA III (2004)",
     "Structure of PAMAM Dendrimers: Generations 1 through 11",
     "Macromolecules 37(16):6236-6254", "https://doi.org/10.1021/ma035629b",
     "Techo del dendrímero"),
    ("Prosa TJ, Bauer BJ, Amis EJ (2001)",
     "From Stars to Spheres: A SAXS Analysis of Dilute Dendrimer Solutions",
     "Macromolecules 34(14):4897-4906", "https://doi.org/10.1021/ma0002186",
     "Tamaño medido del dendrímero"),
    ("de Gennes PG, Hervet H (1983)",
     "Statistics of «starburst» polymers",
     "J Physique Lettres 44(9):351-360",
     "https://doi.org/10.1051/jphyslet:01983004409035100",
     "Generación límite, teoría"),
    ("Parker NG, Mather ML, Morgan SP, Povey MJW (2010)",
     "Longitudinal acoustic properties of poly(lactic acid) and "
     "poly(lactic-co-glycolic acid)",
     "Biomed Mater 5:055004", "https://doi.org/10.1088/1748-6041/5/5/055004",
     "Densidad de PLA y PLGA"),
    ("Devarakonda B, Hill RA, de Villiers MM (2004)",
     "The effect of PAMAM dendrimer generation size and surface functional "
     "group on the aqueous solubility of nifedipine",
     "Int J Pharm 284(1-2):133-140",
     "https://doi.org/10.1016/j.ijpharm.2004.07.006",
     "Alojamiento del fármaco y carga 1:1"),
    ("Shi SM et al. (2025)",
     "Glycocalyx dysregulation impairs blood-brain barrier in ageing and disease",
     "Nature 639(8056):985-994", "https://doi.org/10.1038/s41586-025-08589-9",
     "Glicocálix envejecido"),
    ("Tracy GC et al. (2023)",
     "Intracerebral Nanoparticle Transport Facilitated by Alzheimer Pathology "
     "and Age", "Nano Lett 23(23):10971-10982",
     "https://doi.org/10.1021/acs.nanolett.3c03222", "Validación, caso C1"),
    ("Cheng MJ, Kumar R, Sridhar S, Webster TJ, Ebong EE (2016)",
     "Endothelial glycocalyx conditions influence nanoparticle uptake for "
     "passive targeting", "Int J Nanomedicine 11:3305-3315",
     "https://doi.org/10.2147/IJN.S106299", "Compuerta C.2"),
    ("Cheng MJ et al. (2019)",
     "Ultrasmall gold nanorods: synthesis and glycocalyx-related permeability "
     "in human endothelial cells", "Int J Nanomedicine 14:319-333",
     "https://doi.org/10.2147/IJN.S184455", "Compuerta C.2, bajo flujo"),
    ("González-Carter D et al. (2020)",
     "Targeting nanoparticles to the brain by exploiting the blood-brain "
     "barrier impermeability to selectively label the brain endothelium",
     "PNAS 117(32):19141-19150", "https://doi.org/10.1073/pnas.2002016117",
     "Señal endotelial frente a parenquimatosa"),
    ("Cortés H et al. (2020)",
     "A Reevaluation of Chitosan-Decorated Nanoparticles to Cross the "
     "Blood-Brain Barrier", "Membranes 10(9):212",
     "https://doi.org/10.3390/membranes10090212", "Ruta C, revisado en V-6"),
    ("Mao Y et al. (2014)",
     "A novel liposomal formulation of FTY720 (fingolimod) for promising "
     "enhanced targeted delivery", "Nanomedicine 10(2):393-400",
     "https://doi.org/10.1016/j.nano.2013.08.001",
     "Potencial ζ del liposoma con fingolimod · cota inferior de descarga"),
    ("Tong H-I, Kang W, Davy PMC, Shi Y, Sun S, Allsopp RC, Lu Y (2016)",
     "Monocyte Trafficking, Engraftment, and Delivery of Nanoparticles and an "
     "Exogenous Gene into the Acutely Inflamed Brain Tissue",
     "PLoS ONE 11(4):e0154022",
     "https://doi.org/10.1371/journal.pone.0154022",
     "Tránsito del monocito al cerebro inflamado, compuerta B.3"),
    ("Berry S, Mastorakos P, Zhang C, Song E, Patel H, Suk JS, Hanes J (2016)",
     "Enhancing intracranial delivery of clinically relevant non-viral gene "
     "vectors",
     "RSC Advances 6:41665–41674",
     "https://doi.org/10.1039/c6ra01546h",
     "Único ζ positivo bajo medido en parénquima: +10.0 mV, &lt;10 % difunde"),
    ("Mastorakos P, Song E, Zhang C, Berry S, Park HW, Kim YE, Park JS, "
     "Lee S, Suk JS, Hanes J (2016)",
     "Biodegradable DNA nanoparticles that provide widespread gene delivery "
     "in the brain",
     "Small 12(5):678–685",
     "https://doi.org/10.1002/smll.201502554",
     "ζ +35.3 mV inmovilizado en el parénquima, por MPT ex vivo"),
    ("Gong X, Fan X, He Y, Wang Y, Zhou F, Yang B (2022)",
     "A pH-sensitive liposomal co-delivery of fingolimod and ammonia borane "
     "for treatment of intracerebral hemorrhage",
     "Nanophotonics 11(22):5133–5142",
     "https://doi.org/10.1515/nanoph-2022-0496",
     "Diseño real del catálogo · 145 nm, ζ −28.33 mV"),
    ("Chow SF et al. (2025)",
     "Rational development of fingolimod nano-embedded microparticles as "
     "nose-to-brain neuroprotective therapy for ischemic stroke",
     "Drug Deliv Transl Res 15(6):2022–2047",
     "https://doi.org/10.1007/s13346-024-01721-8",
     "Diseño real del catálogo · 134 nm, ζ −0.24 mV"),
    ("Muselman A, Yu LW, Nguyen KD, Inayathullah M, Liu Q, Brewer KD, "
     "Malkovskiy AV, Rajadas J, Engleman EG (2026)",
     "Macrophage-targeted PEGylated liposomes ameliorate experimental "
     "autoimmune encephalomyelitis",
     "Front Immunol 16:1657131",
     "https://doi.org/10.3389/fimmu.2025.1657131",
     "Diseño real del catálogo · 700 nm, ζ 0.00 mV · autores y DOI "
     "verificados 2026-08-17 (C10), leído entero vía PMC (PMC12852013)"),
    ("Sochor B, Düdükcü Ö, Lübtow MM, Schummer B, Jaksch S, "
     "Luxenhofer R (2020)",
     "Probing the complex loading-dependent structural changes in "
     "ultrahigh drug-loaded polymer micelles by small-angle neutron scattering",
     "Langmuir 36(13):3494–3503",
     "https://doi.org/10.1021/acs.langmuir.9b03460",
     "Suelo de la micela cargada, 13.0 nm por SANS · tareas G.4 y G.5"),
    ("Israelachvili JN, Mitchell DJ, Ninham BW (1976)",
     "Theory of self-assembly of hydrocarbon amphiphiles into micelles and "
     "bilayers",
     "J Chem Soc Faraday Trans 2, 72:1525–1568",
     "https://doi.org/10.1039/F29767201525",
     "Parámetro de empaquetamiento · arquitectura de la micela"),
    ("Bastiani M, Parton RG (2010)",
     "Caveolae at a glance",
     "J Cell Sci 123(22):3831–3836",
     "https://doi.org/10.1242/jcs.070102",
     "Compuerta de caveola · 60–80 nm, pág. 3831"),
    ("Morris CE, Homann U (2001) · CITA DE SEGUNDA MANO",
     "Cell surface area regulation and membrane tension",
     "J Membr Biol 179(2):79–102 · el proyecto NO lo ha leído directo: entra a "
     "través de Shi y Baumgart 2015",
     "https://doi.org/10.1007/s002320010040",
     "Rango canónico de tensión lateral, 0.003–0.3 mN/m · barrido SIGMA_mNm"),
    ("Shi Z, Baumgart T (2015)",
     "Membrane tension and peripheral protein density mediate membrane shape "
     "transitions",
     "Nat Commun 6:5974",
     "https://doi.org/10.1038/ncomms6974",
     "Vía por la que entra el rango de tensión lateral de Morris y Homann"),
    ("Walter FR, Santa-Maria AR, Mészáros M, Veszelka S, Dér A, Deli MA "
     "(2021)",
     "Surface charge, glycocalyx, and blood-brain barrier function",
     "Tissue Barriers 9(3):1904773",
     "https://doi.org/10.1080/21688370.2021.1904773",
     "Grosor del glicocálix como contraste · densidad en microvaso cerebral"),
    ("Kincses A, Santa-Maria AR, Walter FR, Dér L, Horányi N, Lipka DV, "
     "Valkai S, Deli MA, Dér A (2020)",
     "A chip device to determine surface charge properties of confluent cell "
     "monolayers by measuring streaming potential",
     "Lab Chip 20(20):3792–3805",
     "https://doi.org/10.1039/D0LC00558D",
     "Carga superficial de monocapas endoteliales"),
    ("Santa-Maria AR, Walter FR, Figueiredo R, Kincses A, Vigh JP, "
     "Heymans M, Culot M, Winter P, Gosselet F, Dér A, Deli MA (2019)",
     "Lidocaine turns the surface charge of biological membranes more "
     "positive and changes the permeability of blood-brain barrier culture "
     "models",
     "Biochim Biophys Acta Biomembr 1861(9):1579–1591",
     "https://doi.org/10.1016/j.bbamem.2019.07.008",
     "Potencial ζ de la superficie endotelial de la BHE"),
    ("Yona S, Kim K-W, Wolf Y, Mildner A, Varol D, Breker M, "
     "Strauss-Ayali D, Viukov S, Guilliams M, Misharin A, Hume DA, "
     "Perlman H, Malissen B, Zelzer E, Jung S (2013)",
     "Fate mapping reveals origins and dynamics of monocytes and tissue "
     "macrophages under homeostasis",
     "Immunity 38(1):79–91 · RETIRADA del veredicto el 2026-08-12",
     "https://doi.org/10.1016/j.immuni.2012.12.001",
     "Semivida del monocito · se conserva solo como registro histórico"),
    ("Nagle JF (2017)",
     "Experimentally determined tilt and bending moduli of single-component "
     "lipid bilayers",
     "Chem Phys Lipids 205:18–24",
     "https://doi.org/10.1016/j.chemphyslip.2017.04.006",
     "Anclajes del barrido de κ · K_C de diez bicapas de un solo lípido"),
    ("Campbell SD, Regina KJ, Kharasch ED (2014)",
     "Significance of lipid composition in a blood-brain barrier-mimetic "
     "PAMPA assay",
     "J Biomol Screen 19(3):437–444",
     "https://doi.org/10.1177/1087057113497981",
     "Composición lipídica del endotelio cerebral humano · tarea C5"),
    ("Shi SM, Suh RJ, Shon DJ, Garcia FJ, Buff JK, Atkins M, Li L, Lu N, "
     "Sun B, Luo J, To N-S, Cheung TH, McNerney MW, Heiman M, Bertozzi CR, "
     "Wyss-Coray T (2025)",
     "Glycocalyx dysregulation impairs blood–brain barrier in ageing and "
     "disease",
     "Nature 639:985–994",
     "https://doi.org/10.1038/s41586-025-08589-9",
     "Espesor del glicocálix en capilar CEREBRAL: 540 nm joven, 232 nm viejo"),
    ("Larsen R, Kucharz K, Aydin S, Micael MKB, Choudhury B, "
     "Paulchakrabarti M, Lønstrup M, Lin DC, Abeln M, Münster-Kühnel A, "
     "Gomez Toledo A, Lauritzen M, Esko JD, Daneman R (2025) · PREPRINT",
     "Multi-omic analysis reveals the unique glycan landscape of the "
     "blood-brain barrier glycocalyx",
     "bioRxiv 2025.04.07.645297 · SIN revisión por pares",
     "https://doi.org/10.1101/2025.04.07.645297",
     "726 nm de espesor y 93 % de cobertura · el glicocálix NO se degrada "
     "en EAE"),
    ("Thorne RG, Nicholson C (2006)",
     "In vivo diffusion analysis with quantum dots and dextrans predicts the "
     "width of brain extracellular space",
     "PNAS 103(14):5567–5572",
     "https://doi.org/10.1073/pnas.0509425103",
     "Escenario conservador de 38 nm · umbral del fármaco liberado"),
    ("Hisano Y, Kobayashi N, Kawahara A, Yamaguchi A, Nishi T (2011)",
     "The sphingosine 1-phosphate transporter, SPNS2, functions as a "
     "transporter of the phosphorylated form of the immunomodulating agent "
     "FTY720",
     "J Biol Chem 286(3):1758–1766",
     "https://doi.org/10.1074/jbc.M110.171116",
     "Salida del fármaco de la célula, compuerta B.5 · FTY720-fosfato"),
    ("Foster CA, Howard LM, Schweitzer A, Persohn E, Hiestand PC, Balatoni B, "
     "Reuschel R, Beerli C, Schwartz M, Billich A (2007)",
     "Brain penetration of the oral immunomodulatory drug FTY720 and its "
     "phosphorylation in the central nervous system during experimental "
     "autoimmune encephalomyelitis: consequences for mode of action in "
     "multiple sclerosis",
     "J Pharmacol Exp Ther 323(2):469–475",
     "https://doi.org/10.1124/jpet.107.127183",
     "Partición del FTY720-P: 30–80x menos en LCR que en tejido (Tabla 3) · "
     "dependencia de portador declarada por los autores · compuerta de "
     "difusión del fármaco liberado"),
    ("Bucki R, Kulakowska A, Byfield FJ, Zendzian-Piotrowska M, Baranowski M, "
     "Marzec M, Winer JP, Ciccarelli NJ, Górski J, Drozdowski W, Bittman R, "
     "Janmey PA (2010)",
     "Plasma gelsolin modulates cellular response to sphingosine 1-phosphate",
     "Am J Physiol Cell Physiol 299(6):C1516–C1523",
     "https://doi.org/10.1152/ajpcell.00051.2010",
     "La gelsolina une el FTY720-P de forma débil o nula, a diferencia del "
     "S1P (Fig. 1 A-E) · descarta ese portador candidato en LCR"),
    ("Mouzoura P et al. (2025)",
     "Formulación liposomal de FTY720 (fingolimod): relación molar "
     "fármaco:lípido 1:8 y eficiencia de carga 94–97.2 %",
     "Int J Nanomedicine 20:239–265",
     "https://doi.org/10.2147/IJN.S494512",
     "Carga útil del liposoma · compuerta G.2 · FTY720 directo, leído entero"),
    ("Mishima Y, Kurano M, Kobayashi T, Nishikawa M, Ohkawa R, Tozuka M, "
     "Yatomi Y (2018)",
     "Dihydro-sphingosine 1-phosphate interacts with carrier proteins in a "
     "manner distinct from that of sphingosine 1-phosphate",
     "Biosci Rep 38(5):BSR20181288",
     "https://doi.org/10.1042/BSR20181288",
     "Un análogo del S1P NO hereda sus portadores: el S1P se une a HDL vía "
     "apoM y el DH-S1P no · impide asumir la biología de portador del S1P "
     "para el FTY720-P"),
    ("Curtis C, McKenna M, Pontes C, Toghani D, Choe A, Nance E (2019)",
     "Predicting in situ nanoparticle behavior using multiple particle "
     "tracking and artificial neural networks",
     "Nanoscale 11(46):22515–22530",
     "https://doi.org/10.1039/c9nr06327g",
     "Único punto dentro de la banda 114–200 nm: PS-PEG 163.2 nm, ζ −6.2 mV"),
    ("McKenna M, Shackelford D, Pontes C, Ball B, Nance E (2021)",
     "Multiple particle tracking detects changes in brain extracellular "
     "matrix and predicts neurodevelopmental age",
     "ACS Nano 15(5):8559–8573",
     "https://doi.org/10.1021/acsnano.1c00394",
     "Dependencia con la EDAD: razón 5 (P14) a 34 (P70), poro 76.8 a 36.0 nm"),
]

# -----------------------------------------------------------------------------
#  FIGURAS DE LA WEB
#
#  Cada página de clase lleva DOS bloques:
#    1. SET COMÚN — las mismas tres figuras en todas: ventanas, matriz y
#       recorrido, generadas por `rutas.figuras()` con el catálogo de esa
#       página. Mismo orden y mismo pie siempre.
#    2. ESPECÍFICAS DE ESTA CLASE — lo que solo existe para una clase porque
#       solo para ella hay dato. No se replica en las demás: hacerlo obligaría
#       a inventar números.
# -----------------------------------------------------------------------------

# Física del envolvimiento: NO depende de ningún diseño concreto. Son barridos
# de κ y de σ̃, así que pueden ir en la página de la clase.
FIGURAS_ENVOLVIMIENTO = [
    ("envolvimiento_radio_critico.png", "Radio crítico de envolvimiento",
     "Frente a rigidez de membrana y adhesión."),
    ("envolvimiento_barrera.png", "Barrera de energía", ""),
]

# G(D) SÍ depende de diseños concretos: su leyenda son los TRES TEÓRICOS
# (convencional, furtivo/PEG, catiónico). Estaba en la página de liposomas
# «Reales», que por acuerdo solo lleva lo medido, y ver ahí nombres de diseños
# teóricos hacía que las dos páginas se leyeran contradictorias. Va a la
# sub-página de teóricos.
FIGURAS_ENVOLVIMIENTO_TEORICOS = [
    ("envolvimiento_G_de_D.png", "Energía libre G(D) · los tres diseños teóricos",
     "Frente a la separación partícula-membrana. Pozos: convencional −21.5 kT, "
     "furtivo/PEG −7.5 kT, catiónico −19.2 kT. La barrera de entrada solo "
     "aparece con PEG."),
]

def _figuras_catalogo(prefijo):
    """EL SET COMÚN. Las mismas tres figuras, en el mismo orden y con el mismo
    pie, en TODAS las páginas de clase. Solo cambia el catálogo que dibujan."""
    return [
        (f"{prefijo}_ventanas.png",
         "Ventanas de tamaño · solo las compuertas que dependen del diámetro",
         "Barra = rango permitido. Línea vertical = diseño evaluado. NO es un "
         "veredicto: las compuertas que no dependen del tamaño no salen aquí, "
         "así que una línea en verde en todas las barras de una ruta no significa "
         "que la supere. El veredicto está en la matriz."),
        (f"{prefijo}_matriz.png", "Diseño por ruta",
         "! = compuerta que pasa con salvedad."),
        (f"{prefijo}_recorrido.png", "Compuerta a compuerta",
         "Verde pasa · rojo falla · gris sin dato · hueco = salvedad."),
    ]

FIGURAS_LIP_REALES = _figuras_catalogo("lip_reales")
FIGURAS_LIP_TEORICOS = _figuras_catalogo("lip_teoricos")
FIGURAS_DEND_REALES = _figuras_catalogo("dend_reales")
FIGURAS_DEND_TEORICOS = _figuras_catalogo("dend_teoricos")
FIGURAS_POL_TEORICOS = _figuras_catalogo("pol_teoricos")

FIGURAS_POLIMERO_TEORICOS_PROPIAS = [
    ("polimeros_teoricos.png", "Los tres contra las dos ventanas",
     "Con el suelo de su propia clase: glóbulo colapsado o micela cargada."),
]

FIGURAS_POLIMERO = [
    ("polimero_suelo.png", "Suelo frente a masa molar",
     "Una curva por cada densidad medida por Parker 2010."),
    ("clases_ventanas.png", "Las tres clases contra las dos ventanas",
     "La franja rayada no pertenece a ninguna ventana."),
]

FIGURAS_DENDRIMERO_PROPIAS = [
    ("dendrimero_ventana.png", "Ventana geométrica y zoom del techo", ""),
    ("dendrimero_generaciones.png", "Diámetro por generación",
     "SAXS medido frente a dinámica molecular calculada."),
    ("dendrimero_vs_liposoma.png", "Las dos clases contra las dos ventanas", ""),
    ("dendrimero_carga.png", "K(1:1) por generación y pH",
     "Devarakonda 2004, Tabla 2."),
]

FIGURAS_DENDRIMERO_TEORICOS_PROPIAS = [
    ("dendrimeros_teoricos.png", "Los tres sobre la ventana geométrica",
     "Verde pasa · rojo falla · gris sin dato para esa química."),
    ("dend_teoricos_vs_medidos.png", "Fichas frente a las generaciones medidas",
     "El Ø de la ficha del PAMAM G4 cae por encima del G10, no donde el G4."),
    ("dend_teoricos_pmf.png", "El modelo de acople de las fichas",
     "El pozo depende solo del logP del fármaco; la barrera de superficie, del ζ."),
]

# Todas las que hay que copiar a web/img.
TODAS_LAS_FIGURAS = (FIGURAS_ENVOLVIMIENTO + FIGURAS_ENVOLVIMIENTO_TEORICOS
                     + FIGURAS_LIP_REALES + FIGURAS_LIP_TEORICOS
                     + FIGURAS_DEND_REALES + FIGURAS_DEND_TEORICOS
                     + FIGURAS_POL_TEORICOS
                     + FIGURAS_DENDRIMERO_PROPIAS
                     + FIGURAS_DENDRIMERO_TEORICOS_PROPIAS
                     + FIGURAS_POLIMERO_TEORICOS_PROPIAS
                     + FIGURAS_POLIMERO)

FICHAS = [
    ("verificacion_dendrimero_tarea_G_1a.md", "G.1a · Límite geométrico del dendrímero"),
    ("verificacion_polimero_micela_tarea_G_1b.md", "G.1b · Límite del polímero macizo (incluye la micela, ya fuera de alcance)"),
    ("verificacion_compuertas_C2_D3.md", "C.2 y D.3 · Unión al glicocálix y acceso al receptor"),
    ("verificacion_endotelio_tarea_V_6.md", "V-6 · Señal endotelial frente a parenquimatosa"),
    ("verificacion_envejecimiento_tarea_V_1B.md", "V-1B · El caso del envejecimiento"),
    ("verificacion_nance_tarea_V_1.md", "V-1 · Difusión en el espacio extracelular"),
    ("verificacion_transito_tarea_B_3.md", "B.3 · Tránsito del monocito"),
    ("verificacion_zeta_positivo_difusion.md", "ζ positivo · Difusión extracelular"),
    ("verificacion_zona_gris_tamano.md", "Zona gris 114–200 nm · Difusión extracelular"),
    ("verificacion_fty720_fosfato.md", "F.1 · FTY720-fosfato vs fingolimod neutro"),
    ("verificacion_kappa_C5_C6.md", "C5 y C6 · el barrido de κ"),
    ("resultados_validacion.md", "Resultados de la validación"),
    ("protocolo_validacion.md", "Protocolo de validación"),
]


# =============================================================================
#  MARKDOWN MÍNIMO  (solo lo que usan las fichas)
# =============================================================================

def _en_linea(t):
    t = html.escape(t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    # No-avaricioso y admitiendo '*' dentro: las fichas escriben cosas como
    # "**3.4. Es *in vitro*.**", y con [^*]+ la negrita no se emparejaba y los
    # asteriscos salían impresos. La cursiva de dentro la resuelve la regla
    # siguiente, que corre después.
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", t)
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', t)
    t = re.sub(r"(?<![\"=>])(https?://\S+)", r'<a href="\1">\1</a>', t)
    return t


def markdown(texto):
    """Convertidor mínimo: títulos, tablas, listas, citas, reglas y párrafos."""
    salida, buffer_p, en_tabla = [], [], False

    def cerrar_p():
        if buffer_p:
            # Se unen con espacio, NO con <br>: las fichas están escritas a 80
            # columnas y con <br> la web heredaba los cortes del archivo fuente.
            #
            # OJO AL ORDEN: primero se JUNTA y después se formatea. Al revés
            # (formatear línea a línea y luego juntar) una negrita que abre en
            # una línea y cierra en la siguiente no se empareja nunca y los
            # asteriscos salen impresos. Las 17 fichas del proyecto tienen ese
            # caso, porque están escritas a 80 columnas.
            salida.append("<p>" + _en_linea(" ".join(buffer_p)) + "</p>")
            buffer_p.clear()

    def cerrar_tabla():
        nonlocal en_tabla
        if en_tabla:
            salida.append("</tbody></table></div>")
            en_tabla = False

    lineas = texto.split("\n")
    i = 0
    while i < len(lineas):
        ln = lineas[i].rstrip()
        if not ln.strip():
            cerrar_p(); cerrar_tabla(); i += 1; continue
        if re.match(r"^---+$", ln.strip()):
            cerrar_p(); cerrar_tabla(); salida.append("<hr>"); i += 1; continue
        m = re.match(r"^(#{1,6})\s+(.*)$", ln)
        if m:
            cerrar_p(); cerrar_tabla()
            n = len(m.group(1)) + 1
            salida.append(f"<h{min(n,6)}>{_en_linea(m.group(2))}</h{min(n,6)}>")
            i += 1; continue
        # tabla
        if ln.lstrip().startswith("|") and i + 1 < len(lineas) \
                and re.match(r"^\s*\|[\s:|-]+\|\s*$", lineas[i + 1]):
            cerrar_p()
            cabecera = [c.strip() for c in ln.strip().strip("|").split("|")]
            salida.append('<div class="tabla-scroll"><table><thead><tr>'
                          + "".join(f"<th>{_en_linea(c)}</th>" for c in cabecera)
                          + "</tr></thead><tbody>")
            en_tabla = True
            i += 2
            while i < len(lineas) and lineas[i].lstrip().startswith("|"):
                cel = [c.strip() for c in lineas[i].strip().strip("|").split("|")]
                salida.append("<tr>" + "".join(f"<td>{_en_linea(c)}</td>" for c in cel) + "</tr>")
                i += 1
            cerrar_tabla(); continue
        if re.match(r"^\s*[-*]\s+", ln) or re.match(r"^\s*\d+\.\s+", ln):
            cerrar_p(); cerrar_tabla()
            ordenada = bool(re.match(r"^\s*\d+\.\s+", ln))
            # Mismo criterio que en cerrar_p: se acumula el TEXTO CRUDO de cada
            # ítem y se formatea una sola vez al cerrarlo. Un ítem partido en
            # varias líneas con una negrita a caballo se emparejaba mal.
            items = []
            while i < len(lineas) and (re.match(r"^\s*[-*]\s+", lineas[i])
                                       or re.match(r"^\s*\d+\.\s+", lineas[i])
                                       # continuación del ítem: 2 espacios o
                                       # más. Exigir 3 dejaba fuera la sangría
                                       # de 2 que usan casi todas las fichas.
                                       or re.match(r"^\s{2,}\S", lineas[i])):
                if re.match(r"^\s*([-*]|\d+\.)\s+", lineas[i]):
                    items.append(re.sub(r"^\s*([-*]|\d+\.)\s+", "", lineas[i]).strip())
                elif items:
                    items[-1] += " " + lineas[i].strip()
                i += 1
            salida.append("<ol>" if ordenada else "<ul>")
            salida.extend(f"<li>{_en_linea(x)}</li>" for x in items)
            salida.append("</ol>" if ordenada else "</ul>")
            continue
        if ln.lstrip().startswith(">"):
            cerrar_p(); cerrar_tabla()
            # Una cita de varias líneas es UNA cita, no una por línea, y se
            # formatea entera para no partir las negritas.
            cita = []
            while i < len(lineas) and lineas[i].lstrip().startswith(">"):
                cita.append(lineas[i].lstrip().lstrip(">").strip())
                i += 1
            salida.append("<blockquote>" + _en_linea(" ".join(cita).strip())
                          + "</blockquote>")
            continue
        if ln.strip().startswith("<"):
            cerrar_p(); i += 1; continue
        buffer_p.append(ln.strip())
        i += 1
    cerrar_p(); cerrar_tabla()
    return "\n".join(salida)


# =============================================================================
#  DATOS  —  todo sale de ejecutar el simulador, nada está escrito a mano
# =============================================================================

def recoger():
    u = dict(
        glicocalix=2.0 * G.radio_exclusion_nm(),
        envolvimiento=R.g_envolvimiento(R.Diseno("_", 20.0, 0.0)).umbral,
        fagocitosis=550.0,
        liposoma_min=G.diametro_liposoma_minimo_nm(4.0, 4.0),
    )
    d = dict(umbrales=u, rutas=list(R.RUTAS.keys()))

    d["margen_ventanas"] = u["envolvimiento"] - u["glicocalix"]

    def _filas(catalogo):
        filas = []
        for dis in catalogo:
            ev = R.evaluar(dis)
            fila = dict(nombre=dis.nombre, diametro=dis.diametro_nm,
                        zeta=dis.zeta_mV, peg=dis.peg_nm, nota=dis.nota, rutas={})
            for nombre, (ver, res) in ev.items():
                salvedades = (sum(1 for x in res if x.advertencia)
                              if ver == "NO EXCLUIDA" else 0)
                fila["rutas"][nombre] = dict(
                    veredicto=ver, salvedades=salvedades,
                    muere=R._quien_lo_mata(res),
                    faltan=[x.compuerta for x in res if x.estado == R.DESCONOCIDA])
            filas.append(fila)
        return filas

    d["catalogo"] = _filas(R.CATALOGO)
    d["catalogo_real"] = _filas(R.CATALOGO_REAL)
    d["catalogo_teorico"] = _filas(R.CATALOGO_TEORICO)

    # -------------------------------------------------------------------
    #  DATASET DE 50 LIPOSOMAS SINTÉTICOS (2026-08-18). Semilla fija en
    #  dataset_50_liposomas.py: mismos 50 diseños en cada build. Reutiliza
    #  _filas() para la tabla resumen (mismo formato que catalogo_real/
    #  teorico) y además guarda el desglose completo de compuertas + genera
    #  las 3 figuras propias de CADA liposoma para el detalle desplegable.
    # -------------------------------------------------------------------
    import contextlib as _ctxlib
    import io as _io

    _disenos_50 = _generar_dataset_50()
    d["dataset_50"] = _filas(_disenos_50)

    _img_dir_50 = SALIDA / "img" / "dataset_50"
    _img_dir_50.mkdir(parents=True, exist_ok=True)
    _detalle_50 = []
    for _i, _dis in enumerate(_disenos_50, start=1):
        _ev = R.evaluar(_dis)
        _rutas_det = {}
        for _nombre, (_ver, _res) in _ev.items():
            _rutas_det[_nombre] = dict(
                veredicto=_ver,
                compuertas=[dict(nombre=_r.compuerta, estado=_r.estado,
                                 valor=_r.valor, umbral=_r.umbral,
                                 unidad=_r.unidad, fuente=_r.fuente)
                           for _r in _res])
        _prefijo = str(_img_dir_50 / f"liposoma_{_i:02d}")
        with _ctxlib.redirect_stdout(_io.StringIO()):
            R.figuras(prefijo=_prefijo, catalogo=[_dis], incluir_ventanas=True)
        _detalle_50.append(dict(
            indice=_i, nombre=_dis.nombre, diametro=_dis.diametro_nm,
            zeta=_dis.zeta_mV, peg=_dis.peg_nm,
            img=f"img/dataset_50/liposoma_{_i:02d}", rutas=_rutas_det))
    d["dataset_50_detalle"] = _detalle_50

    # -------------------------------------------------------------------
    #  DETALLE POR DISEÑO — LIPOSOMA REALES Y TEÓRICOS (2026-08-18k). Mismo
    #  patrón que el dataset de 50: tarjeta desplegable con las 3 figuras
    #  propias de CADA diseño (no las genéricas del catálogo entero) + el
    #  desglose de compuertas por ruta + insignia "Fabricable" (primera
    #  compuerta de cada ruta, es la misma en las cuatro: g_transportador_
    #  fabricable). Unifica el formato con Dataset, como quedó descrito en
    #  el handoff del 2026-08-18 (sesión a) pero nunca se conectó al código.
    # -------------------------------------------------------------------
    def _detalle_por_diseno(catalogo, carpeta):
        _img_dir = SALIDA / "img" / carpeta
        _img_dir.mkdir(parents=True, exist_ok=True)
        _detalle = []
        for _i, _dis in enumerate(catalogo, start=1):
            _ev = R.evaluar(_dis)
            _rutas_det = {}
            _fabricable = None
            for _nombre, (_ver, _res) in _ev.items():
                _rutas_det[_nombre] = dict(
                    veredicto=_ver,
                    compuertas=[dict(nombre=_r.compuerta, estado=_r.estado,
                                     valor=_r.valor, umbral=_r.umbral,
                                     unidad=_r.unidad, fuente=_r.fuente)
                               for _r in _res])
                if _fabricable is None and _res:
                    _fabricable = _res[0].estado == R.PASA
            _prefijo = str(_img_dir / f"liposoma_{_i:02d}")
            with _ctxlib.redirect_stdout(_io.StringIO()):
                R.figuras(prefijo=_prefijo, catalogo=[_dis], incluir_ventanas=True)
            _detalle.append(dict(
                indice=_i, nombre=_dis.nombre, diametro=_dis.diametro_nm,
                zeta=_dis.zeta_mV, peg=_dis.peg_nm, nota=_dis.nota,
                fabricable=_fabricable,
                img=f"img/{carpeta}/liposoma_{_i:02d}", rutas=_rutas_det))
        return _detalle

    d["catalogo_real_detalle"] = _detalle_por_diseno(R.CATALOGO_REAL, "liposoma_reales")
    d["catalogo_teorico_detalle"] = _detalle_por_diseno(R.CATALOGO_TEORICO,
                                                         "liposoma_teoricos_detalle")

    d["liposoma"] = dict(
        t_bicapa=list(G.T_BICAPA_nm),
        min_con_nucleo=[G.diametro_liposoma_minimo_nm(4.0, t) for t in G.T_BICAPA_nm],
        min_nucleo_nulo=[G.diametro_liposoma_minimo_nm(0.0, t) for t in G.T_BICAPA_nm],
        limite=G.puede_existir_liposoma_que_pase(),
    )

    d["polimero"] = dict(
        densidades=dict(R.POLIMERO_DENSIDAD_g_cm3),
        suelos=[(etq, mw, R.POLIMERO_DENSIDAD_g_cm3[etq],
                 R.diametro_globulo_colapsado_nm(mw * 1000.0,
                                                 R.POLIMERO_DENSIDAD_g_cm3[etq]))
                for etq, mw in (("PLA15", 15), ("PLA24", 24),
                                ("PLGA 85:15", 53), ("PLA60", 60))],
        sensibilidad=[(n, R.diametro_globulo_colapsado_nm(53000.0 + n * 307.48, 1.19))
                      for n in (0, 1, 10)],
    )

    d["dendrimero"] = dict(
        suelo=R.DENDRIMERO_SUELO_nm, techo=R.DENDRIMERO_TECHO_nm,
        techo_medido=R.DENDRIMERO_TECHO_MEDIDO_nm, precision=R.DENDRIMERO_PRECISION,
        margen=u["envolvimiento"] - R.DENDRIMERO_TECHO_nm, generaciones=[])
    for g, dn in R.DENDRIMERO_GENERACIONES_nm.items():
        r = R.g_transportador_fabricable(R.Diseno(f"G{g}", dn, 0.0, clase="dendrimero"))
        d["dendrimero"]["generaciones"].append(dict(
            gen=g, diametro=dn, estado=r.estado,
            glicocalix=dn <= u["glicocalix"], envuelve=dn >= u["envolvimiento"]))

    # Transportadores teóricos del simulador de acople. DATOS SINTÉTICOS.
    def _teoricos(catalogo):
        filas = []
        for dis in catalogo:
            fila = dict(nombre=dis.nombre, diametro=dis.diametro_nm,
                        zeta=dis.zeta_mV, clase=dis.clase, nota=dis.nota,
                        puertas=[])
            for p in (R.g_transportador_fabricable, R.g_glicocalix_tamiz,
                      R.g_envolvimiento):
                r = p(dis)
                fila["puertas"].append(dict(nombre=r.compuerta, estado=r.estado,
                                            valor=r.valor, umbral=r.umbral,
                                            motivo=r.motivo,
                                            advertencia=r.advertencia))
            filas.append(fila)
        return filas

    d["teoricos"] = _teoricos(R.catalogo_dendrimeros_teoricos())
    d["teoricos_polimero"] = _teoricos(R.catalogo_polimeros_teoricos())
    d["plga_mw_kDa"] = R.PLGA_MW_DERIVADA_kDa
    d["micela_suelo_nm"] = R.MICELA_SUELO_CARGADA_nm
    d["micela_vacia_nm"] = R.MICELA_VACIA_nm

    # inventario de compuertas: se descubre recorriendo las rutas, no a mano
    vistas, comps = set(), []
    sonda = R.Diseno("sonda", 40.0, -2.0, 5.0)
    for nombre, cadena in R.RUTAS.items():
        for c in cadena:
            r = c(sonda)
            if r.compuerta in vistas:
                continue
            vistas.add(r.compuerta)
            comps.append(dict(nombre=r.compuerta, fuente=r.fuente or "—",
                              implementada=r.estado != R.DESCONOCIDA or bool(r.fuente),
                              motivo=r.motivo))
    d["compuertas"] = comps
    return d


# =============================================================================
#  HTML
# =============================================================================

CSS = """
:root{
  --tinta:#1c1f24; --suave:#5b6470; --linea:#e3e6ea; --fondo:#fbfcfd;
  --papel:#ffffff; --verde:#2e7d32; --rojo:#c62828; --gris:#78848f;
  --azul:#1565c0; --morado:#6a1b9a; --aviso:#b26a00; --avisofondo:#fff8ec;
  --ancho:min(1120px, 100% - 3rem);
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{
  margin:0; background:var(--fondo); color:var(--tinta);
  font:16px/1.65 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  -webkit-font-smoothing:antialiased; display:flex; flex-direction:column;
  min-height:100vh;
}
.envoltorio{width:var(--ancho); margin-inline:auto}

/* ---- cabeceras ---- */
header.principal{background:linear-gradient(160deg,#12213a,#1f3b63); color:#fff;
  padding:3.2rem 0 2.6rem}
header.pagina{background:linear-gradient(160deg,#12213a,#1f3b63); color:#fff;
  padding:2.1rem 0 1.9rem}
header .expediente{font-size:.75rem; letter-spacing:.16em; text-transform:uppercase;
  color:#9fc4f0; margin:0 0 .7rem}
header.principal h1{font-size:1.85rem; line-height:1.3; margin:0 0 .5rem; font-weight:650}
header.pagina h1{font-size:1.5rem; line-height:1.25; margin:0 0 .35rem; font-weight:650}
header .sub{color:#c9dcf3; margin:0 0 1.4rem; font-size:1rem}
header.pagina .sub{margin:0; font-size:.93rem}
.meta{display:flex; flex-wrap:wrap; gap:.5rem .6rem}
.meta span{background:rgba(255,255,255,.10); border:1px solid rgba(255,255,255,.16);
  padding:.3rem .7rem; border-radius:999px; font-size:.82rem}

/* ---- navegación ---- */
nav.barra{position:sticky; top:0; z-index:20; background:rgba(251,252,253,.94);
  backdrop-filter:blur(8px); border-bottom:1px solid var(--linea)}
nav.barra ul{display:flex; gap:.2rem; list-style:none; margin:0; padding:.5rem 0;
  overflow-x:auto; width:var(--ancho); margin-inline:auto}
nav.barra a{color:var(--suave); text-decoration:none; font-size:.87rem;
  white-space:nowrap; padding:.4rem .8rem; border-radius:7px; display:block}
nav.barra a:hover{background:#eef1f5; color:var(--tinta)}
nav.barra a.activo{background:#12213a; color:#fff; font-weight:600}

/* ---- sub-navegación de clase ---- */
nav.subbarra{background:#eef1f5; border-bottom:1px solid var(--linea)}
nav.subbarra .envoltorio{display:flex; gap:.35rem; padding:.4rem 0}
nav.subbarra a{color:var(--suave); text-decoration:none; font-size:.82rem;
  padding:.25rem .7rem; border-radius:999px; white-space:nowrap;
  border:1px solid transparent}
nav.subbarra a:hover{background:#fff; color:var(--tinta)}
nav.subbarra a.activo{background:#fff; color:var(--tinta); font-weight:600;
  border-color:var(--linea)}

main{flex:1; padding:2.4rem 0 1rem}
main > .envoltorio > h2:first-child{margin-top:0}
h2{font-size:1.3rem; margin:2.6rem 0 .3rem; font-weight:650;
  padding-bottom:.55rem; border-bottom:2px solid var(--linea)}
h3{font-size:1.05rem; margin:2rem 0 .6rem; font-weight:620}
.entradilla{color:var(--suave); margin:.7rem 0 1.5rem; max-width:64ch}

/* ---- tarjetas ---- */
.tarjetas{display:grid; gap:1rem; grid-template-columns:repeat(auto-fit,minmax(240px,1fr))}
.tarjeta{background:var(--papel); border:1px solid var(--linea); border-radius:11px;
  padding:1.15rem 1.25rem}
.tarjeta .cifra{font-size:1.7rem; font-weight:660; letter-spacing:-.02em; line-height:1.15}
.tarjeta .rotulo{font-size:.78rem; text-transform:uppercase; letter-spacing:.09em;
  color:var(--suave); margin-bottom:.45rem}
.tarjeta .pie{font-size:.83rem; color:var(--suave); margin-top:.45rem}
.tarjeta.verde .cifra{color:var(--verde)} .tarjeta.rojo .cifra{color:var(--rojo)}
.tarjeta.azul .cifra{color:var(--azul)}  .tarjeta.morado .cifra{color:var(--morado)}

/* ---- tarjetas de navegación de la portada ---- */
.indice{display:grid; gap:1rem; grid-template-columns:repeat(auto-fit,minmax(270px,1fr));
  margin-top:1.4rem}
a.bloque{display:block; text-decoration:none; color:inherit; background:var(--papel);
  border:1px solid var(--linea); border-radius:11px; padding:1.25rem 1.35rem;
  transition:border-color .15s, transform .15s}
a.bloque:hover{border-color:#9db6d4; transform:translateY(-2px)}
a.bloque .n{font-size:.74rem; letter-spacing:.1em; text-transform:uppercase;
  color:var(--azul); font-weight:700}
a.bloque b{display:block; font-size:1.08rem; margin:.3rem 0 .35rem}
a.bloque span{font-size:.88rem; color:var(--suave)}

/* ---- avisos ---- */
.aviso{background:var(--avisofondo); border:1px solid #f0d9ae;
  border-left:4px solid var(--aviso); border-radius:9px; padding:.95rem 1.15rem;
  margin:1.2rem 0; font-size:.92rem}
.aviso strong{color:var(--aviso)}

/* ---- tablas ---- */
.tabla-scroll{overflow-x:auto; margin:1.1rem 0; border:1px solid var(--linea);
  border-radius:10px; background:var(--papel)}
table{border-collapse:collapse; width:100%; font-size:.9rem}
th,td{padding:.6rem .85rem; text-align:left; border-bottom:1px solid var(--linea);
  vertical-align:top}
th{background:#f4f6f9; font-weight:600; font-size:.8rem; text-transform:uppercase;
  letter-spacing:.05em; color:var(--suave)}
tbody tr:last-child td{border-bottom:none}
td.num,th.num{text-align:right; font-variant-numeric:tabular-nums}

.etq{display:inline-block; padding:.16rem .55rem; border-radius:999px;
  font-size:.76rem; font-weight:600; white-space:nowrap}
.etq.no{background:#fdecea; color:var(--rojo)}
.etq.si{background:#e8f5e9; color:var(--verde)}
.etq.nn{background:#eceff1; color:var(--gris)}
.rev{font-size:.84rem; color:var(--suave)}

/* ---- figuras ---- */
.figura{background:var(--papel); border:1px solid var(--linea); border-radius:11px;
  overflow:hidden; margin:1.4rem 0}
.figura img{width:100%; display:block; background:#fff}
.figura figcaption{padding:.9rem 1.15rem; border-top:1px solid var(--linea)}
.figura figcaption b{display:block; margin-bottom:.2rem; font-size:.95rem}
.figura figcaption span{font-size:.87rem; color:var(--suave)}

/* ---- dataset de 50 liposomas ---- */
details.liposoma{border:1px solid var(--linea); border-radius:10px;
  margin-bottom:.6rem; padding:.5rem 1.1rem; background:var(--papel)}
details.liposoma summary{cursor:pointer; font-weight:600; padding:.4rem 0}
details.liposoma .dataset50-figs{display:grid; grid-template-columns:repeat(3,1fr);
  gap:1rem; margin:1rem 0}
details.liposoma .dataset50-figs .figura{margin:0}
@media (max-width:900px){details.liposoma .dataset50-figs{grid-template-columns:1fr}}

/* ---- listas ---- */
.fuente{padding:.85rem 0; border-bottom:1px solid var(--linea)}
.fuente:last-child{border-bottom:none}
.fuente .aut{font-weight:600; font-size:.92rem}
.fuente .tit{font-size:.92rem}
.fuente a{font-size:.83rem; color:var(--azul); word-break:break-all}
.fuente .uso{display:inline-block; font-size:.74rem; background:#eef2f7;
  color:var(--suave); padding:.13rem .5rem; border-radius:999px; margin-top:.25rem}
.prio{display:inline-block; font-size:.72rem; font-weight:700; padding:.14rem .5rem;
  border-radius:999px; text-transform:uppercase; letter-spacing:.05em}
.prio.alta{background:#fdecea; color:var(--rojo)}
.prio.media{background:#fff4e0; color:var(--aviso)}

code{background:#eef1f5; padding:.1rem .35rem; border-radius:4px; font-size:.87em;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
a{color:var(--azul)}
.formula{text-align:center; font-size:1.02rem; margin:1.1rem 0}

footer{margin-top:3rem; padding:1.8rem 0 2.6rem; border-top:1px solid var(--linea);
  color:var(--suave); font-size:.85rem; background:var(--papel)}

/* ---- fichas ---- */
.ficha{background:var(--papel); border:1px solid var(--linea); border-radius:11px;
  padding:1.8rem 2rem; margin:0 0 1.5rem}
.ficha h2{font-size:1.2rem; margin-top:1.8rem} .ficha h2:first-child{margin-top:0}
.ficha h3{font-size:1rem} .ficha h4{font-size:.95rem; margin:1.3rem 0 .4rem}
.ficha blockquote{border-left:3px solid var(--linea); margin:1rem 0;
  padding:.2rem 0 .2rem 1rem; color:var(--suave)}
.volver{display:inline-block; margin:0 0 1.3rem; font-size:.9rem}

@media (max-width:640px){
  header.principal h1{font-size:1.35rem} header.pagina h1{font-size:1.2rem}
  .ficha{padding:1.2rem}
}
@media print{nav.barra{display:none} body{background:#fff} .figura{break-inside:avoid}}
"""

# Las páginas del sitio. El orden es el de la navegación.
PAGINAS = [
    ("index.html", "Resumen", "Cifras y resultados."),
    ("liposoma.html", "Liposoma", "Suelo geométrico y liposomas reales publicados."),
    ("dendrimero.html", "Dendrímero", "Ventana G3-G10, carga y 6 figuras."),
    ("polimero.html", "Polímero", "Suelo del glóbulo colapsado."),
    ("metodo.html", "Método", "Veredictos, rutas y compuertas."),
    ("validacion.html", "Validación", "Recuento contra la literatura."),
    ("fichas.html", "Fichas", "Fuente primaria de cada resultado."),
    ("pendientes.html", "Pendientes", "Lo que falta."),
    ("fuentes.html", "Fuentes", "Bibliografía."),
]

# Sub-páginas de cada clase de transportador. La página de la clase lleva SOLO
# lo medido en fuente primaria; la sub-página, SOLO los diseños teóricos.
# clave = archivo de la clase; valor = [(archivo, etiqueta, subtítulo), ...]
SUBPAGINAS = {
    "liposoma.html": [
        ("liposoma.html", "Reales", ""),
        ("liposoma_teoricos.html", "Teóricos",
         "Formulaciones propuestas, sin ningún número medido."),
        ("dataset_50.html", "Dataset",
         "50 liposomas sintéticos, parámetros aleatorios dentro del rango real."),
    ],
    "dendrimero.html": [
        ("dendrimero.html", "Reales", ""),
        ("dendrimero_teoricos.html", "Teóricos",
         "Dendrímeros del simulador de acople. Datos sintéticos."),
    ],
    "polimero.html": [
        ("polimero.html", "Reales", ""),
        ("polimero_teoricos.html", "Teóricos",
         "Nanopartícula PLGA y dos micelas. Datos sintéticos."),
    ],
}

# archivo de sub-página -> archivo de su clase, para saber qué sub-barra pintar
# y qué pestaña principal marcar como activa.
_CLASE_DE = {sub: clase
             for clase, subs in SUBPAGINAS.items()
             for sub, _, _ in subs}

_NOMBRE_PAGINA = {a: n for a, n, _ in PAGINAS}

# Las sub-páginas también se generan, pero NO salen en la barra principal. El
# título lleva el nombre de la clase delante para que se sepa dónde estás.
SUBPAGINAS_EXTRA = [(a, f"{_NOMBRE_PAGINA[clase]} · {e}", s)
                    for clase, subs in SUBPAGINAS.items()
                    for a, e, s in subs if a != clase]


def _etq(v):
    clase = {"EXCLUIDA": "no", "NO EXCLUIDA": "si", "NO EVALUABLE": "nn"}[v]
    texto = {"EXCLUIDA": "NO", "NO EXCLUIDA": "SÍ", "NO EVALUABLE": "??"}[v]
    return f'<span class="etq {clase}">{texto}</span>'


def _figuras_html(lista, base=""):
    trozos = []
    for archivo, titulo, lectura in lista:
        if not (AQUI / archivo).exists():
            continue
        trozos.append(
            f'<figure class="figura"><img src="{base}img/{archivo}?v={_BUILD_TS}" '
            f'alt="{html.escape(titulo)}" loading="lazy">'
            f'<figcaption><b>{html.escape(titulo)}</b>'
            f'<span>{html.escape(lectura)}</span></figcaption></figure>')
    return "\n".join(trozos)


def _bloques_figuras(comunes, propias=(), sin_comunes=""):
    """Los dos bloques de figuras de una página de clase, siempre en este orden.

    `sin_comunes` es el motivo que se imprime cuando el set común todavía no se
    puede generar para esa página (falta el catálogo de diseños).
    """
    trozos = ["<h2>Figuras</h2>"]
    if comunes:
        trozos.append('<p class="rev">Set común: las mismas tres figuras en '
                      'todas las clases.</p>')
        trozos.append(_figuras_html(comunes))
    elif sin_comunes:
        trozos.append(f'<div class="aviso"><strong>Set común no disponible.'
                      f'</strong> {sin_comunes}</div>')
    if propias:
        trozos.append("<h2>Específicas de esta clase</h2>")
        trozos.append('<p class="rev">Solo existen para esta clase: no hay dato '
                      'para replicarlas en las demás.</p>')
        trozos.append(_figuras_html(propias))
    return "\n".join(trozos)


def envoltura(archivo, titulo, subtitulo, cuerpo, fecha, base="", portada=False):
    """Cabecera, navegación y pie comunes a todas las páginas."""
    clase_activa = _CLASE_DE.get(archivo)

    def _enlace(a, n):
        activo = (a == archivo) or (a == clase_activa)
        return (f'<li><a href="{base}{a}"'
                + (' class="activo"' if activo else "")
                + f'>{html.escape(n)}</a></li>')

    nav = "".join(_enlace(a, n) for a, n, _ in PAGINAS)

    subnav = ""
    if clase_activa:
        enlaces = "".join(
            f'<a href="{base}{a}"'
            + (' class="activo"' if a == archivo else "")
            + f'>{html.escape(e)}</a>'
            for a, e, _ in SUBPAGINAS[clase_activa])
        subnav = (f'<nav class="subbarra"><div class="envoltorio">{enlaces}'
                  f'</div></nav>')

    if portada:
        cabecera = f"""<header class="principal"><div class="envoltorio">
  <p class="expediente">{html.escape(PROYECTO['expediente'])}</p>
  <h1>{html.escape(PROYECTO['titulo'])}</h1>
  <p class="sub">{html.escape(PROYECTO['subtitulo'])}</p>
  <div class="meta">
    <span>{html.escape(PROYECTO['investigadores'])}</span>
    <span>{html.escape(PROYECTO['farmaco'])}</span>
    <span>{html.escape(PROYECTO['enfermedad'])}</span>
    <span>Generado el {fecha}</span>
  </div>
</div></header>"""
    else:
        cabecera = f"""<header class="pagina"><div class="envoltorio">
  <p class="expediente">{html.escape(PROYECTO['expediente'])}</p>
  <h1>{html.escape(titulo)}</h1>
  <p class="sub">{html.escape(subtitulo)}</p>
</div></header>"""

    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="robots" content="noindex, nofollow">
<title>{html.escape(titulo)} · {html.escape(PROYECTO['expediente'])}</title>
<link rel="stylesheet" href="{base}estilo.css">
</head>
<body>
{cabecera}
<nav class="barra"><ul>{nav}</ul></nav>
{subnav}
<main><div class="envoltorio">
{cuerpo}
</div></main>
<footer><div class="envoltorio">
  <p>{html.escape(PROYECTO['expediente'])} · documento interno · resultados
  preliminares, no publicados.</p>
  <p>Generado el {fecha} · <code>sh correr.sh web</code></p>
</div></footer>
</body></html>
"""


# =============================================================================
#  CUERPO DE CADA PÁGINA
# =============================================================================

def cuerpo_index(d):
    u, den, lip = d["umbrales"], d["dendrimero"], d["liposoma"]
    bloques = "".join(
        f'<a class="bloque" href="{a}"><span class="n">{k:02d}</span>'
        f'<b>{html.escape(n)}</b><span>{html.escape(desc)}</span></a>'
        for k, (a, n, desc) in enumerate(PAGINAS[1:], start=1))
    return f"""
<div class="tarjetas">
  <div class="tarjeta azul"><div class="rotulo">Tamiz del glicocálix</div>
    <div class="cifra">≤ {u['glicocalix']:.2f} nm</div></div>
  <div class="tarjeta rojo"><div class="rotulo">Envolvimiento de membrana</div>
    <div class="cifra">≥ {u['envolvimiento']:.2f} nm</div></div>
  <div class="tarjeta verde"><div class="rotulo">Separación entre las dos</div>
    <div class="cifra">{d['margen_ventanas']:.2f} nm</div></div>
  <div class="tarjeta"><div class="rotulo">Captación fagocítica</div>
    <div class="cifra">≥ {u['fagocitosis']:.0f} nm</div></div>
  <div class="tarjeta verde"><div class="rotulo">Suelo del liposoma</div>
    <div class="cifra">≥ {lip['min_con_nucleo'][1]:.1f} nm</div></div>
  <div class="tarjeta morado"><div class="rotulo">Ventana del dendrímero</div>
    <div class="cifra">{den['suelo']:.2f} – {den['techo']:.2f} nm</div></div>
  <div class="tarjeta"><div class="rotulo">Suelo del polímero macizo</div>
    <div class="cifra">{d['polimero']['suelos'][2][3]:.2f} nm</div></div>
</div>

<h2>Resultados</h2>
<p>Los dos rangos de tamaño teóricos para el paso pasivo hacia el cerebro no se solapan: el tamiz del
glicocálix exige un diámetro de hasta {u['glicocalix']:.2f} nm y el envolvimiento de membrana exige al
menos {u['envolvimiento']:.2f} nm. La separación entre ambos rangos es de {d['margen_ventanas']:.2f} nm
— un resultado geométrico firme dentro del barrido de sensibilidad aplicado. Sin embargo, superar el
tamaño del tamiz ya no se cuenta como exclusión automática: dos medidas reales en barrera
hematoencefálica (una en tejido nativo, otra en cultivo) muestran nanopartículas más grandes que ese
poro cruzando de forma medible, así que el simulador reporta ese caso como dato insuficiente para
decidir, no como bloqueo, salvo que otro criterio (carga superficial, tiempos de tránsito) sí decida
por su cuenta.</p>
<p>De las formulaciones reales evaluadas:</p>
{_resumen_catalogo_prosa(d["catalogo_real"])}
<p>El dendrímero no puede fabricarse por encima de {den['techo']:.2f} nm de diámetro (límite medido de
la química PAMAM): por encima de ese máximo, la membrana ya no puede envolverlo (necesita al menos
{u['envolvimiento']:.2f} nm). El margen es de apenas {den['margen']:.2f} nm.</p>
<p><b>Empates técnicos.</b> Dos de estos resultados quedan dentro del margen de duda del propio
barrido y no deben citarse como cerrados: la separación entre los dos rangos de tamaño cae a 0.26 nm
en el escenario más permisivo del barrido (κ = 15 kT, Hamaker = 6.5·10⁻²¹ J), y el margen del máximo
del dendrímero ({den['margen']:.2f} nm) es menor que la incertidumbre propia de esa medida (±5 %).</p>

<h2>Alcance y limitaciones del simulador</h2>
<p>Este simulador es una herramienta de cribado con trazabilidad completa a fuente científica
publicada — no un modelo predictivo validado. La diferencia importa: cribado significa que cada
criterio aplica un umbral físico medido y publicado a un diseño concreto; predicción validada
significaría que el resultado del modelo se contrastó contra una medida real de adhesión o
permeación de ese mismo tipo de diseño, y ese contraste todavía no existe para ninguna de las
formulaciones reales evaluadas.</p>
<p>Sí existen ya dos medidas reales en barrera hematoencefálica que sirvieron para corregir el criterio
simple de tamaño del glicocálix y de entrada por caveola (ver Resultados): mostraron que el tamaño
solo, sin más información, no basta para excluir un diseño. Lo que todavía no existe es una medida
real de adhesión o permeación específicamente a la escala de un liposoma que confirme o descarte el
modelo de energía más fino que el simulador usa para el glicocálix. Sin ese dato, ese criterio más
fino sigue sin poder conectarse a los veredictos finales.</p>
<p>No fue posible determinar con confianza cuánto fármaco entrega cada diseño una vez dentro del
cerebro. Ese cálculo se investigó a fondo, presentó inconsistencias al verificarlo, y se decidió
reportarlo como dato desconocido en vez de publicar una cifra sin verificar.</p>
<p>La hipótesis original del proyecto incluía el coeficiente de partición octanol-agua (LogP) del
nanotransportador como factor de diseño. Ese factor no se implementó: el LogP es una propiedad bien
definida para una molécula (la del fingolimod, el fármaco, ya usada en otra parte del análisis), pero
no existe una física publicada que conecte el LogP de un nanotransportador completo —como un
liposoma— con su capacidad de adherirse o atravesar la barrera. Añadirlo sin una fuente científica
que lo respalde habría significado inventar un mecanismo, algo que este proyecto decidió no hacer.</p>
<p>El modelo cubre únicamente las etapas de adhesión y envolvimiento del nanotransportador en la
superficie del vaso sanguíneo cerebral. Las etapas posteriores —el paso activo hacia el otro lado de
la barrera, el tránsito celular y la liberación final del fármaco— dependen de energía celular y
quedaron fuera del alcance de este modelo desde su planteamiento inicial.</p>

<h2>Secciones</h2>
<div class="indice">{bloques}</div>
"""


def _resumen_catalogo_prosa(filas_catalogo):
    """Un párrafo por formulación real, en prosa, con su veredicto dominante.

    Reemplaza la vieja tabla de badges. Para cada formulación: si TODAS las
    rutas están EXCLUIDAs, dice por qué (la razón de la primera). Si hay
    mezcla, dice cuántas rutas excluye y por qué, y qué falta en las demás.
    """
    trozos = []
    for f in filas_catalogo:
        veredictos = [f["rutas"][n]["veredicto"] for n in f["rutas"]]
        n_excluidas = veredictos.count("EXCLUIDA")
        n_total = len(veredictos)
        nombre = html.escape(f["nombre"])
        ficha = (f"{f['diametro']:.1f} nm, ζ {f['zeta']:+.2f} mV")
        if n_excluidas == n_total:
            motivo = next((r["muere"] for r in f["rutas"].values() if r["muere"]), "—")
            trozos.append(
                f'<p><b>{nombre}</b> ({ficha}) queda excluida en las {n_total} vías '
                f'evaluadas: {html.escape(motivo)}.</p>')
        elif n_excluidas == 0:
            trozos.append(
                f'<p><b>{nombre}</b> ({ficha}) no tiene resultado en ninguna vía: '
                f'faltan datos en varias compuertas.</p>')
        else:
            razones = sorted({html.escape(r["muere"]) for r in f["rutas"].values()
                               if r["veredicto"] == "EXCLUIDA" and r["muere"]})
            trozos.append(
                f'<p><b>{nombre}</b> ({ficha}) queda excluida en {n_excluidas} de '
                f'{n_total} vías ({"; ".join(razones)}); en las demás no hay resultado '
                f'por falta de dato.</p>')
    return "\n".join(trozos)


def cuerpo_liposoma(d):
    """SOLO lo medido: el suelo de Pan 2008 y los liposomas publicados."""
    u, lip = d["umbrales"], d["liposoma"]
    tabla_veredictos, tabla_muere = _tablas_catalogo(d["catalogo_real"],
                                                     d["rutas"])

    tabla_bicapa = (
        '<div class="tabla-scroll"><table><thead><tr>'
        '<th>t bicapa medido (nm)</th>'
        '<th class="num">Ø mín. con núcleo de 4 nm</th>'
        '<th class="num">Ø mín. con núcleo nulo</th></tr></thead><tbody>'
        + "".join(f'<tr><td>{t:.1f}</td><td class="num">{a:.1f}</td>'
                  f'<td class="num">{b:.1f}</td></tr>'
                  for t, a, b in zip(lip["t_bicapa"], lip["min_con_nucleo"],
                                     lip["min_nucleo_nulo"]))
        + '</tbody></table></div>')

    return f"""
<h2>Suelo geométrico</h2>
<p class="formula"><code>d_externo = d_núcleo + 2 · t_bicapa</code></p>
{tabla_bicapa}
<div class="tarjetas">
  <div class="tarjeta verde"><div class="rotulo">Ø mínimo real</div>
    <div class="cifra">{lip['min_con_nucleo'][1]:.1f} nm</div>
    <div class="pie">Núcleo utilizable de 4 nm.</div></div>
  <div class="tarjeta rojo"><div class="rotulo">Tamiz del glicocálix</div>
    <div class="cifra">{u['glicocalix']:.2f} nm</div>
    <div class="pie">El suelo ya está por encima.</div></div>
  <div class="tarjeta"><div class="rotulo">Caso imposible: núcleo nulo</div>
    <div class="cifra">{lip['min_nucleo_nulo'][0]:.1f} nm</div>
    <div class="pie">No encapsularía nada. A {lip['limite']['margen_nm']:.2f} nm del umbral.</div></div>
</div>
<p class="rev">Fuente: Pan et al. 2008, PRL 100:198103, Fig. 3c
(t bicapa medido, {lip['t_bicapa'][0]:.1f}–{lip['t_bicapa'][-1]:.1f} nm).</p>

<h2>Liposomas publicados</h2>
<p class="rev">Ø y ζ medidos, con la referencia bajo cada nombre.</p>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Primera compuerta que falla</h2>
{tabla_muere}

<div class="aviso"><strong>Ruta B.</strong> Muselman 2026 (700 nm) era el único diseño
no excluido del modelo. Desde el 2026-08-12 sale <b>no evaluable</b>: la compuerta B.3
volvió a abrirse porque el tránsito del monocito al cerebro inflamado
({R.T_TRANSITO_PRIMERA_DETECCION_h:.0f}–{R.T_TRANSITO_PICO_h:.0f} h, Tong 2016) se solapa
con la descarga del liposoma (&gt; {R.T_LIBERACION_COTA_INFERIOR_h:.0f} h sin techo medido,
Mao 2014).</div>

<h2>Detalle por liposoma</h2>
<p class="rev">Clic sobre cada fila para desplegar sus 3 figuras propias
(ventanas, matriz, recorrido) y el desglose de compuertas por ruta.</p>
{_bloques_detalle_liposoma(d["catalogo_real_detalle"])}
"""


def _tablas_catalogo(filas_catalogo, nombres_rutas):
    """Las dos tablas de un catálogo: veredictos y primera compuerta que falla."""
    filas_v, filas_m = [], []
    for f in filas_catalogo:
        cv, cm = [], []
        for n in nombres_rutas:
            r = f["rutas"][n]
            marca = "!" * min(r["salvedades"], 3)
            cv.append(f'<td>{_etq(r["veredicto"])}'
                      + (f' <b>{marca}</b>' if marca else '') + '</td>')
            if r["veredicto"] == "EXCLUIDA":
                cm.append(f'<td><span class="rev">{html.escape(r["muere"] or "—")}</span></td>')
            elif r["veredicto"] == "NO EXCLUIDA":
                cm.append('<td><span class="etq si">no excluida</span></td>')
            else:
                falta = ", ".join(r["faltan"][:2]) or "—"
                cm.append(f'<td><span class="rev">falta: {html.escape(falta)}</span></td>')
        nota = f'<div class="rev">{html.escape(f["nota"])}</div>' if f["nota"] else ""
        filas_v.append(f'<tr><td><b>{html.escape(f["nombre"])}</b>{nota}</td>'
                       f'<td class="num">{f["diametro"]:.1f}</td>'
                       f'<td class="num">{f["zeta"]:+.2f}</td>'
                       f'<td class="num">{f["peg"]:.0f}</td>' + "".join(cv) + '</tr>')
        filas_m.append(f'<tr><td><b>{html.escape(f["nombre"])}</b>{nota}</td>'
                       f'<td class="num">{f["diametro"]:.1f}</td>' + "".join(cm) + '</tr>')

    cab = "".join(f'<th>{html.escape(n.split(" (")[0])}</th>' for n in nombres_rutas)
    tabla_veredictos = (
        '<div class="tabla-scroll"><table><thead><tr><th>Formulación</th>'
        '<th class="num">Ø (nm)</th><th class="num">ζ (mV)</th>'
        f'<th class="num">PEG (nm)</th>{cab}</tr></thead><tbody>'
        + "".join(filas_v) + '</tbody></table></div>')
    tabla_muere = (
        '<div class="tabla-scroll"><table><thead><tr><th>Formulación</th>'
        f'<th class="num">Ø (nm)</th>{cab}</tr></thead><tbody>'
        + "".join(filas_m) + '</tbody></table></div>')
    return tabla_veredictos, tabla_muere


def cuerpo_liposoma_teoricos(d):
    """Los tres diseños teóricos de liposoma, y solo ellos."""
    tabla_veredictos, tabla_muere = _tablas_catalogo(d["catalogo_teorico"],
                                                     d["rutas"])
    return f"""
<div class="aviso"><strong>Diseños teóricos.</strong> Formulaciones propuestas:
ni el Ø ni el ζ están medidos. Los liposomas publicados están en la pestaña
<b>Reales</b>.</div>

<h2>Veredictos</h2>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Primera compuerta que falla</h2>
{tabla_muere}

<div class="aviso"><strong>Zona muerta 31–50 nm.</strong> Los tres: demasiado
grandes para el glicocálix, demasiado pequeños para el macrófago.</div>

<div class="aviso"><strong>Tensión de carga.</strong> ζ positivo hace falta para
entrar (ruta C) e impide difundir (Nance 2012: nada con ζ &lt; −6 mV difunde;
100 % de las carboxiladas inmovilizadas, incluidas las de 40 nm). Por el lado
positivo el estorbo también está medido, pero más arriba: con ζ
+{R.ZETA_ADHESIVO_POSITIVO_mV:.1f} mV difunde menos del 10 % de la población
(Berry 2016) y con +35.3 mV la partícula queda inmovilizada (Mastorakos 2016),
las dos por seguimiento de partículas en cerebro de rata <em>ex vivo</em>. Los
tres diseños (+2.0, +5.0, +6.7 mV) caen en el hueco sin dato que va de 0 a
+{R.ZETA_ADHESIVO_POSITIVO_mV:.0f} mV, así que su compuerta sigue
DESCONOCIDA: no se extrapola.</div>

<div class="aviso"><strong>Salvedad de esos dos datos.</strong> Berry y
Mastorakos son polímero/ADN, no liposomas, y su ζ está medido en NaCl 10 mM a
pH 7.0, no en aCSF. En aCSF las dos formulaciones catiónicas pierden la
estabilidad coloidal, así que la inmovilización puede deberse a adhesión
electrostática o a obstrucción estérica por agregación: los experimentos no
separan las dos causas.</div>

<h2>Detalle por liposoma</h2>
<p class="rev">Clic sobre cada fila para desplegar sus 3 figuras propias
(ventanas, matriz, recorrido) y el desglose de compuertas por ruta.</p>
{_bloques_detalle_liposoma(d["catalogo_teorico_detalle"])}
"""


def _etq_gate(estado):
    clase = {"PASA": "si", "FALLA": "no", "DESCONOCIDA": "nn"}.get(estado, "nn")
    return f'<span class="etq {clase}">{html.escape(estado)}</span>'


def _fmt_num(v):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.3g}"
    return str(v)


def _bloques_detalle_liposoma(detalle):
    """Una tarjeta <details> por diseño: sus 3 figuras propias (ventanas,
    matriz, recorrido) + desglose de compuertas por ruta + insignia
    Fabricable. Patrón único, compartido por Reales/Teóricos/Dataset desde
    el 2026-08-18k (antes solo lo tenía Dataset)."""
    bloques = []
    for it in detalle:
        figs = "".join(
            f'<div class="figura"><img src="{it["img"]}_{suf}.png?v={_BUILD_TS}" '
            f'alt="{suf}" loading="lazy"></div>'
            for suf in ("ventanas", "matriz", "recorrido"))
        secciones_rutas = []
        for nombre_ruta, info in it["rutas"].items():
            filas_c = "".join(
                f'<tr><td>{html.escape(c["nombre"])}</td>'
                f'<td>{_etq_gate(c["estado"])}</td>'
                f'<td class="num">{_fmt_num(c["valor"])}</td>'
                f'<td class="num">{_fmt_num(c["umbral"])}</td>'
                f'<td>{html.escape(c["unidad"] or "")}</td>'
                f'<td class="rev">{html.escape(c["fuente"] or "—")}</td></tr>'
                for c in info["compuertas"])
            secciones_rutas.append(
                f'<h4>{html.escape(nombre_ruta)} — {_etq(info["veredicto"])}</h4>'
                '<div class="tabla-scroll"><table><thead><tr><th>Compuerta</th>'
                '<th>Estado</th><th class="num">Valor</th><th class="num">Umbral</th>'
                f'<th>Unidad</th><th>Fuente</th></tr></thead><tbody>{filas_c}</tbody></table></div>')
        fab = it.get("fabricable")
        insignia_fab = (f' · <span class="etq {"si" if fab else "no"}">'
                        f'Fabricable: {"SÍ" if fab else "NO"}</span>'
                        if fab is not None else "")
        nota = f' <span class="rev">({html.escape(it["nota"])})</span>' if it.get("nota") else ""
        bloques.append(f'''
<details class="liposoma">
<summary>{html.escape(it["nombre"])}{nota} — Ø {it["diametro"]:.1f} nm ·
ζ {it["zeta"]:+.2f} mV · PEG {it["peg"]:.2f} nm{insignia_fab}</summary>
<div class="dataset50-figs">{figs}</div>
{"".join(secciones_rutas)}
</details>''')
    return "".join(bloques)


def cuerpo_dataset_50(d):
    """50 liposomas sintéticos (semilla 42), último día 2026-08-18. Rango de
    diámetro/ζ/PEG tomado de CATALOGO_REAL + teóricos. NO citables como
    predicción — mismo criterio que los diseños teóricos. La carga útil (G.2)
    queda DESCONOCIDA en los 50, por decisión vigente del 2026-08-13."""
    tabla_veredictos, tabla_muere = _tablas_catalogo(d["dataset_50"], d["rutas"])

    return f"""
<div class="aviso"><strong>50 diseños SINTÉTICOS (semilla 42).</strong>
Combinaciones aleatorias de diámetro, ζ y PEG dentro del rango que cubren los
diseños reales ya validados del proyecto (Mao 2014, Gong 2022, Chow 2025,
Muselman 2026) más los teóricos. No son mediciones ni predicciones citables.
La carga útil (G.2) queda DESCONOCIDA en los 50, por decisión vigente del
2026-08-13 — no se fijó ningún valor nuevo bajo presión de tiempo.</div>

<h2>Veredictos</h2>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Detalle por liposoma</h2>
<p class="rev">Clic sobre cada fila para desplegar sus 3 figuras propias
(ventanas, matriz, recorrido) y el desglose de compuertas por ruta.</p>
{_bloques_detalle_liposoma(d["dataset_50_detalle"])}
"""


def cuerpo_dendrimero(d):
    u, den = d["umbrales"], d["dendrimero"]
    fd = "".join(
        f'<tr><td><b>G{g["gen"]}</b></td><td class="num">{g["diametro"]:.2f}</td>'
        f'<td>{_etq("NO EXCLUIDA") if g["estado"] == "PASA" else _etq("NO EVALUABLE")}</td>'
        f'<td>{"sí" if g["glicocalix"] else "no"}</td>'
        f'<td>{"sí" if g["envuelve"] else "no"}</td></tr>'
        for g in den["generaciones"])
    tabla = ('<div class="tabla-scroll"><table><thead><tr><th>Generación</th>'
             '<th class="num">Ø medido (nm)</th><th>Fabricable</th>'
             '<th>¿Pasa el glicocálix?</th><th>¿Se envuelve?</th></tr></thead>'
             '<tbody>' + fd + '</tbody></table></div>')

    return f"""
<h2>Ventana geométrica</h2>
<div class="tarjetas">
  <div class="tarjeta morado"><div class="rotulo">Suelo · G3</div>
    <div class="cifra">{den['suelo']:.2f} nm</div>
    <div class="pie">Ø medido y alojamiento del fármaco demostrado.</div></div>
  <div class="tarjeta morado"><div class="rotulo">Techo · G10</div>
    <div class="cifra">{den['techo']:.2f} nm</div>
    <div class="pie">Última generación completable.</div></div>
  <div class="tarjeta rojo"><div class="rotulo">Margen al envolvimiento</div>
    <div class="cifra">{den['margen']:.2f} nm</div>
    <div class="pie">Exige {u['envolvimiento']:.2f} nm. No lo alcanza.</div></div>
</div>
{tabla}
<p class="rev">Ø medidos: Prosa 2001, Tabla 1 (SAXS en metanol, ±5 %).
Techo: Maiti 2004, Tabla 4 y Figs. 21-22 · de Gennes &amp; Hervet 1983.
Suelo: Devarakonda 2004, Tabla 2.</p>

<div class="aviso"><strong>Empate técnico.</strong> ±5 % sobre los
{den['techo_medido']:.2f} nm de G10 da {den['techo_medido']*(1-den['precision']):.2f}–{den['techo_medido']*(1+den['precision']):.2f} nm,
y {den['techo_medido']*(1+den['precision']):.2f} &gt; {u['envolvimiento']:.2f}.
Medido en metanol, no en agua.</div>

<div class="aviso"><strong>Carga 1:1.</strong> Devarakonda 2004 (G0–G3): una
molécula de fármaco por dendrímero. Sin dato para G4–G10 ni para fingolimod.</div>

<div class="aviso"><strong>ζ.</strong> En el catálogo va como positivo sin valor:
no hay ζ medido en fuente primaria para esta clase.</div>

{_bloques_figuras(FIGURAS_DEND_REALES, FIGURAS_DENDRIMERO_PROPIAS)}
"""


def _tablas_teoricos(filas_teoricos, con_clase=False):
    """Tabla de veredictos + desglose «dónde choca cada uno», para cualquier
    catálogo de transportadores teóricos."""
    def _e(estado):
        v = {"PASA": "NO EXCLUIDA", "FALLA": "EXCLUIDA",
             "DESCONOCIDA": "NO EVALUABLE"}[estado]
        return _etq(v)

    filas = []
    for t in filas_teoricos:
        celdas = "".join(f"<td>{_e(p['estado'])}</td>" for p in t["puertas"])
        col_clase = (f'<td>{html.escape(t["clase"])}</td>' if con_clase else "")
        filas.append(
            f'<tr><td><b>{html.escape(t["nombre"])}</b></td>'
            f'<td class="num">{t["diametro"]:.2f}</td>'
            f'<td class="num">{t["zeta"]:+.1f}</td>{col_clase}{celdas}</tr>')
    cab_clase = "<th>Clase</th>" if con_clase else ""
    tabla = ('<div class="tabla-scroll"><table><thead><tr><th>Transportador</th>'
             '<th class="num">Ø (nm)</th><th class="num">ζ (mV)</th>'
             f'{cab_clase}'
             '<th>Fabricable</th><th>Glicocálix</th><th>Envolvimiento</th>'
             '</tr></thead><tbody>' + "".join(filas) + '</tbody></table></div>')

    porque = []
    for t in filas_teoricos:
        lineas = "".join(
            f'<li><b>{html.escape(p["nombre"])}</b> — {p["estado"]}'
            + (f' ({p["valor"]:g} vs {p["umbral"]:.2f} nm)'
               if p["valor"] is not None and p["umbral"] is not None else "")
            + (f'<br><span class="rev">{html.escape(p["motivo"])}</span>'
               if p["motivo"] else "")
            + (f'<br><span class="rev">Salvedad: '
               f'{html.escape(p.get("advertencia") or "")}</span>'
               if p.get("advertencia") else "")
            + '</li>'
            for p in t["puertas"])
        porque.append(f'<h3>{html.escape(t["nombre"])}</h3><ul>{lineas}</ul>')
    return tabla, "".join(porque)


def cuerpo_dendrimero_teoricos(d):
    tabla, porque = _tablas_teoricos(d["teoricos"])
    return f"""
<div class="aviso"><strong>Datos sintéticos.</strong> Los valores de estas tres
fichas son inventados dentro de un rango teórico plausible. No son medidas, no
cierran ninguna tarea de verificación y no pueden citarse como resultado.</div>

{tabla}

<h2>Dónde choca cada uno</h2>
{porque}

<div class="aviso"><strong>ζ.</strong> Ninguno de los tres entra en el modelo:
sigue sin haber ζ medido en fuente primaria para esta clase.</div>

<div class="aviso"><strong>PPI y carbosilano.</strong> La ventana suelo/techo del
código está derivada solo de PAMAM, así que otra química devuelve DESCONOCIDA en
lugar de heredarla. Tarea G.1a-bis.</div>

<div class="aviso"><strong>El acople no ordena los tres.</strong> ΔG =
−2.303·k<sub>B</sub>T·logP depende solo del logP del fármaco (4.16), así que vale
−9.58 k<sub>B</sub>T en los tres. Las fracciones acopladas que da el simulador de
acople son fracción de volumen de la caja, no una propiedad del transportador.</div>

{_bloques_figuras(FIGURAS_DEND_TEORICOS, FIGURAS_DENDRIMERO_TEORICOS_PROPIAS)}
"""


# Las tres figuras del set común se dibujan sobre un CATÁLOGO de diseños, y para
# el polímero todavía no hay ninguno: haría falta ζ, y no hay ζ de polímero
# medido en fuente primaria. Inventarlo rompería la regla del proyecto.
_SIN_CATALOGO_POLIMERO = (
    "las tres figuras del set común se dibujan sobre un catálogo de diseños y "
    "todavía no hay ninguno de esta clase. Cada diseño necesita Ø, ζ y masa "
    "molar; el ζ del polímero no está medido en fuente primaria y no se "
    "inventa.")


def cuerpo_polimero_teoricos(d):
    tabla, porque = _tablas_teoricos(d["teoricos_polimero"], con_clase=True)
    notas = "".join(
        f'<li><b>{html.escape(t["nombre"])}</b> — '
        f'{html.escape(t["nota"].replace("DATO SINTÉTICO · ", ""))}</li>'
        for t in d["teoricos_polimero"])
    return f"""
<div class="aviso"><strong>Datos sintéticos.</strong> Los valores de estas tres
fichas son inventados dentro de un rango teórico plausible. No son medidas, no
cierran ninguna tarea de verificación y no pueden citarse como resultado.</div>

<div class="aviso"><strong>Dos de las tres son micelas, no polímero macizo.</strong>
Una micela está hecha de polímero pero es un agregado autoensamblado de núcleo y
corona, con concentración micelar crítica por debajo de la cual se deshace. Su
suelo no es el glóbulo de una cadena colapsada sino el que fijan el número de
agregación y la carga de fármaco, así que usa su propia compuerta.</div>

{tabla}
<ul class="rev">{notas}</ul>

<h2>Dónde choca cada uno</h2>
{porque}

<div class="aviso"><strong>Suelo de la micela: {d['micela_suelo_nm']:.1f} nm</strong>
(Sochor 2020, SANS, micela cargada 10/1). Vacía mide {d['micela_vacia_nm']:.1f} nm:
el fármaco no es un pasajero, multiplica el diámetro por 3.6. Por debajo del suelo
la compuerta devuelve DESCONOCIDA, no FALLA, porque el dato viene de otro polímero
y otro fármaco.</div>

<div class="aviso"><strong>Mw del PLGA: {d['plga_mw_kDa']:.1f} kDa, DERIVADA.</strong>
La ficha no la declara. Sale de sus propios datos (462 monómeros por cadena)
suponiendo razón 85:15, la única con densidad medida. La razón casi no importa:
con 50:50 el suelo pasa de 4.42 a 4.31 nm y ningún veredicto cambia.</div>

<div class="aviso"><strong>ζ.</strong> Ninguno de los tres entra en el modelo:
sigue sin haber ζ medido en fuente primaria para estas clases.</div>

{_bloques_figuras(FIGURAS_POL_TEORICOS, FIGURAS_POLIMERO_TEORICOS_PROPIAS)}
"""


def cuerpo_polimero(d):
    u, pol = d["umbrales"], d["polimero"]
    fs = "".join(
        f'<tr><td><b>{html.escape(e)}</b></td><td class="num">{mw}</td>'
        f'<td class="num">{rho:.2f}</td><td class="num"><b>{dd:.2f}</b></td>'
        f'<td>{"sí" if dd <= u["glicocalix"] else "no"}</td></tr>'
        for e, mw, rho, dd in pol["suelos"])
    tabla_suelos = (
        '<div class="tabla-scroll"><table><thead><tr><th>Polímero</th>'
        '<th class="num">Mw (kDa)</th><th class="num">ρ (g/cm³)</th>'
        '<th class="num">Suelo (nm)</th><th>¿Pasa el glicocálix?</th>'
        '</tr></thead><tbody>' + fs + '</tbody></table></div>')

    base = pol["sensibilidad"][0][1]
    fsen = "".join(
        f'<tr><td>{"solo la cadena" if n == 0 else f"+ {n} molécula(s) de fingolimod"}</td>'
        f'<td class="num">{dd:.4f}</td>'
        f'<td class="num">{"—" if n == 0 else f"{100*(dd/base-1):+.2f} %"}</td></tr>'
        for n, dd in pol["sensibilidad"])
    tabla_sen = ('<div class="tabla-scroll"><table><thead><tr><th>Contenido</th>'
                 '<th class="num">d (nm)</th><th class="num">cambio</th>'
                 '</tr></thead><tbody>' + fsen + '</tbody></table></div>')

    return f"""
<h2>Polímero macizo — suelo del glóbulo colapsado</h2>
<p class="formula"><code>V = M / (ρ · N_A)</code> &nbsp;·&nbsp;
<code>d = (6V / π)^(1/3)</code></p>
<p class="rev">Una partícula maciza no puede ser más pequeña que una sola cadena
del polímero colapsada sobre sí misma. Sin techo arquitectónico.</p>
{tabla_suelos}
<p class="rev">Densidades: Parker et al. 2010, Biomed Mater 5:055004, Tabla 2
(derivadas por los autores de su propia velocidad del sonido e impedancia).
La Tabla 1 del mismo artículo da densidades del fabricante y no se usa.</p>

<h3>El fármaco casi no mueve el suelo</h3>
{tabla_sen}
<p class="rev">Sobre PLGA de 53 kDa. Comprobado, no supuesto: por eso el suelo se
calcula solo con la cadena y no arrastra el tamaño del fingolimod, que hasta el
2026-08-13 no tenía
fuente.</p>

{_bloques_figuras([], FIGURAS_POLIMERO, sin_comunes=_SIN_CATALOGO_POLIMERO)}
"""


def cuerpo_metodo(d):
    fc = "".join(
        f'<tr><td><b>{html.escape(c["nombre"])}</b></td>'
        + (f'<td>{html.escape(c["fuente"][:150])}{"…" if len(c["fuente"]) > 150 else ""}</td>'
           if c["fuente"] != "—" else '<td><i>sin fuente primaria verificada</i></td>')
        + f'<td>{html.escape(c["motivo"][:120])}{"…" if len(c["motivo"]) > 120 else ""}</td></tr>'
        for c in d["compuertas"])
    tabla = ('<div class="tabla-scroll"><table><thead><tr><th>Compuerta</th>'
             '<th>Anclaje</th><th>Observación</th></tr></thead><tbody>'
             + fc + '</tbody></table></div>')
    rutas = "".join(f'<li>{html.escape(n)}</li>' for n in d["rutas"])

    return f"""
<h2>Qué significa cada veredicto</h2>
<div class="tarjetas">
  <div class="tarjeta rojo"><div class="rotulo">EXCLUIDA</div>
    <div class="cifra" style="font-size:1.05rem">Afirmación fuerte</div>
    <div class="pie">No puede usar esa ruta.</div></div>
  <div class="tarjeta verde"><div class="rotulo">NO EXCLUIDA</div>
    <div class="cifra" style="font-size:1.05rem">Afirmación débil</div>
    <div class="pie">Candidato. No es predicción de éxito.</div></div>
  <div class="tarjeta"><div class="rotulo">NO EVALUABLE</div>
    <div class="cifra" style="font-size:1.05rem">Falta información</div>
    <div class="pie">Sin dato no se da por superada.</div></div>
</div>

<div class="aviso"><strong>Techo estructural.</strong> La transcitosis es
transporte activo dependiente de ATP: fuera del alcance de un modelo de
equilibrio, devuelve «sin dato» siempre. Está en las rutas A, C y D, así que esas
tres nunca pueden salir «no excluida».</div>

<h2>Rutas</h2>
<ul>{rutas}</ul>

<h2>Compuertas</h2>
{tabla}

<h2>Clases evaluables</h2>
<div class="tabla-scroll"><table><tbody>
<tr><td>Liposoma</td><td><span class="etq si">sí</span></td></tr>
<tr><td>Dendrímero</td><td><span class="etq si">sí</span></td></tr>
<tr><td>Polímero macizo</td><td><span class="etq si">sí</span></td></tr>
<tr><td>Micela</td><td><span class="etq nn">fuera de alcance</span></td></tr>
</tbody></table></div>
"""


def cuerpo_validacion(d):
    fv = "".join(f'<tr><td>{html.escape(k)}</td><td class="num"><b>{v}</b></td></tr>'
                 for k, v in VALIDACION["conteo"])
    return f"""
<h2>Recuento a {VALIDACION['fecha']}</h2>
<div class="tabla-scroll"><table><thead><tr><th>Categoría</th>
<th class="num">Casos</th></tr></thead><tbody>{fv}</tbody></table></div>
<p>{html.escape(VALIDACION['veredicto'])}</p>
<p class="rev">Fuente: <code>{VALIDACION['fuente']}</code> ·
<a href="fichas/resultados_validacion.html">detalle caso por caso</a> ·
<a href="fichas/protocolo_validacion.html">protocolo</a></p>
"""


def cuerpo_fichas(d):
    lfi = "".join(
        f'<div class="fuente"><div class="aut">'
        f'<a href="fichas/{Path(f).stem}.html">{html.escape(t)}</a></div>'
        f'<div class="rev">{html.escape(f)}</div></div>'
        for f, t in FICHAS if (RAIZ / "verificacion" / f).exists())
    return f"""<h2>Fichas de verificación</h2>{lfi}"""


def cuerpo_pendientes(d):
    lp = "".join(
        f'<div class="fuente"><span class="prio {p}">{p}</span> '
        f'<span class="aut">{html.escape(t)}</span>'
        f'<div class="tit">{html.escape(x)}</div></div>'
        for t, p, x in PENDIENTES)
    return f"""<h2>Pendientes</h2>{lp}"""


def cuerpo_fuentes(d):
    lf = "".join(
        f'<div class="fuente"><div class="aut">{html.escape(a)}</div>'
        f'<div class="tit">{html.escape(t)}</div>'
        f'<div class="rev">{html.escape(r)}</div>'
        f'<a href="{u_}">{u_}</a><br><span class="uso">{html.escape(uso)}</span></div>'
        for a, t, r, u_, uso in FUENTES)
    return f"""<h2>Fuentes</h2>{lf}"""


CUERPOS = {
    "index.html": cuerpo_index,
    "liposoma.html": cuerpo_liposoma,
    "liposoma_teoricos.html": cuerpo_liposoma_teoricos,
    "dataset_50.html": cuerpo_dataset_50,
    "dendrimero.html": cuerpo_dendrimero,
    "dendrimero_teoricos.html": cuerpo_dendrimero_teoricos,
    "polimero.html": cuerpo_polimero,
    "polimero_teoricos.html": cuerpo_polimero_teoricos,
    "metodo.html": cuerpo_metodo,
    "validacion.html": cuerpo_validacion,
    "fichas.html": cuerpo_fichas,
    "pendientes.html": cuerpo_pendientes,
    "fuentes.html": cuerpo_fuentes,
}


def construir_ficha(md_path, titulo, fecha):
    cuerpo = ('<a class="volver" href="../fichas.html">← volver al índice de fichas</a>'
              f'<div class="ficha">{markdown(md_path.read_text(encoding="utf-8"))}</div>'
              '<a class="volver" href="../fichas.html">← volver al índice de fichas</a>')
    return envoltura("fichas.html", titulo, md_path.name, cuerpo, fecha, base="../")


# =============================================================================
#  PRINCIPAL
# =============================================================================

def main():
    fecha = datetime.date.today().isoformat()
    d = recoger()

    # GUARDA CONTRA BIBLIOGRAFÍA INCOMPLETA. Va la PRIMERA porque es barata y
    # porque no depende de las figuras. La bibliografía es uno de los tres
    # bloques que van A MANO, así que se desincroniza sola: se añade una fuente
    # al código y nadie se acuerda de la web. El 2026-08-13 se descubrió que
    # Thorne & Nicholson 2006, que gobierna el umbral de 38 nm del fármaco
    # liberado, llevaba desde el principio sin aparecer en la página de fuentes.
    # Esto lo hace imposible: cada apellido que rutas.py cite como `fuente=` o
    # en un `_F_*` tiene que estar en FUENTES, o se aborta.
    # Se miran los TRES módulos, no solo rutas.py: el 2026-08-13 se descubrió que
    # Nagle 2017, que ancla el barrido entero de kappa, vivía en
    # envolvimiento_core.py y tampoco estaba en la bibliografía.
    # Una cadena de `fuente=` puede empezar por APELLIDO ("Nance et al. 2012...")
    # o por REVISTA ("Nanomedicine 10:393..."), así que se acepta si contiene
    # cualquier apellido O cualquier revista de FUENTES.
    _citas = []
    for _f in ("rutas.py", "glicocalix.py", "envolvimiento_core.py"):
        _citas += re.findall(r'(?:fuente=|_F_[A-Z_]+ *= *)\(?"([^"]+)"',
                             (AQUI / _f).read_text(encoding="utf-8"))
    # Se compara por APELLIDO **Y AÑO**, no solo por apellido. Si solo se mira el
    # apellido, una entrada nueva del mismo autor pasa desapercibida: fue
    # exactamente lo que ocurrió con Nagle 2017, que faltaba mientras Nagle 2008
    # sí estaba, y una guarda por apellido no lo habría detectado.
    # Dos comprobaciones SEPARADAS, y el orden importa:
    #  · por REVISTA (campo 2 de FUENTES): basta el nombre. NO se exige año,
    #    porque en «Drug Deliv Transl Res 15:2022» el 2022 es el número de
    #    PÁGINA y confundirlo con el año daba un falso positivo.
    #  · por APELLIDO (campo 0): se exige además el AÑO, porque si no, una
    #    entrada nueva del mismo autor pasa desapercibida. Fue el caso de
    #    Nagle 2017, que faltaba mientras Nagle 2008 sí estaba.
    _en_biblio_txt = " ".join(b[0] + " " + b[2] for b in FUENTES)
    _revistas, _autores = set(), []
    for _b in FUENTES:
        _revistas |= set(re.findall(r"[A-ZÁÉÍÓÚ][A-Za-zÀ-ÿ]{3,}", _b[2]))
        _autores.append((set(re.findall(r"[A-ZÁÉÍÓÚ][A-Za-zÀ-ÿ]{3,}", _b[0])),
                         set(re.findall(r"(?:19|20)\d{2}", _b[0]))))

    def _citada(cita):
        palabras = set(re.findall(r"[A-ZÁÉÍÓÚ][A-Za-zÀ-ÿ]{3,}", cita))
        if palabras & _revistas:
            return True
        anios = set(re.findall(r"(?:19|20)\d{2}", cita))
        return any(palabras & _ap and (not anios or anios & _an)
                   for _ap, _an in _autores)

    _sin_citar = sorted({c[:60] for c in _citas if not _citada(c)})

    # SEGUNDA PASADA, sobre los COMENTARIOS. La comprobación de arriba solo ve
    # las fuentes declaradas con `fuente=`, y el 2026-08-13 se descubrió que NUEVE
    # fuentes vivían solo en comentarios y ninguna estaba en la bibliografía,
    # entre ellas Sochor 2020 (que gobierna el suelo de la micela), Bastiani y
    # Parton 2010 (la compuerta de caveola) y Morris y Homann 2001 (el barrido de
    # tensión). Se buscan patrones «Apellido 2020», «Apellido et al. 2020» y
    # «Apellido & Otro 2020», y basta con que el apellido esté en la bibliografía.
    _sueltas = set()
    for _f in ("rutas.py", "glicocalix.py", "envolvimiento_core.py"):
        for _ap, _an in re.findall(
                r"([A-Z][a-zÀ-ÿ]{3,})\s+(?:et al\.?,?\s*|& *[A-Z][a-zÀ-ÿ]+ *)?"
                r"((?:19|20)\d{2})",
                (AQUI / _f).read_text(encoding="utf-8")):
            if _ap not in _en_biblio_txt:
                _sueltas.add(f"{_ap} {_an}  (citado en un comentario)")
    _sin_citar += sorted(_sueltas)
    if _sin_citar:
        print()
        print("  ABORTADO · FUENTES DEL CÓDIGO QUE NO ESTÁN EN LA BIBLIOGRAFÍA")
        for a in _sin_citar:
            print(f"    · {a}")
        print("  Añádelas a FUENTES en construir_web.py")
        print()
        raise SystemExit(5)

    (SALIDA / "img").mkdir(parents=True, exist_ok=True)
    (SALIDA / "fichas").mkdir(parents=True, exist_ok=True)

    # GUARDA CONTRA FIGURAS VIEJAS. Este script calcula los NÚMEROS ejecutando
    # el simulador, pero las figuras son PNG en disco y solo las copia. Hasta el
    # 2026-08-12 cada grupo de figuras se generaba con su propio comando, así
    # que la web mezclaba PNG de momentos distintos y se contradecían entre sí
    # y con sus propias tablas. Bug detectado por Jhovan. Lo normal es entrar
    # por `correr.sh web`, que ya regenera todo antes; esto cubre el caso de
    # ejecutar construir_web.py a pelo.
    _codigo = max((AQUI / f).stat().st_mtime
                  for f in ("rutas.py", "glicocalix.py", "envolvimiento_core.py",
                            "envolvimiento_script.py", "construir_web.py")
                  if (AQUI / f).exists())

    copiadas, faltan, viejas = 0, [], []
    for archivo, _, _ in TODAS_LAS_FIGURAS:
        origen = AQUI / archivo
        if origen.exists():
            if origen.stat().st_mtime < _codigo:
                viejas.append(archivo)
            shutil.copy2(origen, SALIDA / "img" / archivo)
            copiadas += 1
        else:
            faltan.append(archivo)

    if viejas:
        print()
        print("  ABORTADO · FIGURAS MÁS VIEJAS QUE EL CÓDIGO QUE LAS DIBUJA")
        print("  Se contradirían con las tablas. Regenera con:  sh correr.sh web")
        for a in viejas:
            print(f"    · {a}")
        print()
        raise SystemExit(4)

    (SALIDA / "estilo.css").write_text(CSS, encoding="utf-8")

    print("  PÁGINAS")
    for archivo, titulo, subtitulo in PAGINAS + SUBPAGINAS_EXTRA:
        cuerpo = CUERPOS[archivo](d)
        (SALIDA / archivo).write_text(
            envoltura(archivo, titulo, subtitulo, cuerpo, fecha,
                      portada=(archivo == "index.html")),
            encoding="utf-8")
        print(f"    web/{archivo:20s} {titulo}")

    n_fichas = 0
    for archivo, titulo in FICHAS:
        p = RAIZ / "verificacion" / archivo
        if not p.exists():
            continue
        (SALIDA / "fichas" / f"{p.stem}.html").write_text(
            construir_ficha(p, titulo, fecha), encoding="utf-8")
        n_fichas += 1

    print(f"\n  web/estilo.css         hoja de estilo, sin dependencias externas")
    print(f"  web/fichas/            {n_fichas} fichas en HTML")
    print(f"  web/img/               {copiadas} figuras copiadas")
    if faltan:
        print("\n  FALTAN figuras (regenéralas con 'sh correr.sh todo' y "
              "'sh correr.sh dendrimero'):")
        for f in faltan:
            print(f"    {f}")
    print(f"\n  en {SALIDA}")
    print(f"  Ábrela con:  xdg-open '{SALIDA / 'index.html'}'")


if __name__ == "__main__":
    main()
