#!/usr/bin/env python3
# =============================================================================
#  construir_web.py · proyecto BHE / SENACYT
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
#  "MANTENIDO A MANO": la bibliografía, en sus dos formas: FUENTES (anclaje
#  por compuerta, lo que vigila la guarda de main()) y cuerpo_bibliografia()
#  (las referencias Vancouver que se publican).
#
#  USO:  python3 construir_web.py        (o: sh correr.sh web)
#  SALE: web/index.html, web/estilo.css y los PNG copiados.
# =============================================================================

import html
import re
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
#  MANTENIDO A MANO: lo que no vive en el código
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
#    1. SET COMÚN · las mismas tres figuras en todas: ventanas, matriz y
#       recorrido, generadas por `rutas.figuras()` con el catálogo de esa
#       página. Mismo orden y mismo pie siempre.
#    2. ESPECÍFICAS DE ESTA CLASE · lo que solo existe para una clase porque
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
# Registro de las 9 figuras generales. YA NO se copian a web/img (decisión de
# Jhovan, 2026-08-24) — se siguen generando en envolvimiento/ pero se quedan
# ahí. Las 3 envolvimiento_* son física general del modelo (radio crítico,
# barrera, G(D)) sin equivalente por liposoma: se conservan por si hacen falta
# para el informe SENACYT. Las 6 lip_* son las versiones combinadas de las
# formulaciones juntas; su contenido ya está desglosado por liposoma en
# web/img/liposoma_reales/ y web/img/liposoma_teoricos_detalle/. La sección
# general se retiró del HTML el 2026-08-18 (handoff 18m); esta constante
# queda solo como catálogo/documentación, ya no alimenta ninguna copia.
TODAS_LAS_FIGURAS = (FIGURAS_ENVOLVIMIENTO + FIGURAS_ENVOLVIMIENTO_TEORICOS
                     + FIGURAS_LIP_REALES + FIGURAS_LIP_TEORICOS)

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
    #  DETALLE POR DISEÑO: LIPOSOMA REALES Y TEÓRICOS (2026-08-18k). Mismo
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
.figura img{width:100%; display:block; background:#fff; cursor:zoom-in}
/* ---- lightbox de figuras ---- */
.lightbox{display:none; position:fixed; inset:0; background:rgba(20,22,26,.92);
  z-index:999; padding:3.5vh 3vw; box-sizing:border-box; overflow:auto}
.lightbox.abierto{display:flex; align-items:center; justify-content:center}
.lightbox img{max-width:100%; max-height:100%; border-radius:8px;
  box-shadow:0 10px 40px rgba(0,0,0,.5); cursor:zoom-in;
  transition:transform .15s ease; transform-origin:center center}
.lightbox img.zoom{max-width:none; max-height:none; width:180%; cursor:grab}
.lightbox img.zoom.arrastrando{cursor:grabbing; transition:none}
.lightbox.zoom-activo{justify-content:flex-start; align-items:flex-start}
.lightbox .cerrar{position:fixed; top:1.2rem; right:1.6rem; color:#fff;
  font-size:2rem; line-height:1; cursor:pointer; opacity:.85; z-index:1000}
.lightbox .cerrar:hover{opacity:1}
.lightbox .ayuda{position:fixed; bottom:1.2rem; left:50%; transform:translateX(-50%);
  color:#cfd3d8; font-size:.8rem; z-index:1000}

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

/* ---- citas y bibliografía ---- */
a.cita{text-decoration:none; font-weight:600; white-space:nowrap}
a.cita:hover{text-decoration:underline}
ol.biblio li:target{background:#fff8ec; border-radius:5px;
  box-shadow:0 0 0 .5rem #fff8ec}
ol.biblio{margin:1.2rem 0 0; padding-left:1.9rem; max-width:78ch}
ol.biblio li{margin-bottom:.85rem; font-size:.93rem; line-height:1.55}
ol.biblio a{font-size:.85rem; word-break:break-all}

code{background:#eef1f5; padding:.1rem .35rem; border-radius:4px; font-size:.87em;
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
a{color:var(--azul)}
.formula{text-align:center; font-size:1.02rem; margin:1.1rem 0}

footer{margin-top:3rem; padding:1.8rem 0 2.6rem; border-top:1px solid var(--linea);
  color:var(--suave); font-size:.85rem; background:var(--papel)}


@media (max-width:640px){
  header.principal h1{font-size:1.35rem} header.pagina h1{font-size:1.2rem}
}
@media print{nav.barra{display:none} body{background:#fff} .figura{break-inside:avoid}}
"""

# Las páginas del sitio. El orden es el de la navegación.
PAGINAS = [
    ("index.html", "Resumen", "Cifras y resultados."),
    ("liposoma.html", "Liposoma", "Suelo geométrico y liposomas reales publicados."),
    ("bibliografia.html", "Bibliografía",
     "Referencias citadas en el sitio, estilo Vancouver, en orden de aparición."),
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


# =============================================================================
#  CITAS VANCOUVER  ·  MANTENIDO A MANO junto a cuerpo_bibliografia()
# =============================================================================
#  El sitio cita con número entre corchetes, como anuncia la propia página de
#  bibliografía. Este mapa traduce el «Apellido AÑO» con que el código nombra
#  sus fuentes al número de la referencia. El número es el `value=` del <li>
#  correspondiente en cuerpo_bibliografia(): si se reordena la lista, hay que
#  reordenar esto.
#
#  Se cita por PRIMER APELLIDO + AÑO. Excepción declarada: el código llama
#  «Mishima y Kurano 2018» a la ref. 15, así que las dos formas apuntan al 15.
#
#  Regla dura: no se inventa un número. Si el sitio cita algo que no está
#  aquí, la guarda de main() aborta el build en vez de publicar un [n] mal
#  asignado o una cita sin referencia.
VANCOUVER = {
    "Mao 2014": 1,
    "Gong 2022": 2,
    "Chow 2025": 3,
    "Muselman 2026": 4,
    "Pan 2008": 5,
    "Weinbaum 2003": 6,
    "Deserno 2004": 7,
    "Bastiani 2010": 8,
    "Berry 2016": 9,
    "Mastorakos 2016": 10,
    "Tong 2016": 11,
    "Hisano 2011": 12,
    "Foster 2007": 13,
    "Bucki 2010": 14,
    "Mishima 2018": 15,
    "Kurano 2018": 15,
    "Nance 2012": 16,
    "Lockman 2004": 17,
    "Gromnicova 2016": 18,
    "Mouzoura 2025": 19,
    "Cheng 2016": 20,
    "Shi 2025": 21,
    "Larsen 2025": 22,
}

#  «Apellido [y/& Otro] [et al.] AÑO[, Revista vol:pág]». El fragmento de
#  revista es opcional y se descarta: con la cita numerada, repetir la revista
#  en el cuerpo del texto sobra, está en la referencia. Lo que va DESPUÉS
#  (Tabla 3, Fig. 1 A-E, Sec. III C…) se conserva: no es estilo de cita, es la
#  localización exacta dentro del paper, y sin ella se pierde trazabilidad.
#  El texto llega ESCAPADO, así que el «&» de «Bastiani & Parton» viaja como
#  «&amp;»: se contemplan las dos formas. Sin eso, el segundo autor se leería
#  como una cita propia («Parton 2010») que no existe en la bibliografía.
_RE_CITA = re.compile(
    r"\b([A-Z][a-zà-ÿ]{1,})"
    r"(?:\s+(?:y|&amp;|&)\s+[A-Z][a-zà-ÿ]+)?"
    r"(?:\s+et al\.?)?"
    r",?\s+((?:19|20)\d{2})"
    r"(?:\s*,\s*[^,·;]*?\d+\s*[:(]\s*[\w.]+[^,·;]*?(?=\s*(?:[,·;]|$)))?")

_SIN_REFERENCIA = set()


def _ref(n, base=""):
    return f'<a class="cita" href="{base}bibliografia.html#ref{n}">[{n}]</a>'


def _vancouver(texto, base=""):
    """Cambia «Apellido et al. AÑO, Revista vol:pág» por su [n] enlazado.

    El texto entra YA escapado. Lo que no esté en VANCOUVER se deja tal cual y
    se apunta en _SIN_REFERENCIA para que la guarda de main() aborte: antes
    publicar una cita sin convertir que un número inventado.
    """
    def _sub(m):
        clave = f"{m.group(1)} {m.group(2)}"
        n = VANCOUVER.get(clave)
        if n is None:
            _SIN_REFERENCIA.add(clave)
            return m.group(0)
        return _ref(n, base)
    return _RE_CITA.sub(_sub, texto)


def _cita_de(nombre):
    """El [n] del paper que da nombre a un diseño («Mao 2014 (real)» -> [1]).

    Devuelve "" si el nombre no lleva apellido+año, que es el caso de los
    teóricos y del dataset: esos no salen de ningún paper.
    """
    m = _RE_CITA.search(nombre)
    if not m:
        return ""
    clave = f"{m.group(1)} {m.group(2)}"
    n = VANCOUVER.get(clave)
    if n is None:
        _SIN_REFERENCIA.add(clave)
        return ""
    return " " + _ref(n)


def _etq(v):
    clase = {"EXCLUIDA": "no", "NO EXCLUIDA": "si", "NO EVALUABLE": "nn"}[v]
    texto = {"EXCLUIDA": "NO", "NO EXCLUIDA": "SÍ", "NO EVALUABLE": "??"}[v]
    return f'<span class="etq {clase}">{texto}</span>'


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
{LIGHTBOX}
</body></html>
"""


# =============================================================================
#  LIGHTBOX DE FIGURAS
# =============================================================================
#  Clic en una figura: se abre a pantalla completa. Clic en la imagen abierta:
#  zoom al 180 %, y entonces se puede arrastrar con el ratón o el dedo. Clic
#  fuera, la X, o Escape: cerrar.
#
#  Recuperado literal del commit 2efce14 (`web/dataset_50.html`), donde vivía
#  solo en la página del dataset. Se perdió al pasar esa página a
#  construir_web.py. Ahora va en envoltura(), así que lo tienen TODAS las
#  páginas: cualquier `.figura img` del sitio se amplía.
#
#  Va en una constante aparte y NO en el f-string de envoltura() porque las
#  llaves del JavaScript chocan con las de la interpolación.
LIGHTBOX = """<div class="lightbox" id="lightbox"><span class="cerrar">&times;</span><img id="lightbox-img" src="" alt=""><span class="ayuda">clic en la imagen: zoom · clic afuera o Escape: cerrar</span></div>
<script>
(function(){
  var lb = document.getElementById('lightbox');
  var lbImg = document.getElementById('lightbox-img');
  document.querySelectorAll('.figura img').forEach(function(img){
    img.addEventListener('click', function(){
      lbImg.src = img.src;
      lbImg.alt = img.alt;
      lbImg.classList.remove('zoom');
      lb.classList.remove('zoom-activo');
      lb.classList.add('abierto');
    });
  });
  var arrastrando = false, movio = false, x0 = 0, y0 = 0, sl0 = 0, st0 = 0;
  function empezarArrastre(x, y){
    arrastrando = true; movio = false;
    x0 = x; y0 = y; sl0 = lb.scrollLeft; st0 = lb.scrollTop;
    lbImg.classList.add('arrastrando');
  }
  function moverArrastre(x, y){
    if (!arrastrando) return;
    if (Math.abs(x - x0) > 3 || Math.abs(y - y0) > 3) movio = true;
    lb.scrollLeft = sl0 - (x - x0);
    lb.scrollTop = st0 - (y - y0);
  }
  function terminarArrastre(){ arrastrando = false; lbImg.classList.remove('arrastrando'); }
  lbImg.addEventListener('mousedown', function(e){
    if (!lbImg.classList.contains('zoom')) return;
    e.preventDefault();
    empezarArrastre(e.clientX, e.clientY);
  });
  document.addEventListener('mousemove', function(e){ moverArrastre(e.clientX, e.clientY); });
  document.addEventListener('mouseup', terminarArrastre);
  lbImg.addEventListener('touchstart', function(e){
    if (!lbImg.classList.contains('zoom')) return;
    var t = e.touches[0]; empezarArrastre(t.clientX, t.clientY);
  }, {passive: true});
  lbImg.addEventListener('touchmove', function(e){
    var t = e.touches[0]; moverArrastre(t.clientX, t.clientY);
  }, {passive: true});
  lbImg.addEventListener('touchend', terminarArrastre);
  lbImg.addEventListener('click', function(e){
    e.stopPropagation();
    if (movio) { movio = false; return; }
    lbImg.classList.toggle('zoom');
    lb.classList.toggle('zoom-activo', lbImg.classList.contains('zoom'));
    if (lbImg.classList.contains('zoom')) { lb.scrollTop = 0; lb.scrollLeft = 0; }
  });
  function cerrar(){
    lb.classList.remove('abierto', 'zoom-activo');
    lbImg.classList.remove('zoom');
    lbImg.src = '';
  }
  lb.addEventListener('click', cerrar);
  document.addEventListener('keydown', function(e){ if (e.key === 'Escape') cerrar(); });
})();
</script>"""


# =============================================================================
#  CUERPO DE CADA PÁGINA
# =============================================================================

def cuerpo_index(d):
    u, lip = d["umbrales"], d["liposoma"]
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
</div>

<h2>Resultados</h2>
<p>Los dos rangos de tamaño teóricos para el paso pasivo hacia el cerebro no se solapan: el tamiz del
glicocálix exige un diámetro de hasta {u['glicocalix']:.2f} nm y el envolvimiento de membrana exige al
menos {u['envolvimiento']:.2f} nm. La separación entre ambos rangos es de {d['margen_ventanas']:.2f} nm,
un resultado geométrico firme dentro del barrido de sensibilidad aplicado. Sin embargo, superar el
tamaño del tamiz ya no se cuenta como exclusión automática: dos medidas reales en barrera
hematoencefálica (una en tejido nativo, otra en cultivo) muestran nanopartículas más grandes que ese
poro cruzando de forma medible, así que el simulador reporta ese caso como dato insuficiente para
decidir, no como bloqueo, salvo que otro criterio (carga superficial, tiempos de tránsito) sí decida
por su cuenta.</p>
<p>De las formulaciones reales evaluadas:</p>
{_resumen_catalogo_prosa(d["catalogo_real"])}
<p><b>Empates técnicos.</b> Uno de estos resultados queda dentro del margen de duda del propio
barrido y no debe citarse como cerrado: la separación entre los dos rangos de tamaño cae a 0.26 nm
en el escenario más permisivo del barrido (κ = 15 kT, Hamaker = 6.5·10⁻²¹ J).</p>

<h2>Alcance y limitaciones del simulador</h2>
<p>Este simulador es una herramienta de cribado con trazabilidad completa a fuente científica
publicada, no un modelo predictivo validado. La diferencia importa: cribado significa que cada
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
no existe un parámetro publicado que conecte el LogP de un nanotransportador completo (como un
liposoma) con su capacidad de adherirse o atravesar la barrera. Añadirlo sin una fuente científica
que lo respalde habría significado inventar un mecanismo, algo que este proyecto decidió no hacer.</p>
<p>El modelo cubre únicamente las etapas de adhesión y envolvimiento del nanotransportador en la
superficie del vaso sanguíneo cerebral. Las etapas posteriores (el paso activo hacia el otro lado de
la barrera, el tránsito celular y la liberación final del fármaco) dependen de energía celular y
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
        nombre = html.escape(f["nombre"]) + _cita_de(f["nombre"])
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
<p class="rev">Fuente: {_ref(VANCOUVER["Pan 2008"])}, Fig. 3c
(t bicapa medido, {lip['t_bicapa'][0]:.1f}–{lip['t_bicapa'][-1]:.1f} nm).</p>

<h2>Liposomas publicados</h2>
<p class="rev">Ø y ζ medidos, con su referencia numerada junto al nombre.</p>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Primera compuerta que falla</h2>
{tabla_muere}


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
        cita = _cita_de(f["nombre"])
        nota = ("" if cita else
                (f'<div class="rev">{html.escape(f["nota"])}</div>' if f["nota"] else ""))
        filas_v.append(f'<tr><td><b>{html.escape(f["nombre"])}</b>{cita}{nota}</td>'
                       f'<td class="num">{f["diametro"]:.1f}</td>'
                       f'<td class="num">{f["zeta"]:+.2f}</td>'
                       f'<td class="num">{f["peg"]:.0f}</td>' + "".join(cv) + '</tr>')
        filas_m.append(f'<tr><td><b>{html.escape(f["nombre"])}</b>{cita}{nota}</td>'
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

<h2>Veredictos</h2>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Primera compuerta que falla</h2>
{tabla_muere}




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
                f'<td class="rev">{_vancouver(html.escape(c["fuente"] or "—"))}</td></tr>'
                for c in info["compuertas"])
            secciones_rutas.append(
                f'<h4>{html.escape(nombre_ruta)} · {_etq(info["veredicto"])}</h4>'
                '<div class="tabla-scroll"><table><thead><tr><th>Compuerta</th>'
                '<th>Estado</th><th class="num">Valor</th><th class="num">Umbral</th>'
                f'<th>Unidad</th><th>Fuente</th></tr></thead><tbody>{filas_c}</tbody></table></div>')
        fab = it.get("fabricable")
        insignia_fab = (f' · <span class="etq {"si" if fab else "no"}">'
                        f'Fabricable: {"SÍ" if fab else "NO"}</span>'
                        if fab is not None else "")
        cita = _cita_de(it["nombre"])
        nota = ("" if cita else
                (f' <span class="rev">({html.escape(it["nota"])})</span>'
                 if it.get("nota") else ""))
        bloques.append(f'''
<details class="liposoma">
<summary>{html.escape(it["nombre"])}{cita}{nota} · Ø {it["diametro"]:.1f} nm ·
ζ {it["zeta"]:+.2f} mV · PEG {it["peg"]:.2f} nm{insignia_fab}</summary>
<div class="dataset50-figs">{figs}</div>
{"".join(secciones_rutas)}
</details>''')
    return "".join(bloques)


def cuerpo_dataset_50(d):
    """50 liposomas sintéticos (semilla 42), último día 2026-08-18. Rango de
    diámetro/ζ/PEG tomado de CATALOGO_REAL + teóricos. NO citables como
    predicción; mismo criterio que los diseños teóricos. La carga útil (G.2)
    queda DESCONOCIDA en los 50, por decisión vigente del 2026-08-13."""
    tabla_veredictos, tabla_muere = _tablas_catalogo(d["dataset_50"], d["rutas"])

    return f"""
<div class="aviso"><strong>Dataset de 50 diseños SINTÉTICOS.</strong>
No son mediciones ni predicciones: no deben citarse.</div>

<h2>Veredictos</h2>
{tabla_veredictos}
<p class="rev">Cada <b>!</b> junto a un SÍ = una compuerta que pasa con salvedad
declarada.</p>

<h2>Detalle por liposoma</h2>
<p class="rev">Clic sobre cada fila para desplegar sus 3 figuras propias
(ventanas, matriz, recorrido) y el desglose de compuertas por ruta.</p>
{_bloques_detalle_liposoma(d["dataset_50_detalle"])}
"""


def cuerpo_bibliografia(d):
    """BIBLIOGRAFÍA · VA A MANO, igual que FUENTES.

    Restaurada literal del último HTML publicado que la tenía
    (commit 2efce14, `web/bibliografia.html`, generado el 2026-08-18). Se
    recuperó tal cual: 22 referencias en estilo Vancouver, en orden de
    aparición, con su numeración original intacta. NO se reconstruyó de
    memoria ni se completó a ojo.

    NO confundir con `FUENTES`, que es otra lista y tiene más entradas: FUENTES
    es el anclaje interno del código (una por compuerta, incluidos manuales de
    coloides y fuentes que solo viven en comentarios) y es lo que vigila la
    guarda de `main()`. Esta lista son las referencias CITADAS EN EL SITIO, que
    son menos por definición. Que FUENTES tenga más entradas no es un hueco.

    Verificado el 2026-08-22: las 4 páginas del sitio citan 14 apellido+año, y
    las 14 tienen su entrada aquí.
    """
    return """<h2>Referencias</h2>
<p class="rev">Las citas en el texto usan el número entre corchetes que corresponde a esta lista,
p. ej. [5]. Cuando una nota cita más de una fuente aparecen varios números juntos, p. ej. [1,11]
o [13-15].</p>
<ol class="biblio">
<li id="ref1" value="1">Mao Y, Wang J, Zhao Y, Wu Y, Kwak KJ, Chen CS, et al. A novel liposomal formulation of FTY720 (fingolimod) for promising enhanced targeted delivery. Nanomedicine. 2014;10(2):393-400. <a href="https://doi.org/10.1016/j.nano.2013.08.001" target="_blank" rel="noopener">doi:10.1016/j.nano.2013.08.001</a></li>
<li id="ref2" value="2">Gong X, Fan X, He Y, Wang Y, Zhou F, Yang B. A pH-sensitive liposomal co-delivery of fingolimod and ammonia borane for treatment of intracerebral hemorrhage. Nanophotonics. 2022;11(22):5133-42. <a href="https://doi.org/10.1515/nanoph-2022-0496" target="_blank" rel="noopener">doi:10.1515/nanoph-2022-0496</a></li>
<li id="ref3" value="3">Chow SF, et al. Rational development of fingolimod nano-embedded microparticles as nose-to-brain neuroprotective therapy for ischemic stroke. Drug Deliv Transl Res. 2025;15(6):2022-47. <a href="https://doi.org/10.1007/s13346-024-01721-8" target="_blank" rel="noopener">doi:10.1007/s13346-024-01721-8</a></li>
<li id="ref4" value="4">Muselman A, Yu LW, Nguyen KD, Inayathullah M, Liu Q, Brewer KD, et al. Macrophage-targeted PEGylated liposomes ameliorate experimental autoimmune encephalomyelitis. Front Immunol. 2026;16:1657131. <a href="https://doi.org/10.3389/fimmu.2025.1657131" target="_blank" rel="noopener">doi:10.3389/fimmu.2025.1657131</a></li>
<li id="ref5" value="5">Pan J, Tristram-Nagle S, Kučerka N, Nagle JF. Temperature dependence of structure, bending rigidity, and bilayer interactions of DOPC bilayers. Phys Rev Lett. 2008;100:198103. <a href="https://doi.org/10.1103/PhysRevLett.100.198103" target="_blank" rel="noopener">doi:10.1103/PhysRevLett.100.198103</a></li>
<li id="ref6" value="6">Weinbaum S, Zhang X, Han Y, Vink H, Cowin SC. Mechanotransduction and flow across the endothelial glycocalyx. Proc Natl Acad Sci U S A. 2003;100(13):7988-95. <a href="https://doi.org/10.1073/pnas.1332808100" target="_blank" rel="noopener">doi:10.1073/pnas.1332808100</a></li>
<li id="ref7" value="7">Deserno M. Elastic deformation of a fluid membrane upon colloid binding. Phys Rev E. 2004;69:031903. <a href="https://doi.org/10.1103/PhysRevE.69.031903" target="_blank" rel="noopener">doi:10.1103/PhysRevE.69.031903</a></li>
<li id="ref8" value="8">Bastiani M, Parton RG. Caveolae at a glance. J Cell Sci. 2010;123(22):3831-6. <a href="https://doi.org/10.1242/jcs.070102" target="_blank" rel="noopener">doi:10.1242/jcs.070102</a></li>
<li id="ref9" value="9">Berry S, Mastorakos P, Zhang C, Song E, Patel H, Suk JS, et al. Enhancing intracranial delivery of clinically relevant non-viral gene vectors. RSC Adv. 2016;6:41665-74. <a href="https://doi.org/10.1039/c6ra01546h" target="_blank" rel="noopener">doi:10.1039/c6ra01546h</a></li>
<li id="ref10" value="10">Mastorakos P, Song E, Zhang C, Berry S, Park HW, Kim YE, et al. Biodegradable DNA nanoparticles that provide widespread gene delivery in the brain. Small. 2016;12(5):678-85. <a href="https://doi.org/10.1002/smll.201502554" target="_blank" rel="noopener">doi:10.1002/smll.201502554</a></li>
<li id="ref11" value="11">Tong H-I, Kang W, Davy PMC, Shi Y, Sun S, Allsopp RC, et al. Monocyte trafficking, engraftment, and delivery of nanoparticles and an exogenous gene into the acutely inflamed brain tissue. PLoS One. 2016;11(4):e0154022. <a href="https://doi.org/10.1371/journal.pone.0154022" target="_blank" rel="noopener">doi:10.1371/journal.pone.0154022</a></li>
<li id="ref12" value="12">Hisano Y, Kobayashi N, Kawahara A, Yamaguchi A, Nishi T. The sphingosine 1-phosphate transporter, SPNS2, functions as a transporter of the phosphorylated form of the immunomodulating agent FTY720. J Biol Chem. 2011;286(3):1758-66. <a href="https://doi.org/10.1074/jbc.M110.171116" target="_blank" rel="noopener">doi:10.1074/jbc.M110.171116</a></li>
<li id="ref13" value="13">Foster CA, Howard LM, Schweitzer A, Persohn E, Hiestand PC, Balatoni B, et al. Brain penetration of the oral immunomodulatory drug FTY720 and its phosphorylation in the central nervous system during experimental autoimmune encephalomyelitis: consequences for mode of action in multiple sclerosis. J Pharmacol Exp Ther. 2007;323(2):469-75. <a href="https://doi.org/10.1124/jpet.107.127183" target="_blank" rel="noopener">doi:10.1124/jpet.107.127183</a></li>
<li id="ref14" value="14">Bucki R, Kulakowska A, Byfield FJ, Zendzian-Piotrowska M, Baranowski M, Marzec M, et al. Plasma gelsolin modulates cellular response to sphingosine 1-phosphate. Am J Physiol Cell Physiol. 2010;299(6):C1516-23. <a href="https://doi.org/10.1152/ajpcell.00051.2010" target="_blank" rel="noopener">doi:10.1152/ajpcell.00051.2010</a></li>
<li id="ref15" value="15">Mishima Y, Kurano M, Kobayashi T, Nishikawa M, Ohkawa R, Tozuka M, et al. Dihydro-sphingosine 1-phosphate interacts with carrier proteins in a manner distinct from that of sphingosine 1-phosphate. Biosci Rep. 2018;38(5):BSR20181288. <a href="https://doi.org/10.1042/BSR20181288" target="_blank" rel="noopener">doi:10.1042/BSR20181288</a></li>
<li id="ref16" value="16">Nance EA, Woodworth GF, Sailor KA, Shih TY, Xu Q, Swaminathan G, et al. A dense poly(ethylene glycol) coating improves penetration of large polymeric nanoparticles within brain tissue. Sci Transl Med. 2012;4(149):149ra119. <a href="https://doi.org/10.1126/scitranslmed.3003594" target="_blank" rel="noopener">doi:10.1126/scitranslmed.3003594</a></li>
<li id="ref17" value="17">Lockman PR, Koziara JM, Mumper RJ, Allen DD. Nanoparticle surface charges alter blood-brain barrier integrity and permeability. J Drug Target. 2004;12(9-10):635-41. <a href="https://doi.org/10.1080/10611860400015936" target="_blank" rel="noopener">doi:10.1080/10611860400015936</a></li>
<li id="ref18" value="18">Gromnicova R, Kaya M, Romero IA, Williams P, Satchell S, Sharrack B, et al. Transport of gold nanoparticles by vascular endothelium from different human tissues. PLoS One. 2016;11(8):e0161610. <a href="https://doi.org/10.1371/journal.pone.0161610" target="_blank" rel="noopener">doi:10.1371/journal.pone.0161610</a></li>
<li id="ref19" value="19">Mouzoura P, Marazioti A, Gkartziou F, Metsiou D-N, Antimisiaris SG. Potential of liposomal FTY720 for bone regeneration: proliferative, osteoinductive, chemoattractive, and angiogenic properties compared to free bioactive lipid. Int J Nanomedicine. 2025;20:239-65. <a href="https://doi.org/10.2147/IJN.S494512" target="_blank" rel="noopener">doi:10.2147/IJN.S494512</a></li>
<li id="ref20" value="20">Cheng MJ, Kumar R, Sridhar S, Webster TJ, Ebong EE. Endothelial glycocalyx conditions influence nanoparticle uptake for passive targeting. Int J Nanomedicine. 2016;11:3305-15. <a href="https://doi.org/10.2147/IJN.S106299" target="_blank" rel="noopener">doi:10.2147/IJN.S106299</a></li>
<li id="ref21" value="21">Shi SM, Suh RJ, Shon DJ, Garcia FJ, Buff JK, Atkins M, et al. Glycocalyx dysregulation impairs blood-brain barrier in ageing and disease. Nature. 2025;639:985-94. <a href="https://doi.org/10.1038/s41586-025-08589-9" target="_blank" rel="noopener">doi:10.1038/s41586-025-08589-9</a></li>
<li id="ref22" value="22">Larsen R, Kucharz K, Aydin S, Micael MKB, Choudhury B, Paulchakrabarti M, et al. Multi-omic analysis reveals the unique glycan landscape of the blood-brain barrier glycocalyx [preprint]. bioRxiv. 2025:2025.04.07.645297. <a href="https://doi.org/10.1101/2025.04.07.645297" target="_blank" rel="noopener">doi:10.1101/2025.04.07.645297</a></li>
</ol>
"""


CUERPOS = {
    "index.html": cuerpo_index,
    "liposoma.html": cuerpo_liposoma,
    "liposoma_teoricos.html": cuerpo_liposoma_teoricos,
    "dataset_50.html": cuerpo_dataset_50,
    "bibliografia.html": cuerpo_bibliografia,
}


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

    # TODAS_LAS_FIGURAS (envolvimiento_*, lip_reales_*, lip_teoricos_*) ya NO
    # se copian aquí a web/img — ver comentario junto a esa constante. La
    # guarda contra figuras viejas que existía en este punto (Bug detectado
    # por Jhovan, hasta el 2026-08-12 la web mezclaba PNG de momentos
    # distintos) ya no aplica: nada de disco se copia a ciegas en este build.
    # Las figuras que SÍ van a la web (dataset_50, liposoma_reales,
    # liposoma_teoricos_detalle) se generan frescas en cada build, más arriba
    # en este mismo script, así que no pueden quedar desactualizadas.

    (SALIDA / "estilo.css").write_text(CSS, encoding="utf-8")

    print("  PÁGINAS")
    for archivo, titulo, subtitulo in PAGINAS + SUBPAGINAS_EXTRA:
        cuerpo = CUERPOS[archivo](d)
        (SALIDA / archivo).write_text(
            envoltura(archivo, titulo, subtitulo, cuerpo, fecha,
                      portada=(archivo == "index.html")),
            encoding="utf-8")
        print(f"    web/{archivo:20s} {titulo}")

    # GUARDA CONTRA CITA SIN REFERENCIA. Va DESPUÉS de escribir las páginas
    # porque _vancouver() solo se entera de lo que de verdad se publicó. Si el
    # sitio nombra un «Apellido AÑO» que no está en VANCOUVER, la cita se queda
    # sin número: no se publica así, se aborta. Nunca un [n] inventado.
    if _SIN_REFERENCIA:
        print()
        print("  ABORTADO · CITAS DEL SITIO SIN NÚMERO DE REFERENCIA")
        for c in sorted(_SIN_REFERENCIA):
            print(f"    · {c}")
        print("  Añádelas a cuerpo_bibliografia() y a VANCOUVER en construir_web.py")
        print()
        raise SystemExit(6)

    _n_img = sum(1 for _ in (SALIDA / "img").rglob("*.png"))
    print(f"\n  web/estilo.css         hoja de estilo, sin dependencias externas")
    print(f"  web/img/               {_n_img} figuras (detalle por diseño; las 9 generales")
    print(f"                         de envolvimiento_*/lip_reales_*/lip_teoricos_* ya no")
    print(f"                         se copian aquí, quedan en envolvimiento/)")
    print(f"\n  en {SALIDA}")
    print(f"  Ábrela con:  xdg-open '{SALIDA / 'index.html'}'")


if __name__ == "__main__":
    main()
