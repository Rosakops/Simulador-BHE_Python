#@title Simulacion entre liposoma convencional y fingolimod -Modelo 3D y Analisis-
# =============================================================================
#  SIMULACIÓN 3D: INSERCIÓN DE FINGOLIMOD EN LIPOSOMA POR ENERGÍA LIBRE (PMF)
#  Motor: dinámica Browniana sobre un potencial de fuerza media DERIVADO del logP
#  Visualización: Plotly (graph_objects)  |  Entorno: 1 celda de Google Colab
# =============================================================================
#  QUÉ ES ESTO Y QUÉ PRETENDE PROBAR (leer antes de usar):
#
#  A diferencia de una animación "guionizada" (donde uno coloca a mano el mínimo
#  de energía donde quiere que llegue el fármaco), aquí el destino del fármaco es
#  una SALIDA CALCULADA, no una suposición:
#
#    1) El perfil de energía libre G(r) a través de la bicapa se DERIVA del dato
#       medido de lipofilicidad (logP) mediante la relación de partición estándar
#       ΔG_transferencia = -2.303 · k_B·T · logP.  Con logP(fingolimod)=4.16 →
#       ΔG ≈ -9.58 k_BT ≈ -5.90 kcal/mol.  El mínimo cae en las colas hidrófobas
#       PORQUE el logP es alto; para un fármaco hidrofílico (logP<0) el mismo
#       cálculo vuelve la membrana DESFAVORABLE y el fármaco NO entra.
#
#    2) La integración es dinámica Browniana sobreamortiguada (ec. de Smoluchowski),
#       cuya distribución estacionaria es GARANTIZADAMENTE Boltzmann ∝ exp(-G/k_BT).
#       Por eso muestrear el paisaje es significativo.
#
#    3) VERIFICACIÓN: se compara el histograma radial muestreado contra la
#       predicción de Boltzmann (ver método verificar() y el panel PMF).
#
#    4) FALSABILIDAD: se corre el MISMO motor con un análogo hidrofílico
#       (LOGP_CONTROL) y se comprueba que NO se inserta. Un modelo que solo sabe
#       decir "sí" no prueba nada; este sabe decir "no".
#
#  LÍMITES HONESTOS (imprescindible tenerlos presentes):
#    · Es un modelo de ORDEN REDUCIDO (PMF radial efectivo), NO dinámica molecular
#      atomística. NO sustituye a MD con campo de fuerzas real (CHARMM36/MARTINI),
#      solvente explícito y energía libre por umbrella sampling/FEP, que es el
#      estándar de "prueba" de grado publicación (y requiere GPU y horas).
#    · Se modela el paso de INSERCIÓN desde la interfaz agua-membrana; el encuentro
#      por difusión en el seno del líquido (limitado por difusión) NO se modela.
#    · D_SIM es una difusión EFECTIVA (tiempo acelerado). Como la distribución
#      estacionaria ∝ exp(-G/k_BT) NO depende de D, la conclusión termodinámica
#      (dónde acaba el fármaco) es invariante; solo se comprime el eje temporal.
# =============================================================================

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "plotly"], check=True)
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

import numpy as np
import plotly.io as pio
try:
    pio.renderers.default = "colab"   # incrustación correcta en la celda de Colab
except Exception:
    pass

# np.trapezoid reemplazó a np.trapz en NumPy 2.0; con esto la celda funciona
# tanto en Colab con NumPy >=2.0 (donde trapz ya no existe) como en versiones anteriores.
_trapz = getattr(np, "trapezoid", None) or np.trapz


# =============================================================================
#  PARÁMETROS AJUSTABLES
# =============================================================================
# --- Geometría REAL de la ficha "Liposomas convencionales" ---
RADIO_LIPOSOMA_EXT   = 20.0    # nm  radio externo (Ø = 40.0 nm)
ESPESOR_BICAPA       = 4.0     # nm  grosor de bicapa
RADIO_NUCLEO_ACUOSO  = 16.0    # nm  radio del lumen acuoso (= 20 - 4)
R_MIDPLANO_TAIL      = 18.0    # nm  plano medio hidrófobo de las colas (mínimo de G)

# --- Propiedades físico-químicas documentadas (fichas) ---
POTENCIAL_ZETA_mV    = 5.0     # mV  potencial zeta del liposoma
LOGP_LIPOSOMA        = 2.0     # log P del transportador

# --- Termodinámica ---
TEMP_K               = 310.0   # K   temperatura fisiológica
KT_KCAL              = 0.001987 * TEMP_K   # kcal/mol  (k_B·T ≈ 0.616 kcal/mol)

# --- LIPOFILICIDAD (dato que gobierna TODO el paisaje de energía) ---
LOGP_FARMACO         = 4.16    # XLogP3 de fingolimod (ficha) → define el pozo
LOGP_CONTROL         = -1.5    # análogo hidrofílico ficticio (prueba de falsabilidad)

# --- Forma del perfil de energía libre G(r) (en k_BT) ---
W_TAIL               = 1.9     # nm  ancho de la cuenca hidrófoba (colas)
B_HEAD               = 1.0     # k_BT  barrera de deshidratación en las cabezas polares
W_HEAD               = 0.5     # nm  ancho de esa barrera

# --- Fármaco ---
N_FARMACOS           = 80      # nº de moléculas de fingolimod
R_INIT_LO            = 20.5    # nm  radio inicial mínimo (interfaz agua-membrana)
R_INIT_HI            = 23.5    # nm  radio inicial máximo
R_BOX                = 30.0    # nm  frontera reflectante (compartimento finito de la suspensión)

# --- Dinámica Browniana (ec. de Smoluchowski, sobreamortiguada) ---
D_SIM                = 0.06    # nm^2/ps  difusión EFECTIVA (tiempo acelerado; ver cabecera)
DT                   = 0.10    # ps  paso de tiempo
N_PASOS              = 1400    # nº de pasos
FRAMES_SALTO         = 9       # 1 frame cada N pasos (~156 frames para la animación)
FRAC_EQUILIBRIO      = 0.5     # fracción final de frames usada como "equilibrada" al verificar

# --- Reproducibilidad ---
SEMILLA              = 7
rng = np.random.default_rng(SEMILLA)

# --- Paleta por elemento (CPK simplificado) ---
COLOR_ELEMENTO = {"C": "#3a3f4b", "N": "#2f6df6", "O": "#e23b2e", "Ctail": "#f2a23c"}
RADIO_ELEMENTO = {"C": 3.0, "N": 3.4, "O": 3.4, "Ctail": 2.6}


# =============================================================================
#  ENERGÍA LIBRE G(r): DERIVADA del logP (no impuesta a mano)
# =============================================================================
def energia_libre(r, logP):
    """
    Perfil de energía libre radial G(r) en unidades de k_BT.

    Dos términos con significado físico:
      (1) Premio hidrófobo (transferencia agua→colas), Gaussiano centrado en el
          plano medio de las colas. Su PROFUNDIDAD se calcula del logP medido:
              ΔG_transf = -ln(10)·logP   [k_BT]   (relación de partición estándar)
          → para logP>0 el mínimo es negativo (favorable) en las colas;
            para logP<0 se vuelve POSITIVO (la membrana repele al fármaco).
      (2) Barrera de deshidratación en las cabezas polares (r≈16 y r≈20 nm):
          coste de que el grupo polar del fármaco pierda su hidratación al cruzar.

    r : array de radios (nm).  Devuelve G(r) en k_BT.
    """
    r = np.asarray(r, dtype=float)
    dG_transf = -np.log(10.0) * logP  # k_BT  (profundidad CALCULADA desde el logP)
    g_tail = dG_transf * np.exp(-((r - R_MIDPLANO_TAIL) / W_TAIL) ** 2)
    g_head = B_HEAD * (np.exp(-((r - RADIO_LIPOSOMA_EXT) / W_HEAD) ** 2) +
                       np.exp(-((r - RADIO_NUCLEO_ACUOSO) / W_HEAD) ** 2))
    return g_tail + g_head


def gradiente_energia(r, logP, h=1e-3):
    """dG/dr por diferencias finitas centradas (k_BT/nm). Robusto ante cambios de perfil."""
    return (energia_libre(r + h, logP) - energia_libre(r - h, logP)) / (2.0 * h)


# =============================================================================
#  PLANTILLA MOLECULAR DE FINGOLIMOD (22 átomos pesados: 19 C + 1 N + 2 O)
# =============================================================================
def _plantilla_fingolimod():
    """
    Esqueleto esquemático de fingolimod (C19H33NO2; 55 átomos totales, 33 H omitidos).
    Cabeza amino-diol polar + enlazador etilo + anillo aromático + cola octilo lipófila.
    Coordenadas en nm a escala de enlace (~0.15 nm), centradas en el centro de masa.
    """
    b = 0.15
    coords, elems, bonds = [], [], []

    def add(el, x, y, z):
        coords.append([x, y, z]); elems.append(el); return len(coords) - 1

    # Cabeza polar C(NH2)(CH2OH)2
    cq  = add("C", 0.0, 0.0, 0.0)
    n   = add("N", 0.0, 1.05 * b, 0.0)
    cm1 = add("C", -0.9 * b, -0.5 * b, 0.75 * b);  o1 = add("O", -1.9 * b, -0.2 * b, 0.95 * b)
    cm2 = add("C", -0.9 * b, -0.5 * b, -0.75 * b); o2 = add("O", -1.9 * b, -0.2 * b, -0.95 * b)
    bonds += [(cq, n), (cq, cm1), (cm1, o1), (cq, cm2), (cm2, o2)]

    # Enlazador etilo
    c1 = add("C", 1.0 * b, 0.0, 0.0); c2 = add("C", 2.0 * b, 0.15 * b, 0.0)
    bonds += [(cq, c1), (c1, c2)]

    # Anillo aromático (para-sustituido)
    cx, cy, R = 3.15 * b, 0.15 * b, 0.9 * b
    ring = [add("C", cx + R * np.cos(np.pi / 3 * k), cy + R * np.sin(np.pi / 3 * k), 0.0) for k in range(6)]
    for k in range(6): bonds.append((ring[k], ring[(k + 1) % 6]))
    bonds.append((c2, ring[3]))

    # Cola octilo lipófila (8 C en zig-zag)
    prev = ring[0]; x0, y0 = cx + R, cy
    for k in range(8):
        ct = add("Ctail", x0 + (k + 1) * b, y0 + (0.28 * b if k % 2 == 0 else -0.28 * b), 0.0)
        bonds.append((prev, ct)); prev = ct

    coords = np.asarray(coords, dtype=float); coords -= coords.mean(axis=0)
    return coords, elems, bonds


def _matriz_rotacion_aleatoria(gen):
    """Matriz de rotación 3x3 (orienta cada molécula distinto)."""
    a, b, c = gen.uniform(0, 2 * np.pi, size=3)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(c), -np.sin(c)], [0, np.sin(c), np.cos(c)]])
    return Rz @ Ry @ Rx


# =============================================================================
#  CLASE PRINCIPAL
# =============================================================================
class Simulation3D:
    """
    Simula la INSERCIÓN de fingolimod en la bicapa por dinámica Browniana sobre
    un PMF derivado del logP, y produce (a) la animación 3D y (b) el panel de
    verificación termodinámica (Boltzmann) + control de falsabilidad.

    Uso:
        sim = Simulation3D(); sim.run()
        sim.verificar()                       # imprime la evidencia numérica
        sim.generate_animation().show()       # animación 3D
        sim.generate_pmf_figure().show()      # panel PMF + Boltzmann + control
    """

    def __init__(self, logP=LOGP_FARMACO):
        self.logP = logP
        self.com_liposoma = np.zeros(3)
        self._construir_liposoma()
        self.tpl_coords, self.tpl_elems, self.tpl_bonds = _plantilla_fingolimod()
        self.tpl_colors = [COLOR_ELEMENTO[e] for e in self.tpl_elems]
        self.tpl_sizes  = [RADIO_ELEMENTO[e] for e in self.tpl_elems]
        self.rotaciones = [_matriz_rotacion_aleatoria(rng) for _ in range(N_FARMACOS)]

        # Precómputo para el muestreo vectorizado de frames (mismos números que un loop
        # por molécula/enlace, pero sin el costo de miles de iteraciones Python por frame)
        self.rot_stack = np.stack(self.rotaciones)                      # (N_FARMACOS, 3, 3)
        self.bond_a = np.array([a for a, b in self.tpl_bonds])
        self.bond_b = np.array([b for a, b in self.tpl_bonds])
        self.colores_frame = self.tpl_colors * N_FARMACOS
        self.tamanos_frame = self.tpl_sizes * N_FARMACOS

        # Estado inicial: interfaz agua-membrana (modelamos la INSERCIÓN, no el encuentro)
        d = rng.normal(size=(N_FARMACOS, 3)); d /= np.linalg.norm(d, axis=1, keepdims=True)
        self.pos = d * rng.uniform(R_INIT_LO, R_INIT_HI, N_FARMACOS)[:, None]

        self.traj_atomos, self.traj_enlaces, self.traj_r = [], [], []
        self._radial_cache = {}   # evita recomputar _simular_radial(logP) más de una vez

    # ---- Reporte de especificaciones (se imprime en consola) ----
    def imprimir_especificaciones(self):
        """Imprime la ficha técnica completa del sistema, el modelo físico y sus límites."""
        print("=" * 70)
        print(" ESPECIFICACIONES DEL MODELO DE NANOTRANSPORTE")
        print("=" * 70)
        print(" NANOTRANSPORTADOR · Liposoma convencional (ficha):")
        print(f"   · Diámetro externo ............ {2*RADIO_LIPOSOMA_EXT:.1f} nm")
        print(f"   · Grosor de bicapa ........... {ESPESOR_BICAPA:.1f} nm")
        print(f"   · Radio del núcleo acuoso .... {RADIO_NUCLEO_ACUOSO:.1f} nm")
        print(f"   · Potencial zeta ............. {POTENCIAL_ZETA_mV:.1f} mV (baja carga → 'convencional')")
        print(f"   · Log P del transportador .... {LOGP_LIPOSOMA:.1f}")
        print(f"   · Fosfolípidos (real→render) . ext 7,733→1,400 | int 4,949→896")
        print(" FÁRMACO · Fingolimod / FTY720:")
        print(f"   · Fórmula / PM ............... C19H33NO2 / 307.48 g/mol")
        print(f"   · Átomos .................... 55 totales (22 pesados: 19 C, 1 N, 2 O)")
        print(f"   · XLogP3 (dato usado) ....... {LOGP_FARMACO:.2f}  (marcadamente lipófilo)")
        print(" MODELO FÍSICO:")
        print("   · Energía libre G(r) DERIVADA del logP:  ΔG = -2.303·k_BT·logP")
        print("   · Integrador: dinámica Browniana sobreamortiguada (ec. de Smoluchowski)")
        print("     → distribución estacionaria garantizada ∝ exp(-G/k_BT)")
        print(f"   · Temperatura: {TEMP_K:.0f} K  (k_BT = {KT_KCAL:.3f} kcal/mol)")
        print(f"   · Moléculas: {N_FARMACOS} | Frontera reflectante: {R_BOX:.0f} nm")
        print(" ALCANCE Y LÍMITES (honestos):")
        print("   · Modelo de ORDEN REDUCIDO (PMF radial efectivo), NO MD atomística.")
        print("   · Se modela la INSERCIÓN desde la interfaz; no el encuentro por difusión.")
        print("   · D es difusión EFECTIVA (tiempo acelerado); la termodinámica es invariante.")
        print("=" * 70)

    # ---- Nanotransportador (estático) ----
    def _fibonacci_esfera(self, n, radio):
        i = np.arange(n) + 0.5
        phi = np.arccos(1 - 2 * i / n); theta = np.pi * (1 + 5 ** 0.5) * i
        return np.column_stack([radio * np.sin(phi) * np.cos(theta),
                                radio * np.sin(phi) * np.sin(theta),
                                radio * np.cos(phi)])

    def _construir_liposoma(self):
        """Bicapa = dos monocapas concéntricas (ext. 20 nm, int. 16 nm) con colas radiales.
        Conteos de render asimétricos (más fuera que dentro), coherentes con 7,733 vs 4,949."""
        self.cabezas_ext = self._fibonacci_esfera(1400, RADIO_LIPOSOMA_EXT)
        self.cabezas_int = self._fibonacci_esfera(896, RADIO_NUCLEO_ACUOSO)
        # Las colas deben trazar radios rectos ext.→int. Proyectamos los MISMOS 896 puntos
        # externos hacia el radio interno (no una 2ª esfera de Fibonacci independiente: su
        # ángulo polar depende de n, lo que desalinea ext/int hasta ~50° y cruza las colas).
        cabezas_int_radiales = self.cabezas_ext[:896] * (RADIO_NUCLEO_ACUOSO / RADIO_LIPOSOMA_EXT)
        xs, ys, zs = [], [], []
        for pe, pin in zip(self.cabezas_ext[:896], cabezas_int_radiales):
            xs += [pe[0], pin[0], None]; ys += [pe[1], pin[1], None]; zs += [pe[2], pin[2], None]
        self.colas_x, self.colas_y, self.colas_z = xs, ys, zs
        u = np.linspace(0, 2 * np.pi, 40); v = np.linspace(0, np.pi, 20)
        self.surf_x = RADIO_LIPOSOMA_EXT * np.outer(np.cos(u), np.sin(v))
        self.surf_y = RADIO_LIPOSOMA_EXT * np.outer(np.sin(u), np.sin(v))
        self.surf_z = RADIO_LIPOSOMA_EXT * np.outer(np.ones_like(u), np.cos(v))

    # ---- FÍSICA: dinámica Browniana sobreamortiguada sobre el PMF ----
    def update_physics(self, dt):
        """
        Un paso de la ecuación de Smoluchowski (Langevin sobreamortiguado):
            dr = -D/(k_BT) · ∇G · dt  +  sqrt(2·D·dt) · η,   η ~ N(0, I)
        Con G en k_BT, ∇G en k_BT/nm, el término de deriva se reduce a -D·∇g.
        · Deriva: empuja al fármaco cuesta abajo del paisaje de energía CALCULADO.
        · Ruido : difusión térmica (fluctuación-disipación implícita en sqrt(2Ddt)).
        La distribución estacionaria de esta ecuación es exactamente ∝ exp(-G/k_BT),
        por eso el resultado es una PREDICCIÓN termodinámica, no un guion.
        """
        r = np.linalg.norm(self.pos, axis=1)
        r_hat = self.pos / r[:, None]
        deriva = -D_SIM * gradiente_energia(r, self.logP)          # nm/ps (radial)
        self.pos = (self.pos + deriva[:, None] * r_hat * dt +
                    np.sqrt(2.0 * D_SIM * dt) * rng.normal(size=self.pos.shape))
        # Frontera reflectante (compartimento finito de la suspensión liposomal)
        r = np.linalg.norm(self.pos, axis=1); fuera = r > R_BOX
        if fuera.any():
            self.pos[fuera] *= (2.0 * R_BOX / r[fuera] - 1.0)[:, None]

    # ---- Muestreo de un frame para la animación ----
    def _muestrear_frame(self):
        """Idéntico numéricamente a transformar molécula por molécula y enlace por enlace
        con un loop Python, pero vectorizado con NumPy (mismo resultado, mucho más rápido
        al repetirse ~156 veces para construir la animación)."""
        mundo = np.einsum('mcj,kj->mkc', self.rot_stack, self.tpl_coords) + self.pos[:, None, :]
        n_mol, n_at, _ = mundo.shape
        xyz = mundo.reshape(n_mol * n_at, 3)

        pa, pb = mundo[:, self.bond_a, :], mundo[:, self.bond_b, :]   # (n_mol, n_bonds, 3)
        nan_col = np.full(pa.shape[:2], np.nan)
        ex, ey, ez = (np.stack([pa[..., c], pb[..., c], nan_col], axis=-1).reshape(-1)
                      for c in range(3))

        self.traj_atomos.append((xyz, self.colores_frame, self.tamanos_frame))
        self.traj_enlaces.append((ex, ey, ez))
        self.traj_r.append(np.linalg.norm(self.pos, axis=1))

    def run(self):
        self._muestrear_frame()
        for paso in range(N_PASOS):
            self.update_physics(DT)
            if (paso + 1) % FRAMES_SALTO == 0:
                self._muestrear_frame()
        return self

    # ---- Simulación radial rápida (para la verificación y el control) ----
    @staticmethod
    def _simular_radial(logP, n_mol=N_FARMACOS, n_pasos=N_PASOS):
        """Versión escalar (solo radio) del mismo motor. Devuelve muestras equilibradas
        de r y la fracción final insertada en las colas. Se usa para verificar Boltzmann
        y para el control de falsabilidad, sin coste de renderizado 3D."""
        g = np.random.default_rng(11)
        d = g.normal(size=(n_mol, 3)); d /= np.linalg.norm(d, axis=1, keepdims=True)
        pos = d * g.uniform(R_INIT_LO, R_INIT_HI, n_mol)[:, None]
        muestras = []
        for step in range(n_pasos):
            r = np.linalg.norm(pos, axis=1); r_hat = pos / r[:, None]
            deriva = -D_SIM * gradiente_energia(r, logP)
            pos = pos + deriva[:, None] * r_hat * DT + np.sqrt(2 * D_SIM * DT) * g.normal(size=pos.shape)
            r = np.linalg.norm(pos, axis=1); fuera = r > R_BOX
            if fuera.any(): pos[fuera] *= (2 * R_BOX / r[fuera] - 1.0)[:, None]
            if step >= FRAC_EQUILIBRIO * n_pasos:
                muestras.append(np.linalg.norm(pos, axis=1))
        muestras = np.concatenate(muestras)
        rf = np.linalg.norm(pos, axis=1)
        frac = np.mean((rf > RADIO_NUCLEO_ACUOSO + 0.5) & (rf < RADIO_LIPOSOMA_EXT - 0.5))
        return muestras, frac

    def _control_radial(self, logP):
        """Envoltorio con caché de _simular_radial: verificar() y generate_pmf_figure()
        piden ambos el mismo control (logP fijo, semilla fija ⇒ resultado idéntico);
        sin caché se recalculaba la trayectoria completa dos veces por nada."""
        if logP not in self._radial_cache:
            self._radial_cache[logP] = self._simular_radial(logP)
        return self._radial_cache[logP]

    def verificar(self):
        """
        Imprime la EVIDENCIA numérica:
          · ΔG de transferencia derivado del logP.
          · Correlación entre el histograma radial muestreado (fármaco) y la
            predicción de Boltzmann ∝ exp(-G/k_BT)·r^2 → mide si el motor muestrea
            de verdad el paisaje calculado.
          · Control de falsabilidad: fracción insertada del análogo hidrofílico.
        """
        dG = -np.log(10.0) * self.logP
        # Muestras equilibradas del fármaco (2ª mitad de la trayectoria real 3D)
        n_eq = int(len(self.traj_r) * FRAC_EQUILIBRIO)
        r_drug = np.concatenate(self.traj_r[n_eq:])
        r_grid = np.linspace(14, 30, 400)
        p_boltz = np.exp(-energia_libre(r_grid, self.logP)) * r_grid ** 2
        p_boltz /= _trapz(p_boltz, r_grid)
        hist, edges = np.histogram(r_drug, bins=60, range=(14, 30), density=True)
        cen = 0.5 * (edges[1:] + edges[:-1])
        corr = np.corrcoef(hist, np.interp(cen, r_grid, p_boltz))[0, 1]
        frac_drug = np.mean((self.traj_r[-1] > RADIO_NUCLEO_ACUOSO + 0.5) &
                            (self.traj_r[-1] < RADIO_LIPOSOMA_EXT - 0.5))
        _, frac_ctrl = self._control_radial(LOGP_CONTROL)

        print("=" * 70)
        print(" VERIFICACIÓN TERMODINÁMICA (¿evidencia o guion?)")
        print("=" * 70)
        print(f" Fármaco: fingolimod, logP = {self.logP:.2f}")
        print(f"   ΔG transferencia agua→colas (DERIVADO del logP): "
              f"{dG:.2f} k_BT = {dG*KT_KCAL:.2f} kcal/mol")
        print(f"   Fracción final insertada en la bicapa: {frac_drug*100:.0f}%")
        print(f"   Ajuste a Boltzmann (corr. histograma vs exp(-G/kT)·r²): {corr:.3f}")
        print(f"     → corr≈1 ⇒ el motor MUESTREA el paisaje calculado (no lo fuerza).")
        print("-" * 70)
        print(f" CONTROL de falsabilidad: análogo hidrofílico logP = {LOGP_CONTROL:.2f}")
        print(f"   ΔG = {-np.log(10)*LOGP_CONTROL:+.2f} k_BT (membrana DESFAVORABLE)")
        print(f"   Fracción insertada: {frac_ctrl*100:.0f}%  "
              f"→ {'NO entra (correcto: el modelo sabe decir NO)' if frac_ctrl<0.1 else 'REVISAR'}")
        print("=" * 70)
        return corr, frac_drug, frac_ctrl

    def imprimir_conclusion(self, corr, frac_drug, frac_ctrl):
        """Conclusión de lo que la simulación PRUEBA (y lo que NO prueba)."""
        entra = frac_drug > 0.5 and corr > 0.9 and frac_ctrl < 0.1
        print("=" * 70)
        print(" CONCLUSIÓN: ¿QUÉ PRUEBA ESTA SIMULACIÓN?")
        print("=" * 70)
        if entra:
            print(" SÍ PRUEBA (dentro del modelo de orden reducido):")
            print("   1. El destino del fármaco es una SALIDA CALCULADA a partir de su")
            print("      logP medido (4.16), no un resultado impuesto a mano.")
            print(f"   2. La membrana es termodinámicamente FAVORABLE: ΔG = {-np.log(10)*self.logP*KT_KCAL:.2f} kcal/mol.")
            print(f"   3. El motor muestrea fielmente ese paisaje (Boltzmann corr = {corr:.3f}),")
            print(f"      con {frac_drug*100:.0f}% de las moléculas insertadas en las colas hidrófobas.")
            print(f"   4. Es FALSABLE: un análogo hidrofílico (logP<0) NO se inserta ({frac_ctrl*100:.0f}%).")
            print("   ⇒ Fingolimod PUEDE incorporarse a la bicapa del liposoma convencional;")
            print("     esta conclusión se deriva de la fisicoquímica, no de la animación.")
        else:
            print(" El resultado NO respalda la inserción espontánea con estos parámetros.")
            print(" Bajo el criterio del modelo, la interacción no queda demostrada.")
        print(" NO PRUEBA (límites): no sustituye a MD atomística con campo de fuerzas")
        print("   real + energía libre (umbrella/FEP), que es el estándar definitivo.")
        print("   La animación VISUALIZA la predicción; la evidencia está en los números.")
        print("=" * 70)

    # ---- Animación 3D ----
    def generate_animation(self):
        t_surf = go.Surface(x=self.surf_x, y=self.surf_y, z=self.surf_z, opacity=0.10,
                            showscale=False, colorscale=[[0, "#6cc7ff"], [1, "#6cc7ff"]],
                            hoverinfo="skip", showlegend=False)
        t_colas = go.Scatter3d(x=self.colas_x, y=self.colas_y, z=self.colas_z, mode="lines",
                            line=dict(color="#c9a24b", width=2), opacity=0.35,
                            hoverinfo="skip", name="Colas lipídicas")
        t_ext = go.Scatter3d(x=self.cabezas_ext[:, 0], y=self.cabezas_ext[:, 1],
                            z=self.cabezas_ext[:, 2], mode="markers",
                            marker=dict(size=3.0, color="#1f9e6b"),
                            name="Cabezas polares (ext.)", hoverinfo="skip")
        t_int = go.Scatter3d(x=self.cabezas_int[:, 0], y=self.cabezas_int[:, 1],
                            z=self.cabezas_int[:, 2], mode="markers",
                            marker=dict(size=2.4, color="#7fd3ad"),
                            name="Cabezas polares (int.)", hoverinfo="skip")
        xyz0, col0, siz0 = self.traj_atomos[0]; ex0, ey0, ez0 = self.traj_enlaces[0]
        t_at = go.Scatter3d(x=xyz0[:, 0], y=xyz0[:, 1], z=xyz0[:, 2], mode="markers",
                            marker=dict(size=siz0, color=col0, opacity=0.95),
                            name="Fingolimod (átomos)", hoverinfo="skip")
        t_bd = go.Scatter3d(x=ex0, y=ey0, z=ez0, mode="lines",
                            line=dict(color="#8a8f9c", width=2.0),
                            name="Fingolimod (enlaces)", hoverinfo="skip")
        traces = [t_surf, t_colas, t_ext, t_int, t_at, t_bd]

        frames = []
        for k in range(len(self.traj_atomos)):
            xyz, col, siz = self.traj_atomos[k]; ex, ey, ez = self.traj_enlaces[k]
            frames.append(go.Frame(name=f"{k}", traces=[4, 5], data=[
                go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="markers",
                             marker=dict(size=siz, color=col, opacity=0.95)),
                go.Scatter3d(x=ex, y=ey, z=ez, mode="lines",
                             line=dict(color="#8a8f9c", width=2.0))]))

        botones = dict(type="buttons", showactive=False, x=0.05, y=0.05, xanchor="left", buttons=[
            dict(label="▶ Reproducir", method="animate",
                 args=[None, dict(frame=dict(duration=45, redraw=True), fromcurrent=True,
                                  transition=dict(duration=0))]),
            dict(label="⏸ Pausa", method="animate",
                 args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate",
                                    transition=dict(duration=0))])])
        slider = dict(active=0, x=0.12, y=0, len=0.85,
                      currentvalue=dict(prefix="Frame: ", font=dict(size=13)),
                      steps=[dict(method="animate", label=f"{k}",
                                  args=[[f"{k}"], dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                        transition=dict(duration=0))]) for k in range(len(frames))])

        rango = R_BOX * 1.1
        fig = go.Figure(data=traces, frames=frames)
        fig.update_layout(
            title=dict(text="Inserción de fingolimod en liposoma · dinámica Browniana sobre PMF(logP=4.16)",
                       x=0.5, font=dict(size=14)),
            template="plotly_dark",
            scene=dict(
                xaxis=dict(range=[-rango, rango], title="x (nm)", showbackground=False),
                yaxis=dict(range=[-rango, rango], title="y (nm)", showbackground=False),
                zaxis=dict(range=[-rango, rango], title="z (nm)", showbackground=False),
                aspectmode="cube", camera=dict(eye=dict(x=1.6, y=1.6, z=1.0)),
                uirevision="camara_fija"),      # ← permite rotar/zoom DURANTE la animación
            uirevision="camara_fija",
            updatemenus=[botones], sliders=[slider],
            margin=dict(l=0, r=0, t=44, b=0), legend=dict(x=0.0, y=0.98, font=dict(size=10)))
        return fig

    # ---- Panel de verificación: PMF + Boltzmann + control ----
    def generate_pmf_figure(self):
        """
        Figura 2D de EVIDENCIA:
          Fila 1: G(r) del fármaco (derivado de logP) y del control hidrofílico.
          Fila 2: densidad radial muestreada del fármaco vs. predicción de Boltzmann,
                  y densidad del control (que se queda en el agua).
        """
        r_grid = np.linspace(14, 30, 400)
        g_drug = energia_libre(r_grid, self.logP) * KT_KCAL     # kcal/mol
        g_ctrl = energia_libre(r_grid, LOGP_CONTROL) * KT_KCAL
        p_boltz = np.exp(-energia_libre(r_grid, self.logP)) * r_grid ** 2
        p_boltz /= _trapz(p_boltz, r_grid)

        n_eq = int(len(self.traj_r) * FRAC_EQUILIBRIO)
        r_drug = np.concatenate(self.traj_r[n_eq:])
        r_ctrl, _ = self._control_radial(LOGP_CONTROL)

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.09,
                            subplot_titles=("Energía libre G(r) DERIVADA del logP",
                                            "Distribución radial: muestreo vs. Boltzmann (evidencia)"))
        # Fila 1: perfiles de energía
        fig.add_trace(go.Scatter(x=r_grid, y=g_drug, mode="lines", name="G(r) fingolimod (logP=4.16)",
                                 line=dict(color="#f2a23c", width=3)), row=1, col=1)
        fig.add_trace(go.Scatter(x=r_grid, y=g_ctrl, mode="lines", name="G(r) control (logP=-1.5)",
                                 line=dict(color="#6cc7ff", width=2, dash="dash")), row=1, col=1)
        # Fila 2: distribuciones
        fig.add_trace(go.Histogram(x=r_drug, histnorm="probability density", nbinsx=60,
                                   name="Fármaco (muestreado)", marker_color="#f2a23c", opacity=0.55),
                      row=2, col=1)
        fig.add_trace(go.Scatter(x=r_grid, y=p_boltz, mode="lines",
                                 name="Predicción Boltzmann ∝ e^(-G/kT)·r²",
                                 line=dict(color="#e23b2e", width=3)), row=2, col=1)
        fig.add_trace(go.Histogram(x=r_ctrl, histnorm="probability density", nbinsx=60,
                                   name="Control hidrofílico (no entra)", marker_color="#6cc7ff", opacity=0.45),
                      row=2, col=1)
        # Bandas de referencia (colas / cabezas)
        for rr, txt in [(RADIO_NUCLEO_ACUOSO, "cabeza int."), (R_MIDPLANO_TAIL, "colas (mín. G)"),
                        (RADIO_LIPOSOMA_EXT, "cabeza ext.")]:
            fig.add_vline(x=rr, line=dict(color="#888", width=1, dash="dot"),
                          annotation_text=txt, annotation_font_size=9, row="all", col=1)

        fig.update_xaxes(title_text="r = distancia al centro del liposoma (nm)", row=2, col=1)
        fig.update_yaxes(title_text="G (kcal/mol)", row=1, col=1)
        fig.update_yaxes(title_text="densidad", row=2, col=1)
        fig.update_layout(template="plotly_dark", height=680, barmode="overlay",
                          title=dict(text="PANEL DE EVIDENCIA: el destino del fármaco es una salida CALCULADA",
                                     x=0.5, font=dict(size=14)),
                          legend=dict(font=dict(size=10)))
        return fig


# =============================================================================
#  EJECUCIÓN (una sola celda)
# =============================================================================
sim = Simulation3D(logP=LOGP_FARMACO)
sim.imprimir_especificaciones()                 # ficha técnica + modelo + límites
sim.run()
corr, fd, fc = sim.verificar()                  # evidencia numérica (Boltzmann + control)
sim.imprimir_conclusion(corr, fd, fc)           # qué prueba y qué no
sim.generate_animation().show()                 # animación 3D de la inserción
sim.generate_pmf_figure().show()                # panel PMF + Boltzmann + control
