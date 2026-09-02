# -*- coding: utf-8 -*-
"""
nk_lib.py — Núcleo compartido de los ejercicios de la exposición.

Tema: Estrategias para la convergencia global — Métodos de Newton-Krylov
Materia: Métodos Numéricos II (DAT-252) — UMSA, Carrera de Informática
Expositores: Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

Este módulo implementa, de forma didáctica pero honesta, las piezas de un
solver de Newton-Krylov globalizado:

    1. Producto Jacobiano-vector SIN construir el Jacobiano (diferencias finitas).
    2. Términos de forzado de Eisenstat-Walker (control del "cuánto resolver").
    3. Búsqueda de línea con retroceso y condición de Armijo.
    4. Región de confianza con CG truncado de Steihaug.
    5. Continuación pseudo-transitoria (Psi-tc).
    6. El bucle de Newton-Krylov que las combina todas.

Referencias principales:
    - C. T. Kelley, "Iterative Methods for Linear and Nonlinear Equations", SIAM 1995.
    - C. T. Kelley, "Solving Nonlinear Equations with Newton's Method", SIAM 2003.
    - Dembo, Eisenstat & Steihaug, "Inexact Newton Methods", SINUM 19(2), 1982.
    - Eisenstat & Walker, "Choosing the forcing terms in an inexact Newton
      method", SISC 17(1), 1996.
    - Knoll & Keyes, "Jacobian-free Newton-Krylov methods: a survey of
      approaches and applications", JCP 193, 2004.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

EPS_MAQ = np.finfo(float).eps          # ~2.22e-16
SQRT_EPS = math.sqrt(EPS_MAQ)          # ~1.49e-8


# =============================================================================
# 1. Contador de trabajo
# =============================================================================
class ContadorF:
    """Envuelve F(x) y cuenta cuántas veces se la evalúa.

    En un método de Newton-Krylov matriz-libre, CADA producto J·v cuesta una
    evaluación de F. Por eso la métrica honesta de costo no es "iteraciones de
    Newton" sino "evaluaciones de F". Es el mensaje central del ejercicio 2.
    """

    def __init__(self, F):
        self.F = F
        self.n_F = 0        # evaluaciones totales de F
        self.n_Jv = 0       # de esas, cuántas vinieron de un producto J·v

    def __call__(self, x):
        self.n_F += 1
        return self.F(np.asarray(x, dtype=float))

    def reset(self):
        self.n_F = 0
        self.n_Jv = 0


# =============================================================================
# 2. Producto Jacobiano-vector matriz-libre
# =============================================================================
def jv_diferencias_finitas(F, x, Fx, v):
    """Aproxima J(x)·v por diferencias finitas hacia adelante, sin formar J.

        J(x)·v  ≈  [ F(x + ε v) − F(x) ] / ε

    La elección de ε balancea dos errores opuestos:
      · error de truncamiento  ~ ε‖v‖‖F''‖/2      (crece con ε)
      · error de redondeo      ~ eps_maq‖F‖/(ε‖v‖) (crece al bajar ε)
    El óptimo está cerca de ε ≈ √eps_maq ≈ 1.49e-8, escalado por la magnitud
    de x. Usamos la fórmula de Kelley (2003, §3.2.1), que además preserva el
    signo para no cancelar cifras significativas.

    Devuelve el vector J·v. Cuesta UNA evaluación de F.
    """
    v = np.asarray(v, dtype=float)
    norma_v = np.linalg.norm(v)
    if norma_v == 0.0:
        return np.zeros_like(x)

    # Escala: proyección de x sobre la dirección v (Kelley).
    xs = float(np.dot(x, v)) / norma_v
    eps = SQRT_EPS
    if xs != 0.0:
        eps = eps * max(abs(xs), 1.0) * math.copysign(1.0, xs)
    eps = eps / norma_v

    return (F(x + eps * v) - Fx) / eps


def operador_jacobiano(F, x, Fx, contador=None):
    """Construye el LinearOperator de SciPy que representa J(x) sin almacenarla.

    Esto es lo que permite usar GMRES sobre una matriz que nunca existe en
    memoria: GMRES solo pide productos J·v, y cada uno se fabrica al vuelo.
    """
    n = x.size

    def matvec(v):
        if contador is not None:
            contador.n_Jv += 1
        return jv_diferencias_finitas(F, x, Fx, v)

    return LinearOperator((n, n), matvec=matvec, dtype=float)


# =============================================================================
# 3. Términos de forzado (Eisenstat-Walker)
# =============================================================================
GAMMA_EW = 0.9
ALPHA_EW = (1.0 + math.sqrt(5.0)) / 2.0     # razón áurea ≈ 1.618


def forzado_eisenstat_walker(tipo, normaF, normaF_prev, eta_prev,
                             residuo_lineal_prev=None,
                             eta_max=0.9, tol_abs=None):
    """Calcula η_k: cuán inexactamente se puede resolver el sistema de Newton.

    El criterio de Newton inexacto es   ‖J s + F‖ ≤ η ‖F‖ .
    Teorema (Dembo-Eisenstat-Steihaug 1982):
        η_k ≤ η < 1 constante  ⇒ convergencia LINEAL local
        η_k → 0                ⇒ convergencia SUPERLINEAL
        η_k = O(‖F_k‖)         ⇒ convergencia CUADRÁTICA (recupera Newton)

    Eisenstat-Walker proponen elegir η_k adaptativamente para no "resolver de
    más" (oversolving) lejos de la raíz, donde el modelo lineal no vale nada:

      Choice 1:  η_k = | ‖F_k‖ − ‖F_{k−1} + J_{k−1}s_{k−1}‖ | / ‖F_{k−1}‖
      Choice 2:  η_k = γ (‖F_k‖ / ‖F_{k−1}‖)^α ,  γ=0.9, α=1.618

    La salvaguarda evita que η caiga demasiado rápido tras un buen paso suelto.
    """
    if isinstance(tipo, (int, float)):
        return float(tipo)

    if normaF_prev is None or normaF_prev == 0.0:
        return eta_max                        # primer paso: resolver flojo

    if tipo == "ew1":
        if residuo_lineal_prev is None:
            eta = GAMMA_EW * (normaF / normaF_prev) ** ALPHA_EW
        else:
            eta = abs(normaF - residuo_lineal_prev) / normaF_prev
    elif tipo == "ew2":
        eta = GAMMA_EW * (normaF / normaF_prev) ** ALPHA_EW
    else:
        raise ValueError(f"Término de forzado desconocido: {tipo!r}")

    # Salvaguarda de Eisenstat-Walker: no bajar η más rápido de lo sensato.
    if eta_prev is not None:
        piso = GAMMA_EW * eta_prev ** ALPHA_EW
        if piso > 0.1:
            eta = max(eta, piso)

    eta = min(eta_max, max(eta, 0.0))

    # Última salvaguarda: no resolver mucho más fino que la tolerancia final.
    if tol_abs is not None and normaF > 0:
        eta = min(eta_max, max(eta, 0.5 * tol_abs / normaF))

    return eta


# =============================================================================
# 4. Globalización A — búsqueda de línea con retroceso (Armijo)
# =============================================================================
def linea_armijo(F, x, s, normaF, alpha=1e-4, lam_min=1e-10, max_retro=40):
    """Busca λ ∈ (0,1] tal que el paso x + λs reduzca lo suficiente ‖F‖.

    Función de mérito:  f(x) = ½‖F(x)‖².
    Condición de Armijo (forma en norma del residuo, Kelley 2003 §1.6):

        ‖F(x + λ s)‖ ≤ (1 − α λ) ‖F(x)‖ ,   α = 1e-4

    Es legítima porque el paso de Newton inexacto es dirección de descenso:
    si J s = −F + r con ‖r‖ ≤ η‖F‖ y η < 1, entonces
        ∇f(x)ᵀ s = FᵀJ s = −‖F‖² + Fᵀr ≤ −(1 − η)‖F‖² < 0 .

    El retroceso usa interpolación cuadrática con salvaguardas σ0=0.1, σ1=0.5:
    en vez de partir λ a la mitad a ciegas, ajusta una parábola a la
    información que ya se tiene y salta al mínimo, acotado al intervalo seguro.

    Devuelve (x_nuevo, F_nuevo, normaF_nueva, lam, n_retrocesos, exito).
    """
    sigma0, sigma1 = 0.1, 0.5
    lam = 1.0
    x_t = x + lam * s
    F_t = F(x_t)
    n_t = np.linalg.norm(F_t)

    n_retro = 0
    lam_prev, n_prev = lam, n_t

    # Ojo con la comparación: si el paso completo desborda, n_t es inf o nan.
    # Escrita como "n_t > cota" un nan daría False y ACEPTARÍAMOS el paso
    # (nan no es mayor que nada). Escrita como "not (n_t <= cota)" un nan
    # entra al retroceso, que es lo correcto.
    while not (n_t <= (1.0 - alpha * lam) * normaF):
        if lam < lam_min or n_retro >= max_retro:
            return x_t, F_t, n_t, lam, n_retro, False

        # Modelo cuadrático de g(λ) = ‖F(x+λs)‖² con g(0), g'(0) y g(λ).
        if n_retro == 0 or not np.isfinite(n_t):
            lam_nuevo = 0.5 * lam
        else:
            g0 = normaF ** 2
            gl = n_t ** 2
            # derivada de g en 0 acotada por la condición de descenso
            dg0 = -2.0 * g0
            denom = 2.0 * (gl - g0 - dg0 * lam)
            lam_nuevo = -dg0 * lam * lam / denom if denom != 0 else 0.5 * lam

        # Salvaguardas: el nuevo λ debe estar en [σ0·λ, σ1·λ].
        lam_nuevo = min(max(lam_nuevo, sigma0 * lam), sigma1 * lam)

        lam_prev, n_prev = lam, n_t
        lam = lam_nuevo
        x_t = x + lam * s
        F_t = F(x_t)
        n_t = np.linalg.norm(F_t)
        n_retro += 1

    return x_t, F_t, n_t, lam, n_retro, True


# =============================================================================
# 5. Globalización B — región de confianza con CG truncado (Steihaug)
# =============================================================================
def steihaug_cg(Jop, Fx, delta, tol_rel=1e-4, max_iter=None):
    """CG truncado de Steihaug sobre el modelo de mínimos cuadrados de Gauss-Newton.

    Modelo:  m(s) = ½‖F + J s‖²  sujeto a ‖s‖ ≤ Δ.
    Aplicamos CG a las ecuaciones normales  JᵀJ s = −JᵀF, deteniéndonos si:
      (a) se alcanza la tolerancia,
      (b) se detecta curvatura ≤ 0  → se sigue la dirección hasta el borde,
      (c) el iterado se sale de la región → se corta en el borde.

    La gracia: cada iteración necesita solo J·v y Jᵀ·w. Como aquí J es
    matriz-libre, esto encaja exactamente con la misma maquinaria de Krylov.
    Nota: usamos Jᵀ vía el operador adjunto de SciPy cuando existe; si no,
    se cae a la variante con J·v solamente (Gauss-Newton aproximado).
    """
    n = Fx.size
    if max_iter is None:
        max_iter = min(n, 60)

    s = np.zeros(n)
    # Gradiente del modelo en s=0:  g = Jᵀ F
    g = Jop.rmatvec(Fx) if hasattr(Jop, "rmatvec") else Jop.matvec(Fx)
    norma_g0 = np.linalg.norm(g)
    if norma_g0 == 0.0:
        return s, "gradiente nulo"

    r = -g
    p = r.copy()
    rr = float(np.dot(r, r))
    umbral = tol_rel * norma_g0

    for _ in range(max_iter):
        Jp = Jop.matvec(p)
        curv = float(np.dot(Jp, Jp))          # pᵀ(JᵀJ)p ≥ 0 siempre
        if curv <= 1e-300:
            # Curvatura nula: ir al borde en la dirección p.
            tau = _tau_borde(s, p, delta)
            return s + tau * p, "curvatura nula"

        alpha = rr / curv
        s_nuevo = s + alpha * p
        if np.linalg.norm(s_nuevo) >= delta:
            tau = _tau_borde(s, p, delta)
            return s + tau * p, "borde de la región"

        s = s_nuevo
        JtJp = Jop.rmatvec(Jp) if hasattr(Jop, "rmatvec") else Jp
        r = r - alpha * JtJp
        rr_nuevo = float(np.dot(r, r))
        if math.sqrt(rr_nuevo) < umbral:
            return s, "tolerancia alcanzada"
        p = r + (rr_nuevo / rr) * p
        rr = rr_nuevo

    return s, "máx. iteraciones"


def _tau_borde(s, p, delta):
    """τ ≥ 0 tal que ‖s + τp‖ = Δ (raíz positiva de la cuadrática)."""
    a = float(np.dot(p, p))
    b = 2.0 * float(np.dot(s, p))
    c = float(np.dot(s, s)) - delta ** 2
    disc = max(b * b - 4 * a * c, 0.0)
    return (-b + math.sqrt(disc)) / (2 * a) if a > 0 else 0.0


def paso_dogleg(Jadj, Fx, s_newton, delta):
    """Paso de dogleg (Powell) dentro de la región de confianza ‖s‖ ≤ Δ.

    Modelo de Gauss-Newton:  m(s) = ½‖F + J s‖² .
    Se interpola entre dos direcciones ya conocidas:

      · el punto de Cauchy  s_C = −(‖g‖²/‖J g‖²) g ,  con g = JᵀF = ∇f(x)
        (el mínimo del modelo a lo largo del máximo descenso), y
      · el paso de Newton inexacto s_N que ya devolvió GMRES.

    Reglas:
      ‖s_N‖ ≤ Δ    → se toma s_N (el paso de Newton cabe en la región)
      ‖s_C‖ ≥ Δ    → se recorta el máximo descenso al borde
      en otro caso → s = s_C + τ(s_N − s_C) con ‖s‖ = Δ,  τ ∈ [0,1]

    Por qué dogleg y no Steihaug-CG en el contexto matriz-libre: el CG de
    Steihaug resuelve las ecuaciones normales JᵀJ s = −JᵀF, lo que ELEVA AL
    CUADRADO el número de condición; con una J conocida solo por diferencias
    finitas (ruido relativo ~1e-8) el CG se estanca alrededor de 1e-6 y el
    método deja de progresar. Lo comprobamos experimentalmente. El dogleg
    necesita apenas dos productos extra y es lo que se usa en la práctica
    (Pawlowski, Shadid, Simonis & Walker, SIAM Review 48(4), 2006).

    Devuelve (s, etiqueta_del_tramo).
    """
    norma_sN = np.linalg.norm(s_newton)
    if norma_sN <= delta:
        return s_newton, "Newton"

    g = Jadj.rmatvec(Fx)                     # g = JᵀF = ∇f(x)
    Jg = Jadj.matvec(g)
    den = float(np.dot(Jg, Jg))
    if den <= 1e-300:
        return delta * s_newton / norma_sN, "Newton recortado"

    s_c = -(float(np.dot(g, g)) / den) * g   # punto de Cauchy
    norma_sC = np.linalg.norm(s_c)

    if norma_sC >= delta:
        return delta * s_c / norma_sC, "Cauchy (borde)"

    d = s_newton - s_c
    a = float(np.dot(d, d))
    b = 2.0 * float(np.dot(s_c, d))
    c = float(np.dot(s_c, s_c)) - delta ** 2
    disc = max(b * b - 4 * a * c, 0.0)
    tau = (-b + math.sqrt(disc)) / (2 * a) if a > 0 else 0.0
    tau = min(max(tau, 0.0), 1.0)
    return s_c + tau * d, "dogleg"


# =============================================================================
# 6. Historial de la iteración
# =============================================================================
@dataclass
class Historial:
    """Todo lo que hace falta para graficar y explicar lo que pasó."""
    metodo: str = ""
    x: np.ndarray | None = None
    convergio: bool = False
    motivo: str = ""
    residuales: list = field(default_factory=list)   # ‖F‖ por iteración
    etas: list = field(default_factory=list)         # término de forzado
    lambdas: list = field(default_factory=list)      # paso de la línea
    krylov_iters: list = field(default_factory=list) # iters de GMRES por paso
    trayectoria: list = field(default_factory=list)  # x_k (solo problemas 2D)
    n_F: int = 0
    n_Jv: int = 0

    @property
    def n_newton(self):
        return max(len(self.residuales) - 1, 0)

    @property
    def trabajo_acumulado(self):
        """Productos J·v acumulados tras cada iteración de Newton.

        Es el eje horizontal honesto para comparar estrategias: mide el
        trabajo real, no el número de iteraciones externas.
        """
        return np.cumsum([0] + self.krylov_iters)[:len(self.residuales)]


# =============================================================================
# 7. El solver: Newton-Krylov globalizado
# =============================================================================
def newton_krylov(F, x0, globalizacion="linea", forzado="ew2",
                  tol=1e-10, max_iter=60, precond=None,
                  gmres_restart=30, gmres_maxiter=200,
                  delta0=1.0, delta_max=1e3, ptc_delta0="auto",
                  jacobiano_simetrico=False, max_evals_F=None,
                  guardar_trayectoria=False, verbose=False):
    """Resuelve F(x) = 0 con Newton inexacto + GMRES matriz-libre + globalización.

    Parámetros
    ----------
    globalizacion : "ninguna" | "linea" | "region" | "ptc"
        "ninguna" → Newton inexacto puro (paso completo). Sirve para MOSTRAR
                    que sin globalización el método falla lejos de la raíz.
        "linea"   → búsqueda de línea con retroceso de Armijo.
        "region"  → región de confianza con CG truncado de Steihaug.
        "ptc"     → continuación pseudo-transitoria: resuelve
                    (I/δ + J)s = −F con δ creciendo según la regla SER
                    δ_k = δ_{k−1}·‖F_{k−1}‖/‖F_k‖. Equivale a integrar
                    dx/dt = −F(x) con Euler implícito y paso creciente.
    forzado : float | "ew1" | "ew2"
        Término de forzado η. Un float lo fija; "ew1"/"ew2" lo adaptan.
    precond : LinearOperator | None
        Precondicionador por la derecha para GMRES (aquí se pasa como M).

    Devuelve un Historial.
    """
    cont = ContadorF(F)
    x = np.asarray(x0, dtype=float).copy()
    Fx = cont(x)
    normaF = np.linalg.norm(Fx)

    h = Historial(metodo=f"{globalizacion}/η={forzado}")
    h.residuales.append(normaF)
    if guardar_trayectoria:
        h.trayectoria.append(x.copy())

    normaF0 = normaF
    eta_prev = None
    normaF_prev = None
    residuo_lineal_prev = None
    delta = delta0
    # δ0 de Ψtc: si no se fija, se escala con el residuo inicial. Ψtc es
    # sensible a este valor (se estudia en el ejercicio 3): demasiado chico
    # sobre-amortigua y el método se arrastra; demasiado grande equivale a
    # Newton puro y se pierde la robustez.
    ptc_delta = (1.0 / max(normaF, 1e-300)) if ptc_delta0 == "auto" else float(ptc_delta0)

    for k in range(max_iter):
        if normaF <= tol or normaF <= tol * max(normaF0, 1.0):
            h.motivo = "convergió"
            h.convergio = True
            break
        if not np.isfinite(normaF) or normaF > 1e14:
            h.motivo = "divergió (‖F‖ desbordó)"
            break
        if max_evals_F is not None and cont.n_F > max_evals_F:
            # Presupuesto de trabajo. Sin esto, una corrida que no converge
            # puede gastar minutos moliendo GMRES sin ir a ninguna parte.
            h.motivo = "presupuesto de evaluaciones agotado"
            break

        eta = forzado_eisenstat_walker(forzado, normaF, normaF_prev, eta_prev,
                                       residuo_lineal_prev, tol_abs=tol)
        h.etas.append(eta)

        Jop = operador_jacobiano(cont, x, Fx, contador=cont)

        # ---------------- cálculo del paso ----------------
        if globalizacion == "region":
            # 1) El paso de Newton inexacto se calcula igual que siempre.
            n_antes = cont.n_Jv
            n_kry_box = [0]

            def cb_tr(_res, box=n_kry_box):
                box[0] += 1

            s_newton, info = gmres(Jop, -Fx, rtol=eta, atol=0.0,
                                   restart=gmres_restart, maxiter=gmres_maxiter,
                                   M=precond, callback=cb_tr,
                                   callback_type="pr_norm")
            if info < 0 or not np.all(np.isfinite(s_newton)):
                h.motivo = "GMRES falló"
                break
            residuo_lineal_prev = np.linalg.norm(Fx + Jop.matvec(s_newton))
            # 2) Y el dogleg decide cuánto de ese paso cabe en la región.
            Jadj = operador_con_adjunto(cont, x, Fx, contador=cont,
                                        simetrico=jacobiano_simetrico)
            s, _tramo = paso_dogleg(Jadj, Fx, s_newton, delta)
            n_kry = max(cont.n_Jv - n_antes, 1)
        else:
            if globalizacion == "ptc":
                # Continuación pseudo-transitoria. Se integra el flujo
                #     dx/dt = −F(x)
                # con Euler implícito y paso δ. La ecuación por paso es
                #     (I + δ J) s = −δ F
                # (forma multiplicada por δ: mejor condicionada que I/δ + J
                #  en los dos extremos). δ→0 da el paso explícito s = −δF, o
                # sea máximo descenso amortiguado; δ→∞ recupera Newton exacto.
                base = Jop
                def matvec_ptc(v, base=base, d=ptc_delta):
                    return v + d * base.matvec(v)
                Aop = LinearOperator(Jop.shape, matvec=matvec_ptc, dtype=float)
            else:
                Aop = Jop

            n_kry_box = [0]

            def cb(_res, box=n_kry_box):
                box[0] += 1

            lado_derecho = -ptc_delta * Fx if globalizacion == "ptc" else -Fx
            s, info = gmres(Aop, lado_derecho, rtol=eta, atol=0.0,
                            restart=gmres_restart, maxiter=gmres_maxiter,
                            M=precond, callback=cb, callback_type="pr_norm")
            n_kry = n_kry_box[0]
            if info < 0 or not np.all(np.isfinite(s)):
                h.motivo = "GMRES falló"
                break
            # Residuo lineal real, para Choice 1 de Eisenstat-Walker.
            if globalizacion != "ptc":
                residuo_lineal_prev = np.linalg.norm(Fx + Jop.matvec(s))

        h.krylov_iters.append(max(n_kry, 1))

        # ---------------- globalización ----------------
        if globalizacion == "linea":
            x, Fx, normaF_nueva, lam, n_retro, ok = linea_armijo(cont, x, s, normaF)
            h.lambdas.append(lam)
            if not ok:
                h.motivo = "la búsqueda de línea falló (λ < λ_min)"
                h.residuales.append(normaF_nueva)
                break
            normaF_prev, normaF = normaF, normaF_nueva

        elif globalizacion == "region":
            # Aceptación por razón de reducción real vs. predicha.
            pred = 0.5 * normaF ** 2 - 0.5 * np.linalg.norm(Fx + Jop.matvec(s)) ** 2
            x_t = x + s
            F_t = cont(x_t)
            n_t = np.linalg.norm(F_t)
            real = 0.5 * normaF ** 2 - 0.5 * n_t ** 2
            rho = real / pred if pred > 0 else -1.0

            if rho < 0.25:
                delta = 0.25 * np.linalg.norm(s)
            elif rho > 0.75 and abs(np.linalg.norm(s) - delta) < 1e-8 * delta:
                delta = min(2.0 * delta, delta_max)

            if rho > 1e-4:
                x, Fx = x_t, F_t
                normaF_prev, normaF = normaF, n_t
                h.lambdas.append(1.0)
            else:
                h.lambdas.append(0.0)       # paso rechazado
                if delta < 1e-12:
                    h.motivo = "región de confianza colapsó"
                    break
                continue

        else:   # "ninguna" y "ptc": paso completo
            x = x + s
            Fx = cont(x)
            normaF_prev, normaF = normaF, np.linalg.norm(Fx)
            h.lambdas.append(1.0)
            if globalizacion == "ptc" and normaF > 0:
                # Regla SER: el paso pseudo-temporal crece al bajar el residuo.
                ptc_delta = ptc_delta * normaF_prev / normaF
                ptc_delta = min(ptc_delta, 1e12)

        eta_prev = eta
        h.residuales.append(normaF)
        if guardar_trayectoria:
            h.trayectoria.append(x.copy())
        if verbose:
            print(f"    k={k:2d}  ‖F‖={normaF:.3e}  η={eta:.2e}  "
                  f"λ={h.lambdas[-1]:.3f}  GMRES={h.krylov_iters[-1]}")
    else:
        h.motivo = "máx. iteraciones alcanzado"

    if normaF <= tol or normaF <= tol * max(normaF0, 1.0):
        h.convergio = True
        h.motivo = h.motivo or "convergió"
    h.x = x
    h.n_F = cont.n_F
    h.n_Jv = cont.n_Jv
    return h


def operador_con_adjunto(F, x, Fx, contador=None, simetrico=False, n_max_denso=60):
    """Operador J que además sabe aplicar Jᵀ, necesario para la región de confianza.

    El CG de Steihaug trabaja sobre las ecuaciones normales JᵀJ s = −JᵀF, así que
    necesita Jᵀ·w. Y ahí aparece una limitación real del enfoque matriz-libre:
    la diferencia finita hacia adelante da J·v con UNA evaluación de F, pero NO
    da Jᵀ·w. Hay dos salidas:

      (a) Si J es simétrica, Jᵀ·w = J·w y no hay nada que hacer. Este es
          exactamente nuestro caso en Bratu: el laplaciano con diferencias
          centradas es simétrico y el término de reacción es diagonal, así que
          J es simétrica de forma exacta (se verifica en el ejercicio 3).
      (b) Si no lo es y la dimensión es chica, se forma J columna a columna con
          n productos J·v y se usa J.T. Cuesta n evaluaciones de F: aceptable
          para n ≤ 60, impensable para n = 10 000.

    Esta función implementa ambas y —importante para las gráficas de costo—
    contabiliza todas las evaluaciones que consume.
    """
    n = x.size

    def matvec(v):
        if contador is not None:
            contador.n_Jv += 1
        return jv_diferencias_finitas(F, x, Fx, v)

    if simetrico:
        return LinearOperator((n, n), matvec=matvec, rmatvec=matvec, dtype=float)

    if n <= n_max_denso:
        J = np.empty((n, n))
        base = np.eye(n)
        for j in range(n):
            if contador is not None:
                contador.n_Jv += 1
            J[:, j] = jv_diferencias_finitas(F, x, Fx, base[:, j])
        return LinearOperator((n, n), matvec=lambda v: J @ v,
                              rmatvec=lambda w: J.T @ w, dtype=float)

    # Última opción: suponer simetría (válido en difusión-reacción centrada).
    return LinearOperator((n, n), matvec=matvec, rmatvec=matvec, dtype=float)


# =============================================================================
# 8. Utilidades de presentación
# =============================================================================
def tabla(filas, encabezados, anchos=None):
    """Imprime una tabla ASCII legible desde el proyector."""
    if anchos is None:
        anchos = [max(len(str(encabezados[i])),
                      max((len(str(f[i])) for f in filas), default=0)) + 2
                  for i in range(len(encabezados))]
    linea = "+".join("-" * a for a in anchos)
    print("  " + linea)
    print("  " + "|".join(str(encabezados[i]).center(anchos[i])
                          for i in range(len(encabezados))))
    print("  " + linea)
    for f in filas:
        print("  " + "|".join(str(f[i]).center(anchos[i])
                              for i in range(len(encabezados))))
    print("  " + linea)


def titulo(texto):
    print()
    print("=" * 78)
    print(f"  {texto}")
    print("=" * 78)


def estilo_figuras():
    """Estilo común: legible desde el fondo del aula."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 130,
        "font.size": 11,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "legend.framealpha": 0.9,
        "figure.autolayout": True,
    })
    return plt
