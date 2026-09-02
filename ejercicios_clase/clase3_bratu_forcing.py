# -*- coding: utf-8 -*-
"""
=============================================================================
 EJERCICIO 3 PARA LA CLASE — Newton-Krylov completo sobre la ecuación de Bratu
=============================================================================
 Materia : Métodos Numéricos II (DAT-252) — UMSA
 Tema    : Estrategias para la convergencia global · Métodos de Newton-Krylov

 EL PROBLEMA
   Ecuación de Bratu (modelo de ignición térmica) en 1D:

       u''(x) + λ e^{u(x)} = 0 ,      u(0) = u(1) = 0

   Discretizada con diferencias finitas centradas queda un sistema no lineal
   de N ecuaciones. Con N = 200 el Jacobiano tendría 40 000 entradas.

 QUÉ SE ESTUDIA AQUÍ
   El método completo, con sus tres decisiones de diseño puestas una al lado
   de la otra:

     1. Cómo se elige el término de forzado η (fijo o de Eisenstat-Walker).
     2. Si se usa precondicionador o no.
     3. Qué pasa cuando λ se acerca al punto de retorno λ* ≈ 3.5138, donde el
        problema deja de tener solución.

   Todo se mide en PRODUCTOS J·v, que es el trabajo real, porque cada
   producto J·v cuesta exactamente una evaluación de F.

 REQUISITOS
   pip install numpy scipy matplotlib

 EJECUTAR
   python3 clase3_bratu_forcing.py        (~30 s)
=============================================================================
"""

import math
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import LinearOperator, gmres

warnings.filterwarnings("ignore")

N = 200
H = 1.0 / (N + 1)
X = np.linspace(0.0, 1.0, N + 2)[1:-1]
SQRT_EPS = math.sqrt(np.finfo(float).eps)

# El punto de retorno del problema continuo. Para λ > λ* no existe solución.
LAMBDA_CRITICO = 3.513830719


def hacer_F(lam):
    """Devuelve el residuo discreto de Bratu para un λ dado.

    F_i(u) = −[ (u_{i−1} − 2u_i + u_{i+1})/h² + λ e^{u_i} ] · h²

    El signo menos y el factor h² están puestos a propósito: con ellos el
    Jacobiano queda definido positivo y con entradas de orden 1. Newton no
    nota ninguna de las dos cosas, pero la función de mérito ½‖F‖² sí.
    """
    def F(u):
        ub = np.zeros(N + 2)
        ub[1:-1] = u
        lap = (ub[:-2] - 2.0 * ub[1:-1] + ub[2:]) / H ** 2
        return -(lap + lam * np.exp(u)) * H ** 2
    return F


def precondicionador_laplaciano():
    """M ≈ laplaciano discreto, factorizado UNA vez con LU disperso.

    El Jacobiano de Bratu es  J = L − h²λ diag(e^u),  donde L es el laplaciano
    discreto (tridiagonal). El término no lineal es una perturbación diagonal
    pequeña, así que L solo ya es una aproximación excelente de J. Y L es
    tridiagonal: su LU cuesta O(N) y se calcula una sola vez para toda la
    corrida.
    """
    L = sp.diags([-np.ones(N - 1), 2.0 * np.ones(N), -np.ones(N - 1)],
                 [-1, 0, 1], format="csc")
    lu = spla.splu(L)
    return LinearOperator((N, N), matvec=lu.solve, dtype=float)


# =============================================================================
# EL SOLVER
# =============================================================================
def newton_krylov(F, u0, forzado, precond=None, tol=1e-9, max_iter=50):
    """Newton inexacto + GMRES matriz-libre + búsqueda de línea de Armijo.

    `forzado` puede ser:
        un número  → η fijo en todas las iteraciones
        "ew1"      → Choice 1 de Eisenstat-Walker
        "ew2"      → Choice 2 de Eisenstat-Walker

    Devuelve un diccionario con el historial y los contadores de trabajo.
    """
    u = np.array(u0, dtype=float)
    contador = {"F": 0, "Jv": 0}

    def evaluar(x):
        contador["F"] += 1
        return F(x)

    Fu = evaluar(u)
    normaF = np.linalg.norm(Fu)
    historia = [normaF]
    etas_usados = []
    trabajo = [0]

    normaF_prev = None
    eta_prev = None
    res_lineal_prev = None
    gamma, alpha_ew = 0.9, (1.0 + math.sqrt(5.0)) / 2.0

    for k in range(max_iter):
        if normaF <= tol:
            break
        if not np.isfinite(normaF) or normaF > 1e14:
            break

        # ---------------- 1. elegir η ----------------
        if isinstance(forzado, (int, float)):
            eta = float(forzado)
        elif normaF_prev is None:
            eta = 0.9                       # primera vez: resolver flojo
        else:
            if forzado == "ew1" and res_lineal_prev is not None:
                eta = abs(normaF - res_lineal_prev) / normaF_prev
            else:                            # ew2 (y ew1 en el primer paso)
                eta = gamma * (normaF / normaF_prev) ** alpha_ew
            # Salvaguarda: no dejar que η caiga más rápido de lo razonable.
            if eta_prev is not None:
                piso = gamma * eta_prev ** alpha_ew
                if piso > 0.1:
                    eta = max(eta, piso)
            eta = min(0.9, max(eta, 0.5 * tol / normaF))
        etas_usados.append(eta)

        # ---------------- 2. el operador J, sin matriz ----------------
        def matvec(v, u=u, Fu=Fu):
            contador["Jv"] += 1
            nv = np.linalg.norm(v)
            if nv == 0.0:
                return np.zeros_like(v)
            xs = float(np.dot(u, v)) / nv
            e = SQRT_EPS
            if xs != 0.0:
                e = e * max(abs(xs), 1.0) * math.copysign(1.0, xs)
            e = e / nv
            return (evaluar(u + e * v) - Fu) / e

        Jop = LinearOperator((N, N), matvec=matvec, dtype=float)

        # ---------------- 3. resolver J s = −F solo hasta η ----------------
        s, info = gmres(Jop, -Fu, rtol=eta, atol=0.0, restart=40,
                        maxiter=200, M=precond)
        if info < 0 or not np.all(np.isfinite(s)):
            break
        res_lineal_prev = np.linalg.norm(Fu + Jop.matvec(s))

        # ---------------- 4. globalizar con Armijo ----------------
        lam = 1.0
        aceptado = False
        while lam > 1e-10:
            u_t = u + lam * s
            F_t = evaluar(u_t)
            n_t = np.linalg.norm(F_t)
            if n_t <= (1.0 - 1e-4 * lam) * normaF:
                aceptado = True
                break
            lam *= 0.5
        if not aceptado:
            break

        eta_prev = eta
        normaF_prev = normaF
        u, Fu, normaF = u_t, F_t, n_t
        historia.append(normaF)
        trabajo.append(contador["Jv"])

    return dict(u=u, convergio=normaF <= tol, normaF=normaF,
                iters=len(historia) - 1, n_F=contador["F"],
                n_Jv=contador["Jv"], historia=historia, etas=etas_usados,
                trabajo=trabajo)


# =============================================================================
# EXPERIMENTO 1 — Términos de forzado, con y sin precondicionador
# =============================================================================
CONFIGS = [(1e-1, "η = 1e-1 fijo"),
           (1e-3, "η = 1e-3 fijo"),
           (1e-8, "η = 1e-8 fijo"),
           ("ew1", "Eisenstat-Walker 1"),
           ("ew2", "Eisenstat-Walker 2")]


def experimento_1():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 1 · Cinco maneras de elegir η, con y sin precondicionador")
    print("=" * 78)

    lam = 3.0
    F = hacer_F(lam)
    M = precondicionador_laplaciano()
    u0 = np.zeros(N)

    print(f"""
  λ = {lam}   (el punto de retorno está en λ* ≈ {LAMBDA_CRITICO:.4f})
  N = {N} nodos interiores

  Se compara la MISMA corrida cambiando solo dos cosas: qué η se usa y si hay
  precondicionador. Lo que hay que mirar es la columna de productos J·v.
""")
    print(f"  {'forzado':>22s} | {'iters':>6s} | {'J·v sin M':>10s} | "
          f"{'iters':>6s} | {'J·v con M':>10s} | {'ganancia':>9s}")
    print("  " + "-" * 78)

    resultados = {}
    for forzado, nombre in CONFIGS:
        r_sin = newton_krylov(F, u0, forzado)
        r_con = newton_krylov(F, u0, forzado, precond=M)
        resultados[nombre] = (r_sin, r_con)
        ganancia = r_sin["n_Jv"] / max(r_con["n_Jv"], 1)
        print(f"  {nombre:>22s} | {r_sin['iters']:6d} | {r_sin['n_Jv']:10,d} | "
              f"{r_con['iters']:6d} | {r_con['n_Jv']:10,d} | {ganancia:8.1f}×")

    mejor = min(resultados.items(), key=lambda kv: kv[1][0]["n_Jv"])
    peor = max(resultados.items(), key=lambda kv: kv[1][0]["n_Jv"])
    print(f"""
  Sin precondicionador, el más barato es {mejor[0]} con {mejor[1][0]['n_Jv']:,} productos J·v
  y el más caro es {peor[0]} con {peor[1][0]['n_Jv']:,}. Son {peor[1][0]['n_Jv']/max(mejor[1][0]['n_Jv'],1):.1f}× de diferencia,
  y el caro es justamente el que MENOS iteraciones de Newton usa.

  Eso es el OVERSOLVING: resolver el sistema lineal con más precisión de la
  que el modelo lineal merece, sobre todo lejos de la raíz. Eisenstat-Walker
  lo evita sin que uno le diga nada: empieza flojo y aprieta solo cuando el
  residuo baja.
""")

    # --- figura ---
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
    colores = plt.cm.viridis(np.linspace(0, 0.85, len(CONFIGS)))
    for (nombre, (r_sin, r_con)), c in zip(resultados.items(), colores):
        a1.semilogy(np.maximum(r_sin["trabajo"], 1), r_sin["historia"], "o-",
                    color=c, ms=4, lw=1.5, label=nombre)
        a2.semilogy(np.maximum(r_con["trabajo"], 1), r_con["historia"], "o-",
                    color=c, ms=4, lw=1.5, label=nombre)
    for ax, t in [(a1, "Sin precondicionador"), (a2, "Con precondicionador")]:
        ax.set_xscale("log")
        ax.set_xlabel("productos J·v acumulados")
        ax.set_ylabel("‖F(u_k)‖")
        ax.set_title(t)
        ax.axhline(1e-9, color="0.5", ls=":", lw=1)
    a1.legend(fontsize=8.5)
    fig.suptitle("Lo que importa es el trabajo, no las iteraciones",
                 fontsize=13, fontweight="bold")
    fig.savefig("clase3_forzado.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  → figura guardada en clase3_forzado.png")
    return resultados


# =============================================================================
# EXPERIMENTO 2 — Cómo evoluciona η
# =============================================================================
def experimento_2(resultados):
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 2 · La trayectoria de η")
    print("=" * 78)

    r = resultados["Eisenstat-Walker 2"][1]
    print("""
  Eisenstat-Walker 2 usa   η_k = 0.9 · (‖F_k‖ / ‖F_{k−1}‖)^1.618 .

  O sea: mira cuánto bajó el residuo en el último paso. Si bajó mucho, el
  modelo lineal está funcionando y conviene apretar. Si bajó poco, no vale la
  pena gastar en resolver bien algo que no describe la realidad.
""")
    print(f"  {'k':>3s} | {'η elegido':>11s} | {'‖F_k‖':>11s} | {'razón ‖F_k‖/‖F_{k-1}‖':>22s}")
    print("  " + "-" * 60)
    for k, eta in enumerate(r["etas"]):
        razon = (r["historia"][k] / r["historia"][k - 1]) if k > 0 else float("nan")
        print(f"  {k:3d} | {eta:11.3e} | {r['historia'][k]:11.3e} | "
              f"{razon if k > 0 else 0:22.3e}")


# =============================================================================
# EXPERIMENTO 3 — Acercándose al punto de retorno
# =============================================================================
def experimento_3():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 3 · Qué pasa cuando λ se acerca a λ* ≈ 3.5138")
    print("=" * 78)

    print("""
  La ecuación de Bratu tiene solución solo para λ ≤ λ* ≈ 3.5138. En λ* las
  dos ramas de soluciones se juntan y el Jacobiano se vuelve SINGULAR.
  Acercarse a λ* es la manera limpia de ver qué le pasa a un método de Newton
  cuando el problema se vuelve mal condicionado.
""")
    M = precondicionador_laplaciano()
    print(f"  {'λ':>8s} | {'converge':>9s} | {'iters':>6s} | {'J·v':>8s} | "
          f"{'u máximo':>9s} | {'cond(J) estimado':>17s}")
    print("  " + "-" * 74)

    lams, conds, trabajos = [], [], []
    for lam in (1.0, 2.0, 3.0, 3.4, 3.50, 3.513, 3.6):
        F = hacer_F(lam)
        r = newton_krylov(F, np.zeros(N), "ew2", precond=M, max_iter=80)
        if r["convergio"]:
            J = sp.diags([-np.ones(N - 1),
                          2.0 - H ** 2 * lam * np.exp(r["u"]),
                          -np.ones(N - 1)], [-1, 0, 1], format="csc").toarray()
            cond = np.linalg.cond(J)
            umax = r["u"].max()
            lams.append(lam)
            conds.append(cond)
            trabajos.append(r["n_Jv"])
            print(f"  {lam:8.3f} | {'sí':>9s} | {r['iters']:6d} | {r['n_Jv']:8,d} | "
                  f"{umax:9.4f} | {cond:17.3e}")
        else:
            print(f"  {lam:8.3f} | {'NO':>9s} | {r['iters']:6d} | {r['n_Jv']:8,d} | "
                  f"{'—':>9s} | {'—':>17s}")

    print("""
  Para λ = 3.6 no hay solución que encontrar: el método no falla por ser malo,
  falla porque el problema no tiene respuesta. Es una distinción que conviene
  tener clara antes de culpar al solver.
""")

    if len(lams) >= 3:
        fig, ax = plt.subplots(figsize=(6.8, 4.2))
        ax.semilogy(lams, conds, "o-", color="#c0392b", ms=7, lw=2,
                    label="número de condición de J")
        ax.axvline(LAMBDA_CRITICO, color="black", ls="--", lw=1.5,
                   label=f"λ* ≈ {LAMBDA_CRITICO:.4f}")
        ax.set_xlabel("λ")
        ax.set_ylabel("cond(J) en la solución")
        ax.set_title("Al acercarse al punto de retorno,\nel Jacobiano se vuelve singular")
        ax.legend(fontsize=9)
        fig.savefig("clase3_punto_retorno.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
        print("  → figura guardada en clase3_punto_retorno.png")


# =============================================================================
# PREGUNTAS
# =============================================================================
def preguntas():
    print("\n" + "=" * 78)
    print("  PREGUNTAS PARA ANALIZAR")
    print("=" * 78)
    print("""
  1. En el Experimento 1, η = 1e-8 usa MENOS iteraciones de Newton que
     η = 1e-1 pero MÁS productos J·v. Explique por qué esas dos columnas no
     miden lo mismo y cuál de las dos es la que hay que optimizar.

  2. El precondicionador reduce el trabajo en un factor enorme. Explique por
     qué el laplaciano solo es tan buena aproximación de J en este problema
     (mire la expresión J = L − h²λ diag(e^u)) y diga en qué caso dejaría de
     serlo.

  3. Ponga λ = 3.5 y vuelva a correr el Experimento 1. ¿La ventaja de
     Eisenstat-Walker sobre η fijo crece o se achica? Proponga una explicación.

  4. Eisenstat-Walker 2 tiene la salvaguarda
         if piso > 0.1:  eta = max(eta, piso)
     con piso = 0.9·η_{k−1}^1.618. Coméntela (ponga la línea entre comillas)
     y vuelva a correr. ¿Cambia el número total de productos J·v? ¿Por qué
     hace falta esa salvaguarda?

  5. En el Experimento 3, para λ = 3.6 el método no converge. ¿Cómo
     distinguiría, viendo solo la salida de un programa, entre "el método
     falló" y "el problema no tiene solución"? Nombre al menos dos señales.

  6. Cambie N de 200 a 800 y corra otra vez el Experimento 1. El número de
     productos J·v SIN precondicionador debería crecer bastante, y CON
     precondicionador casi nada. Relacione esto con el hecho de que
     cond(L) = O(N²) y explique por qué precondicionar es lo que hace que el
     método escale a mallas finas.
""")


if __name__ == "__main__":
    print(__doc__)
    resultados = experimento_1()
    experimento_2(resultados)
    experimento_3()
    preguntas()
    print("=" * 78)
    print("  FIN DEL EJERCICIO 3")
    print("=" * 78)
