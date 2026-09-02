# -*- coding: utf-8 -*-
"""
=============================================================================
 EJERCICIO 2 PARA LA CLASE — El Jacobiano que nunca se construye
=============================================================================
 Materia : Métodos Numéricos II (DAT-252) — UMSA
 Tema    : Estrategias para la convergencia global · Métodos de Newton-Krylov

 LA IDEA CENTRAL
   Newton necesita resolver  J(x) s = −F(x).  Si el sistema es grande, formar
   J cuesta n evaluaciones de F, guardarla cuesta n² números y factorizarla
   cuesta O(n³). Para n = 10⁶ eso es imposible.

   Pero los métodos de Krylov (GMRES, BiCGSTAB, ...) NO necesitan la matriz:
   solo necesitan poder calcular productos J·v. Y ese producto se puede
   aproximar con UNA evaluación de F:

        J(x)·v  ≈  [ F(x + ε v) − F(x) ] / ε

   A eso se le llama "Jacobian-free Newton-Krylov" (JFNK). Este programa
   estudia las dos preguntas que decide todo:
        1. ¿Qué ε hay que usar?
        2. ¿Cuán bien hay que resolver el sistema lineal?

 REQUISITOS
   pip install numpy scipy matplotlib

 EJECUTAR
   python3 clase2_matrix_free.py
=============================================================================
"""

import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, gmres

warnings.filterwarnings("ignore")

# --- Problema de prueba: reacción-difusión en 1D -------------------------
#     −ν u'' + u³ = f(x),   u(0) = u(1) = 0,   ν = 0.01
# Elegimos f de modo que la solución exacta sea u(x) = sin(πx), así podemos
# medir el error de verdad y no solo el residuo.
#
# Atención a un detalle que el Experimento 1 va a poner en evidencia: en este
# problema el término difusivo ν/h² ≈ 900 pesa MUCHO más que el no lineal
# u³ ≈ 1. O sea que F es casi lineal, y su segunda derivada es pequeña
# comparada con la propia F. Eso corre el ε óptimo hacia arriba respecto de la
# receta clásica √eps_maq, que supone F, F' y F'' del mismo orden. Lo medimos
# y lo discutimos ahí.
N = 300
NU = 0.01
H = 1.0 / (N + 1)
X = np.linspace(0.0, 1.0, N + 2)[1:-1]
U_EXACTA = np.sin(np.pi * X)
FUENTE = NU * np.pi ** 2 * np.sin(np.pi * X) + np.sin(np.pi * X) ** 3


def F(u):
    """Residuo discreto de −ν u'' + u³ − f."""
    ub = np.zeros(N + 2)
    ub[1:-1] = u
    lap = (ub[:-2] - 2.0 * ub[1:-1] + ub[2:]) / H ** 2
    return -NU * lap + u ** 3 - FUENTE


def J_analitico(u):
    """Jacobiano exacto y disperso: tridiagonal.

    ∂F_i/∂u_i     =  2ν/h² + 3 u_i²
    ∂F_i/∂u_{i±1} = −ν/h²
    Solo lo usamos para MEDIR el error de la versión sin matriz. Un programa
    JFNK de verdad nunca lo escribiría.
    """
    fuera = -NU / H ** 2 * np.ones(N - 1)
    return sp.diags([fuera,
                     2.0 * NU / H ** 2 + 3.0 * u ** 2,
                     fuera], [-1, 0, 1], format="csr")


def Jv_diferencias(u, Fu, v, eps):
    """El producto J·v con UNA evaluación de F. Este es el corazón del método."""
    return (F(u + eps * v) - Fu) / eps


# =============================================================================
# EXPERIMENTO 1 — ¿Qué ε elegir?
# =============================================================================
def experimento_1():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 1 · El paso ε de la diferencia finita")
    print("=" * 78)

    print("""
  Hay dos errores tirando en direcciones opuestas:

    · Error de TRUNCAMIENTO ≈ (ε/2)·‖F''‖·‖v‖². La fórmula [F(x+εv) − F(x)]/ε
      es la derivada más un término O(ε). Baja cuando ε baja.

    · Error de REDONDEO (cancelación) ≈ eps_maq·‖F‖/(ε‖v‖). F(x+εv) y F(x) se
      parecen cada vez más; al restarlos se pierden cifras significativas y
      luego se divide entre un número diminuto, que amplifica el desastre.
      Crece cuando ε baja.

  Igualando ambos:      ε_óptimo ≈ √( eps_maq · ‖F‖ / ‖F''‖ )

  De ahí sale la receta famosa ε ≈ √eps_maq ≈ 1.5e-8: es el caso en que ‖F‖ y
  ‖F''‖ son del mismo orden. En ESTE problema no lo son (la difusión domina y
  F es casi lineal, F'' es chica), así que el óptimo medido queda por encima.
  Vale la pena ver las dos cosas: dónde está el óptimo de verdad, y cuánto se
  pierde por usar la receta de todos modos.
""")

    u = 0.5 * np.sin(np.pi * X) + 0.2 * np.sin(3 * np.pi * X)
    Fu = F(u)
    rng = np.random.default_rng(42)
    v = rng.standard_normal(N)
    v = v / np.linalg.norm(v)
    exacto = J_analitico(u) @ v

    eps_maq = np.finfo(float).eps
    epsilons = np.logspace(-16, -1, 76)
    errores = np.array([np.linalg.norm(Jv_diferencias(u, Fu, v, e) - exacto)
                        / np.linalg.norm(exacto) for e in epsilons])

    i_mejor = int(np.argmin(errores))
    eps_opt = epsilons[i_mejor]
    err_receta = (np.linalg.norm(Jv_diferencias(u, Fu, v, np.sqrt(eps_maq)) - exacto)
                  / np.linalg.norm(exacto))
    print(f"  eps de la máquina              : {eps_maq:.3e}")
    print(f"  √eps_maq  (la receta)          : {np.sqrt(eps_maq):.3e}")
    print(f"  ε que minimiza el error MEDIDO : {eps_opt:.3e}")
    print(f"  error mínimo alcanzable        : {errores[i_mejor]:.3e}")
    print(f"  error usando la receta √eps_maq: {err_receta:.3e}")
    print()
    print(f"  {'ε':>10s} | {'error relativo':>15s} | quién manda")
    print("  " + "-" * 50)
    for e in (1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2):
        err = np.linalg.norm(Jv_diferencias(u, Fu, v, e) - exacto) / np.linalg.norm(exacto)
        quien = "cancelación" if e < eps_opt else "truncamiento"
        print(f"  {e:10.0e} | {err:15.3e} | {quien}")
    print(f"""
  Léase la última línea de la cabecera: la receta √eps_maq da un error de
  {err_receta:.1e}, unas {err_receta/errores[i_mejor]:.0f} veces peor que el óptimo… y da igual.
  Newton inexacto ya tolera un residuo lineal de η‖F‖ con η ~ 1e-2, así que
  un Jacobiano con 7 cifras correctas le sobra. Por eso √eps_maq se usa
  siempre: no es óptima, pero nunca es catastrófica y no exige conocer ‖F''‖.""")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(epsilons, errores, "o-", color="#1f4e79", ms=3, lw=1.4,
              label="error medido")
    ax.loglog(epsilons, epsilons * 5, "--", color="#e67e22", lw=1.3,
              label="modelo del truncamiento  ~ ε")
    ax.loglog(epsilons, eps_maq / epsilons * 0.5, "--", color="#c0392b", lw=1.3,
              label="modelo de la cancelación ~ eps_maq/ε")
    ax.axvline(np.sqrt(eps_maq), color="black", ls=":", lw=1.5,
               label=r"$\sqrt{\epsilon_{maq}}$")
    ax.set_ylim(1e-12, 1e3)
    ax.set_xlabel("paso ε")
    ax.set_ylabel("error relativo en J·v")
    ax.set_title("El ε óptimo no es 'el más chico posible'")
    ax.legend(fontsize=8.5)
    fig.savefig("clase2_epsilon.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("\n  → figura guardada en clase2_epsilon.png")


# =============================================================================
# EXPERIMENTO 2 — Costo: matriz explícita contra matriz-libre
# =============================================================================
def experimento_2():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 2 · ¿Cuánto cuesta formar J?")
    print("=" * 78)

    u = 0.5 * np.sin(np.pi * X)
    Fu = F(u)

    # Formar J columna a columna por diferencias finitas: n evaluaciones de F.
    eps = np.sqrt(np.finfo(float).eps)
    base = np.eye(N)
    J_fd = np.column_stack([Jv_diferencias(u, Fu, base[:, j], eps) for j in range(N)])
    evals_para_formar = N

    print(f"""
  Dimensión del sistema                      : n = {N}
  Evaluaciones de F para formar J por columnas: {evals_para_formar}
  Números que hay que guardar (matriz densa)  : {N*N:,}
  Evaluaciones de F para UN producto J·v      : 1

  Y aquí está el punto: GMRES resuelve el sistema con unas pocas decenas de
  productos J·v. Formar la matriz para después resolverla cuesta {N} — es
  decir, unas {N//30} veces más caro, y encima hay que guardarla.

  Error de la J formada por diferencias vs. la analítica: {np.linalg.norm(J_fd - J_analitico(u).toarray()) / np.linalg.norm(J_analitico(u).toarray()):.3e}
""")


# =============================================================================
# EXPERIMENTO 3 — El oversolving
# =============================================================================
def experimento_3():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 3 · ¿Cuán bien hay que resolver J s = −F?")
    print("=" * 78)

    print("""
  Newton INEXACTO acepta cualquier s que cumpla

        ‖J s + F‖ ≤ η ‖F‖ .

  El parámetro η se llama TÉRMINO DE FORZADO. El teorema de Dembo, Eisenstat
  y Steihaug (1982) dice:

        η_k ≤ η < 1 constante  →  convergencia lineal
        η_k → 0                →  convergencia superlineal
        η_k = O(‖F_k‖)         →  convergencia cuadrática (Newton de verdad)

  Suena a que conviene η lo más chico posible. Pero cada dígito extra en el
  sistema lineal cuesta iteraciones de GMRES, o sea productos J·v, o sea
  evaluaciones de F. Y lejos de la raíz el modelo lineal ni siquiera describe
  bien a F, así que ese esfuerzo se desperdicia. Eso es el OVERSOLVING.

  Contamos el trabajo REAL: número total de productos J·v.
""")

    print(f"  {'η':>10s} | {'iters. Newton':>14s} | {'productos J·v':>14s} | "
          f"{'‖F‖ final':>11s} | {'error vs exacta':>15s}")
    print("  " + "-" * 78)

    resultados = []
    for eta in (0.5, 1e-1, 1e-2, 1e-4, 1e-8, 1e-12):
        u = np.zeros(N)
        Fu = F(u)
        n_jv = 0
        historia = [np.linalg.norm(Fu)]

        for k in range(60):
            if np.linalg.norm(Fu) <= 1e-10:
                break

            # --- El operador lineal que GMRES verá. La matriz no existe. ---
            def matvec(v, u=u, Fu=Fu):
                nonlocal n_jv
                n_jv += 1
                nv = np.linalg.norm(v)
                if nv == 0.0:
                    return np.zeros_like(v)
                # ε escalado con la magnitud de u: la receta de Kelley.
                e = np.sqrt(np.finfo(float).eps) * max(np.linalg.norm(u), 1.0) / nv
                return Jv_diferencias(u, Fu, v, e)

            Jop = LinearOperator((N, N), matvec=matvec, dtype=float)
            s, info = gmres(Jop, -Fu, rtol=eta, atol=0.0, restart=40, maxiter=200)

            # Búsqueda de línea de Armijo, para que la comparación sea justa:
            # lo único que cambia entre corridas es η.
            lam, nF = 1.0, np.linalg.norm(Fu)
            while lam > 1e-10:
                u_t = u + lam * s
                F_t = F(u_t)
                if np.linalg.norm(F_t) <= (1 - 1e-4 * lam) * nF:
                    break
                lam *= 0.5
            u, Fu = u_t, F_t
            historia.append(np.linalg.norm(Fu))

        err = np.linalg.norm(u - U_EXACTA, np.inf)
        print(f"  {eta:10.0e} | {len(historia)-1:14d} | {n_jv:14,d} | "
              f"{np.linalg.norm(Fu):11.2e} | {err:15.3e}")
        resultados.append((eta, len(historia) - 1, n_jv, historia))

    barato = min(resultados, key=lambda r: r[2])
    caro = max(resultados, key=lambda r: r[2])
    print(f"""
  Más barato : η = {barato[0]:.0e}  con {barato[2]:,} productos J·v
  Más caro   : η = {caro[0]:.0e}  con {caro[2]:,} productos J·v   ({caro[2]/barato[2]:.1f}× más)

  Nótese que el más caro NO es el que más iteraciones de Newton usa. Contar
  iteraciones externas engaña: mide lo que se ve, no lo que se paga.

  El error contra la solución exacta es prácticamente el mismo en todas las
  filas. O sea: pagar por resolver el sistema lineal con 12 cifras no compró
  ni una cifra más de precisión en la respuesta final.
""")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    colores = plt.cm.viridis(np.linspace(0, 0.9, len(resultados)))
    for (eta, its, njv, hist), c in zip(resultados, colores):
        a1.semilogy(hist, "o-", color=c, ms=4, lw=1.5, label=f"η = {eta:.0e}")
    a1.set_xlabel("iteración de Newton k")
    a1.set_ylabel("‖F(u_k)‖")
    a1.set_title("Por iteración: η chico parece ganar")
    a1.legend(fontsize=8)

    etas = [r[0] for r in resultados]
    njvs = [r[2] for r in resultados]
    a2.semilogx(etas, njvs, "o-", color="#c0392b", ms=8, lw=2)
    a2.set_xlabel("término de forzado η")
    a2.set_ylabel("productos J·v totales")
    a2.set_title("Por trabajo real: η chico se paga caro")
    a2.invert_xaxis()
    fig.savefig("clase2_oversolving.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  → figura guardada en clase2_oversolving.png")


# =============================================================================
# PREGUNTAS
# =============================================================================
def preguntas():
    print("\n" + "=" * 78)
    print("  PREGUNTAS PARA ANALIZAR")
    print("=" * 78)
    print("""
  1. En el Experimento 1, el error con ε = 1e-16 es PEOR que con ε = 1e-8,
     aunque 1e-16 sea "más pequeño". Explique por qué, y calcule con qué ε se
     igualan los dos modelos ε·C₁ y eps_maq/ε·C₂ (tome C₁ = C₂ = 1).

  1b. El ε óptimo MEDIDO no coincide con √eps_maq: queda unos dos órdenes por
     encima. Usando ε_óptimo ≈ √(eps_maq·‖F‖/‖F''‖) y sabiendo que en este
     problema la parte lineal vale ν/h² ≈ 900 mientras que F'' = 6u ≈ 6,
     estime ε_óptimo a mano y compárelo con el valor que imprime el programa.
     ¿Justifica esto abandonar la receta √eps_maq en la práctica?

  2. En la función matvec del Experimento 3, el paso es
         e = sqrt(eps_maq) * max(‖u‖, 1) / ‖v‖
     y no simplemente e = sqrt(eps_maq). ¿Qué pasaría si u tuviera componentes
     del orden de 10⁶ y usáramos ε fijo? ¿Y si v tuviera norma 10⁻³?

  3. El Experimento 2 muestra que formar J cuesta n = 300 evaluaciones de F.
     Suponga ahora un problema 3D con malla 100×100×100. ¿Cuántas
     evaluaciones de F costaría formar J? ¿Cuánta memoria ocuparía la matriz
     densa en gigabytes (8 bytes por número)?

  4. Ponga N = 1000 al principio del archivo y vuelva a correr. ¿Cómo cambia
     el número de productos J·v que necesita cada η? Relacione la respuesta
     con el número de condición del laplaciano discreto, que crece como 1/h².

  5. En el Experimento 3, ¿cuál sería la estrategia ideal para η? Debería ser
     grande al principio (cuando el modelo lineal no vale) y chica al final
     (para recuperar la convergencia cuadrática). Escriba una fórmula que haga
     eso automáticamente en función de ‖F_k‖ y ‖F_{k−1}‖. Compárela después
     con la de Eisenstat-Walker:  η_k = 0.9 (‖F_k‖/‖F_{k−1}‖)^1.618 .

  6. ¿Por qué el método se llama "matriz-libre" si en el Experimento 2 sí
     construimos la matriz? Identifique en el código qué parte es realmente
     JFNK y qué parte está ahí solo para medir.
""")


if __name__ == "__main__":
    print(__doc__)
    experimento_1()
    experimento_2()
    experimento_3()
    preguntas()
    print("=" * 78)
    print("  FIN DEL EJERCICIO 2")
    print("=" * 78)
