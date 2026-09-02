# -*- coding: utf-8 -*-
"""
=============================================================================
 EJERCICIO 1 PARA LA CLASE — Búsqueda de línea de Armijo: qué garantiza y qué no
=============================================================================
 Materia : Métodos Numéricos II (DAT-252) — UMSA
 Tema    : Estrategias para la convergencia global · Métodos de Newton-Krylov

 QUÉ HACE ESTE PROGRAMA
   Resuelve dos sistemas no lineales de 2×2 con el método de Newton, en dos
   versiones: con el paso completo (Newton puro) y con búsqueda de línea de
   Armijo. Después dibuja las cuencas de convergencia de ambos.

 LO IMPORTANTE
   No es un ejercicio para confirmar que "globalizar siempre es mejor".
   Es un ejercicio para ver EXACTAMENTE qué promete la búsqueda de línea:

       promete que ‖F‖ nunca aumenta,
       NO promete llegar a una raíz.

   El Problema A muestra el lado bueno. El Problema B muestra el lado malo,
   y es el más instructivo de los dos.

 REQUISITOS
   pip install numpy matplotlib

 EJECUTAR
   python3 clase1_armijo.py
=============================================================================
"""

import warnings

import matplotlib
matplotlib.use("Agg")                       # guardar figuras sin abrir ventana
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# LOS DOS PROBLEMAS
# =============================================================================
def F_A(v):
    """Problema A — circunferencia contra exponencial empinada.

        F1 = x1² + x2² − 4
        F2 = e^(3 x1) + x2 − 1

    El e^(3x1) crece muy rápido: si Newton da un paso largo hacia la derecha,
    F2 se dispara y el residuo desborda.
    """
    return np.array([v[0] ** 2 + v[1] ** 2 - 4.0,
                     np.exp(3.0 * v[0]) + v[1] - 1.0])


def J_A(v):
    return np.array([[2.0 * v[0], 2.0 * v[1]],
                     [3.0 * np.exp(3.0 * v[0]), 1.0]])


def F_B(v):
    """Problema B — función de Freudenstein y Roth (1963).

        F1 = −13 + x1 + ((5 − x2) x2 − 2) x2
        F2 = −29 + x1 + ((x2 + 1) x2 − 14) x2

    Raíz única: x* = (5, 4).

    Es un caso clásico de la literatura de optimización porque la función de
    mérito f(x) = ½‖F(x)‖² tiene un MÍNIMO LOCAL que no es raíz, situado
    aproximadamente en x2 ≈ −0.8968, con f ≈ 48.98 ≠ 0. Cualquier método que
    solo sepa "bajar" cae ahí y se queda.
    """
    return np.array([-13.0 + v[0] + ((5.0 - v[1]) * v[1] - 2.0) * v[1],
                     -29.0 + v[0] + ((v[1] + 1.0) * v[1] - 14.0) * v[1]])


def J_B(v):
    return np.array([[1.0, (10.0 - 3.0 * v[1]) * v[1] - 2.0],
                     [1.0, (3.0 * v[1] + 2.0) * v[1] - 14.0]])


# =============================================================================
# EL MÉTODO
# =============================================================================
def newton(F, J, x0, usar_armijo, max_iter=40, tol=1e-10,
           alpha=1e-4, lam_min=1e-10, guardar=False):
    """Newton para sistemas 2×2, con o sin búsqueda de línea.

    Devuelve un diccionario con el resultado y, si se pide, la trayectoria.

    ---- LA PARTE QUE HAY QUE LEER CON CUIDADO ----

    El paso de Newton s resuelve   J(x) s = −F(x).
    La función de mérito es        f(x) = ½ ‖F(x)‖² .
    Su gradiente es                ∇f(x) = J(x)ᵀ F(x) .

    Entonces, a lo largo de la dirección de Newton:

        ∇f(x)ᵀ s = F(x)ᵀ J(x) s = F(x)ᵀ (−F(x)) = −‖F(x)‖² < 0

    O sea: el paso de Newton SIEMPRE apunta cuesta abajo para f. Por lo tanto
    existe algún λ > 0 lo bastante chico que hace bajar a f. Buscar ese λ es
    justamente la búsqueda de línea. La condición que se le exige es la
    CONDICIÓN DE ARMIJO, escrita aquí en términos de la norma del residuo:

        ‖F(x + λ s)‖ ≤ (1 − α λ) ‖F(x)‖ ,     α = 1e-4

    El α es minúsculo a propósito: se pide una mejora casi simbólica, solo
    para impedir que la sucesión se estanque bajando cantidades cada vez más
    pequeñas sin llegar a nada.
    """
    x = np.array(x0, dtype=float)
    Fx = F(x)
    normaF = np.linalg.norm(Fx)
    traza = [x.copy()]
    lambdas = []
    residuos = [normaF]

    for k in range(max_iter):
        if normaF <= tol:
            return dict(estado="raiz", iters=k, x=x, normaF=normaF,
                        traza=traza, lambdas=lambdas, residuos=residuos)
        if not np.isfinite(normaF) or normaF > 1e12:
            return dict(estado="desbordo", iters=k, x=x, normaF=normaF,
                        traza=traza, lambdas=lambdas, residuos=residuos)

        Jx = J(x)
        if abs(np.linalg.det(Jx)) < 1e-14:
            return dict(estado="J singular", iters=k, x=x, normaF=normaF,
                        traza=traza, lambdas=lambdas, residuos=residuos)
        s = np.linalg.solve(Jx, -Fx)

        if not usar_armijo:
            # ---- Newton puro: se acepta el paso completo, pase lo que pase ----
            lam = 1.0
            x = x + s
            Fx = F(x)
            normaF = np.linalg.norm(Fx)
        else:
            # ---- Búsqueda de línea: se prueba λ = 1, ½, ¼, ⅛, ... ----
            lam = 1.0
            aceptado = False
            while lam >= lam_min:
                x_t = x + lam * s
                F_t = F(x_t)
                n_t = np.linalg.norm(F_t)
                # OJO con esta comparación. Está escrita como
                #     if n_t <= cota:  aceptar
                # y NO como
                #     if not (n_t > cota):  aceptar
                # que parecen lo mismo pero no lo son. Si el paso completo
                # desborda, n_t vale nan, y en Python "nan > cota" es False:
                # la segunda versión ACEPTARÍA un paso inválido. La primera
                # da False para nan y el retroceso continúa, que es lo
                # correcto. Es la pregunta 4 del final.
                if n_t <= (1.0 - alpha * lam) * normaF:
                    aceptado = True
                    break
                lam = lam * 0.5
            if not aceptado:
                # λ se hizo tan chico que ya no tiene sentido seguir. Esto NO
                # significa que estemos en una raíz: significa que en esta
                # dirección ya no se puede bajar.
                return dict(estado="linea fallo", iters=k, x=x, normaF=normaF,
                            traza=traza, lambdas=lambdas, residuos=residuos)
            x, Fx, normaF = x_t, F_t, n_t

        lambdas.append(lam)
        residuos.append(normaF)
        if guardar:
            traza.append(x.copy())

    return dict(estado="max iters", iters=max_iter, x=x, normaF=normaF,
                traza=traza, lambdas=lambdas, residuos=residuos)


# =============================================================================
# EXPERIMENTO 1 — Trayectorias desde puntos concretos
# =============================================================================
def experimento_1():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 1 · Qué pasa desde cada punto inicial")
    print("=" * 78)

    casos = [("A", F_A, J_A, [(1.0, 1.5), (-0.6, -1.5), (-0.5, -0.9), (-0.7, -2.2)]),
             ("B", F_B, J_B, [(6.0, 3.0), (0.5, -2.0), (15.0, -2.0), (0.0, 0.0)])]

    for letra, F, J, x0s in casos:
        print(f"\n  --- PROBLEMA {letra} ---")
        print(f"  {'x0':>14s} | {'Newton puro':>26s} | {'Newton + Armijo':>26s}")
        print("  " + "-" * 72)
        for x0 in x0s:
            rp = newton(F, J, x0, usar_armijo=False)
            ra = newton(F, J, x0, usar_armijo=True)
            def resumen(r):
                if r["estado"] == "raiz":
                    return f"raíz ({r['x'][0]:+.3f},{r['x'][1]:+.3f}) {r['iters']:2d}it"
                return f"{r['estado']:12s} ‖F‖={r['normaF']:.2e}"
            print(f"  {str(x0):>14s} | {resumen(rp):>26s} | {resumen(ra):>26s}")

    print("""
  En el Problema B, mírense los ‖F‖ finales de la búsqueda de línea: no son
  cero, son del orden de 7 u 8. Y sin embargo el método se detuvo "porque no
  podía bajar más". Eso es el mínimo local de f, no una raíz.
""")


# =============================================================================
# EXPERIMENTO 2 — Cuencas de convergencia
# =============================================================================
def experimento_2():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 2 · Cuencas de convergencia (41×41 = 1681 puntos)")
    print("=" * 78)

    config = [("A", F_A, J_A, (-2.5, 1.0, -2.5, 2.5)),
              ("B", F_B, J_B, (-10.0, 20.0, -8.0, 8.0))]

    fig, ejes = plt.subplots(2, 2, figsize=(11.5, 8.5))
    for fila, (letra, F, J, caja) in enumerate(config):
        a, b, c, d = caja
        g1 = np.linspace(a, b, 41)
        g2 = np.linspace(c, d, 41)
        print(f"\n  --- PROBLEMA {letra} ---")
        for col, usar_armijo in enumerate((False, True)):
            M = np.full((41, 41), np.nan)
            for i, q in enumerate(g2):
                for j, p in enumerate(g1):
                    r = newton(F, J, (p, q), usar_armijo)
                    if r["estado"] == "raiz":
                        M[i, j] = r["iters"]
            pct = 100 * np.count_nonzero(~np.isnan(M)) / M.size
            nombre = "Newton + Armijo" if usar_armijo else "Newton puro"
            print(f"    {nombre:18s}: llega a una raíz desde {pct:5.1f} % de los puntos")

            ax = ejes[fila, col]
            im = ax.imshow(M, origin="lower", cmap="viridis_r", aspect="auto",
                           extent=[a, b, c, d], vmin=0, vmax=25)
            ax.set_facecolor("#7b241c")
            ax.set_title(f"Problema {letra} — {nombre}\n{pct:.1f} % converge",
                         fontsize=11)
            ax.set_xlabel("x₁ inicial")
            ax.set_ylabel("x₂ inicial")
            fig.colorbar(im, ax=ax, label="iteraciones")

    fig.suptitle("Rojo oscuro = NO llega a ninguna raíz", fontsize=13,
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig("clase1_cuencas.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("\n  → figura guardada en clase1_cuencas.png")


# =============================================================================
# EXPERIMENTO 3 — ¿Y si cambiamos α?
# =============================================================================
def experimento_3():
    print("\n" + "=" * 78)
    print("  EXPERIMENTO 3 · El parámetro α de la condición de Armijo")
    print("=" * 78)

    print("""
  La condición es  ‖F(x + λs)‖ ≤ (1 − αλ)‖F(x)‖.
  α mide cuánta mejora se exige a cambio de aceptar el paso.
""")
    print(f"  {'α':>10s} | {'iteraciones':>12s} | {'‖F‖ final':>12s} | {'λ promedio':>11s} | estado")
    print("  " + "-" * 68)
    for alpha in (1e-8, 1e-4, 1e-2, 0.1, 0.5, 0.9):
        r = newton(F_A, J_A, (1.0, 1.5), usar_armijo=True, alpha=alpha)
        lam_prom = np.mean(r["lambdas"]) if r["lambdas"] else float("nan")
        print(f"  {alpha:10.0e} | {r['iters']:12d} | {r['normaF']:12.2e} | "
              f"{lam_prom:11.3f} | {r['estado']}")

    print("""
  Con α grande se exige tanta mejora que casi ningún λ la cumple: el método
  retrocede mucho, avanza poco y necesita más iteraciones. Por eso el valor
  estándar en la literatura es α = 1e-4: apenas lo justo para que el teorema
  de convergencia funcione, sin frenar al método.
""")


# =============================================================================
# PREGUNTAS
# =============================================================================
def preguntas():
    print("\n" + "=" * 78)
    print("  PREGUNTAS PARA ANALIZAR")
    print("=" * 78)
    print("""
  1. En el Problema A, la búsqueda de línea convierte varios "desbordo" en
     raíces encontradas. Explique, mirando F2 = e^(3x1) + x2 − 1, por qué el
     paso completo de Newton desborda y por qué medio paso no lo hace.

  2. En el Problema B pasa lo contrario: Newton puro encuentra la raíz desde
     muchos más puntos que Newton con Armijo. Esto NO es un error del
     programa. Explique por qué, usando estas dos ideas:
         · la búsqueda de línea obliga a que ‖F‖ baje en cada paso;
         · f(x) = ½‖F(x)‖² tiene un mínimo local en x2 ≈ −0.8968 donde f ≈ 49.
     ¿Qué le pasa a una sucesión que solo sabe bajar cuando cae en ese valle?

  3. Enuncie con sus palabras qué garantiza exactamente el teorema de
     convergencia global de la búsqueda de línea. Complete la frase:
         "Si F es continuamente diferenciable y los iterados permanecen en un
          conjunto acotado, la búsqueda de línea garantiza que la sucesión
          converge a __________________ de f, que puede o no ser una raíz."

  4. En el código, la condición de aceptación está escrita como
         if n_t <= (1.0 - alpha * lam) * normaF:
     y no como
         if not (n_t > (1.0 - alpha * lam) * normaF):
     Suponga que el paso completo produce n_t = nan. Verifique cuál de las dos
     versiones aceptaría el paso inválido. (Pruebe en el intérprete:
     `float('nan') > 1` y `float('nan') <= 1`.)

  5. Ponga lam_min = 1e-1 en la función newton y vuelva a correr el
     Experimento 2. ¿Cómo cambia la cuenca del Problema A? ¿Por qué?

  6. La búsqueda de línea solo se mueve a lo largo de la dirección de Newton.
     La región de confianza, en cambio, puede elegir OTRA dirección cuando el
     paso de Newton no es de fiar. A la vista del Problema B, ¿por qué esa
     libertad podría ayudar? ¿Y por qué tampoco es una garantía?
""")


if __name__ == "__main__":
    print(__doc__)
    experimento_1()
    experimento_2()
    experimento_3()
    preguntas()
    print("=" * 78)
    print("  FIN DEL EJERCICIO 1")
    print("=" * 78)
