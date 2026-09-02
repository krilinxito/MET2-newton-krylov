# -*- coding: utf-8 -*-
"""
EJERCICIO 1 (exposición) — ¿Por qué hace falta globalizar?

Materia: Métodos Numéricos II (DAT-252) — UMSA
Tema:    Estrategias para la convergencia global · Métodos de Newton-Krylov

OBJETIVO
    Mostrar, con las manos, que el método de Newton es cuadráticamente
    convergente pero SOLO LOCALMENTE, y que una estrategia de globalización
    barata (búsqueda de línea o región de confianza) convierte un método que
    diverge en uno que converge.

QUÉ MIRAR DURANTE LA EXPOSICIÓN
    Parte A — arctan(x) = 0 desde x0 = 2. El paso de Newton se dispara porque
              la tangente cruza el eje muy lejos. La sucesión oscila y crece.
    Parte B — Sistema 2×2 de Kelley desde x0 = (2, 0.5). Newton puro desborda;
              Armijo y la región de confianza llegan a la raíz (1,1).
    Parte C — Mapa de cuencas de convergencia: el mismo problema, 2601 puntos
              iniciales. Se ve de un vistazo cuánto se agranda la cuenca.
    Parte D — La letra chica: la globalización tampoco es magia. Si la
              iteración cae sobre la variedad donde J es singular (aquí
              x2 = 0), ninguna búsqueda de línea la salva.

EJECUTAR
    python3 ej1_newton_vs_globalizado.py
"""

import warnings

import numpy as np

import nk_lib
from nk_lib import estilo_figuras, tabla, titulo

warnings.filterwarnings("ignore", category=RuntimeWarning)
plt = estilo_figuras()
SALIDA = "figuras"


# =============================================================================
# PARTE A — El caso más simple posible: arctan(x) = 0
# =============================================================================
def parte_A():
    titulo("PARTE A · arctan(x) = 0 — Newton se dispara desde x0 = 2")

    f = np.arctan
    fp = lambda x: 1.0 / (1.0 + x * x)

    print("""
  La raíz es x* = 0 y arctan es suave y monótona: parece inofensiva.
  Pero arctan se aplana lejos del origen, así que la tangente cruza el eje
  MUY lejos y el paso de Newton x - f(x)/f'(x) sobrepasa la raíz cada vez
  más. El umbral es |x0| ≈ 1.3917: por encima, Newton puro diverge.
""")

    filas = []
    for x0 in (1.0, 1.3, 1.5, 2.0):
        x = x0
        for _ in range(12):
            x = x - f(x) / fp(x)
            if not np.isfinite(x) or abs(x) > 1e12:
                break
        estado = "converge" if abs(x) < 1e-8 else "DIVERGE"
        filas.append([f"{x0:.2f}", f"{x:+.4e}" if np.isfinite(x) else "overflow", estado])
    tabla(filas, ["x0", "x tras 12 pasos", "resultado"])

    # --- Iteraciones desde x0 = 2, con y sin búsqueda de línea ---
    def F(v):
        return np.arctan(v)

    h_puro = nk_lib.newton_krylov(F, np.array([2.0]), globalizacion="ninguna",
                                  forzado=1e-12, tol=1e-12, max_iter=8,
                                  guardar_trayectoria=True)
    h_line = nk_lib.newton_krylov(F, np.array([2.0]), globalizacion="linea",
                                  forzado=1e-12, tol=1e-12, max_iter=8,
                                  guardar_trayectoria=True)
    it_puro = [p[0] for p in h_puro.trayectoria]
    it_line = [p[0] for p in h_line.trayectoria]

    print(f"\n  Newton puro  desde x0=2 : {['%.3g' % v for v in it_puro[:6]]} ...")
    print(f"  Newton+Armijo desde x0=2 : {['%.3g' % v for v in it_line[:6]]} ...")
    print(f"  λ usados por Armijo      : {['%.3f' % v for v in h_line.lambdas]}")
    print(f"\n  Newton puro : {'convergió' if h_puro.convergio else 'NO convergió'}")
    print(f"  Con Armijo  : {'convergió' if h_line.convergio else 'NO convergió'}"
          f" en {h_line.n_newton} iteraciones")

    # --- Figura: la construcción geométrica del disparo ---
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    xs = np.linspace(-30, 30, 800)
    a1.plot(xs, np.arctan(xs), color="#1f4e79", lw=2, label="arctan(x)")
    a1.axhline(0, color="0.4", lw=0.8)
    colores = ["#c0392b", "#e67e22", "#8e44ad"]
    for i in range(3):
        xk = it_puro[i]
        xk1 = it_puro[i + 1]
        a1.plot([xk, xk], [0, f(xk)], color=colores[i], ls=":", lw=1.2)
        a1.plot([xk, xk1], [f(xk), 0], color=colores[i], lw=1.6,
                label=f"tangente en x{i} = {xk:.2f}")
        a1.plot([xk], [f(xk)], "o", color=colores[i], ms=5)
    a1.set_xlim(-30, 30)
    a1.set_ylim(-1.8, 1.8)
    a1.set_title("Newton puro: cada tangente cruza más lejos")
    a1.set_xlabel("x")
    a1.legend(fontsize=8, loc="lower right")

    a2.semilogy(np.abs(np.array(it_puro)) + 1e-18, "o-", color="#c0392b",
                label="Newton puro (paso completo)")
    a2.semilogy(np.abs(np.array(it_line)) + 1e-18, "s-", color="#1e8449",
                label="Newton + búsqueda de línea")
    a2.set_xlabel("iteración k")
    a2.set_ylabel("|x_k − x*|")
    a2.set_title("El mismo x0 = 2, dos destinos")
    a2.legend(fontsize=9)

    fig.savefig(f"{SALIDA}/fig1a_arctan.png", bbox_inches="tight")
    plt.close(fig)
    print(f"\n  → {SALIDA}/fig1a_arctan.png")


# =============================================================================
# PARTE B — Sistema 2×2 y trayectorias en el plano
# =============================================================================
def F_kelley(v):
    """Sistema de prueba (Kelley, 2003):

        F1(x1,x2) = x1² + x2² − 2
        F2(x1,x2) = e^(x1−1) + x2³ − 2

    Raíz: x* = (1, 1). La segunda componente crece exponencialmente en x1,
    así que un paso de Newton demasiado largo hacia la derecha desborda.
    Además J es SINGULAR sobre la recta x2 = 0 (la segunda columna se anula),
    lo que usamos en la Parte D.
    """
    return np.array([v[0] ** 2 + v[1] ** 2 - 2.0,
                     np.exp(v[0] - 1.0) + v[1] ** 3 - 2.0])


def J_kelley(v):
    return np.array([[2.0 * v[0], 2.0 * v[1]],
                     [np.exp(v[0] - 1.0), 3.0 * v[1] ** 2]])


def parte_B():
    titulo("PARTE B · Sistema 2×2 desde x0 = (2.0, 0.5)")

    print("""
    F1 = x1² + x2² − 2                       raíz:  x* = (1, 1)
    F2 = e^(x1−1) + x2³ − 2

  Las tres estrategias arrancan del MISMO punto y usan el MISMO paso de
  Newton. Lo único que cambia es cuánto de ese paso se acepta.
""")

    x0 = np.array([2.0, 0.5])
    corridas = {}
    for etiqueta, kw in [("Newton puro", dict(globalizacion="ninguna")),
                         ("Búsqueda de línea (Armijo)", dict(globalizacion="linea")),
                         ("Región de confianza (dogleg)", dict(globalizacion="region",
                                                               delta0=1.0))]:
        corridas[etiqueta] = nk_lib.newton_krylov(
            F_kelley, x0, forzado=1e-12, tol=1e-11, max_iter=30,
            guardar_trayectoria=True, **kw)

    filas = []
    for etiqueta, h in corridas.items():
        filas.append([etiqueta,
                      "sí" if h.convergio else "NO",
                      h.n_newton,
                      f"{h.residuales[-1]:.2e}",
                      h.n_F,
                      h.motivo])
    tabla(filas, ["estrategia", "¿converge?", "iters", "‖F‖ final",
                  "evals. de F", "motivo de parada"])

    # --- Figura: trayectorias sobre las curvas de nivel de ‖F‖ ---
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))

    g1 = np.linspace(-0.5, 3.2, 320)
    g2 = np.linspace(-1.6, 2.4, 320)
    G1, G2 = np.meshgrid(g1, g2)
    Z = np.log10(np.sqrt((G1 ** 2 + G2 ** 2 - 2) ** 2 +
                         (np.exp(G1 - 1) + G2 ** 3 - 2) ** 2) + 1e-12)
    cf = a1.contourf(G1, G2, Z, levels=28, cmap="Blues_r", alpha=0.85)
    fig.colorbar(cf, ax=a1, label="log₁₀ ‖F(x)‖")
    # las dos curvas F1 = 0 y F2 = 0; se cortan en la raíz
    a1.contour(G1, G2, G1 ** 2 + G2 ** 2 - 2, levels=[0], colors="white",
               linewidths=1.4, linestyles="--")
    a1.contour(G1, G2, np.exp(G1 - 1) + G2 ** 3 - 2, levels=[0], colors="white",
               linewidths=1.4, linestyles=":")

    estilos = {"Newton puro": ("#c0392b", "o-"),
               "Búsqueda de línea (Armijo)": ("#1e8449", "s-"),
               "Región de confianza (dogleg)": ("#8e44ad", "^-")}
    for etiqueta, h in corridas.items():
        color, marca = estilos[etiqueta]
        T = np.array(h.trayectoria)
        T = T[np.all(np.isfinite(T), axis=1)]
        a1.plot(T[:, 0], T[:, 1], marca, color=color, ms=5, lw=1.6,
                label=etiqueta, alpha=0.95)
    a1.plot(1, 1, "*", color="gold", ms=20, mec="black", mew=0.8, label="raíz (1,1)")
    a1.plot(*x0, "X", color="black", ms=11, label="x0 = (2, 0.5)")
    a1.set_xlim(g1[0], g1[-1])
    a1.set_ylim(g2[0], g2[-1])
    a1.set_xlabel("x₁")
    a1.set_ylabel("x₂")
    a1.set_title("Trayectorias sobre las curvas de nivel de ‖F‖")
    a1.legend(fontsize=8, loc="upper left")

    for etiqueta, h in corridas.items():
        color, marca = estilos[etiqueta]
        r = np.array(h.residuales, dtype=float)
        r = np.where(np.isfinite(r), r, np.nan)
        a2.semilogy(r, marca, color=color, ms=5, lw=1.6, label=etiqueta)
    a2.set_xlabel("iteración de Newton k")
    a2.set_ylabel("‖F(x_k)‖")
    a2.set_title("Historia del residuo (mismo x0)")
    a2.legend(fontsize=8)

    fig.savefig(f"{SALIDA}/fig1b_trayectorias.png", bbox_inches="tight")
    plt.close(fig)
    print(f"\n  → {SALIDA}/fig1b_trayectorias.png")
    return corridas


# =============================================================================
# PARTE C — Mapa de cuencas de convergencia
# =============================================================================
def newton_2x2(x0, estrategia, max_iter=12, tol=1e-10):
    """Newton para n = 2 con el sistema lineal resuelto de forma exacta.

    Aquí NO usamos Krylov a propósito: en dimensión 2 resolver J s = −F es
    trivial, y lo que queremos aislar en esta parte es el efecto de la
    GLOBALIZACIÓN, no el del solver lineal. Devuelve (convergió, iteraciones).
    """
    x = np.array(x0, dtype=float)
    Fx = F_kelley(x)
    nF = np.linalg.norm(Fx)

    for k in range(max_iter):
        if nF <= tol:
            return True, k
        if not np.isfinite(nF) or nF > 1e12:
            return False, k
        J = J_kelley(x)
        if abs(np.linalg.det(J)) < 1e-13:
            return False, k
        s = np.linalg.solve(J, -Fx)

        if estrategia == "puro":
            x = x + s
            Fx = F_kelley(x)
            nF = np.linalg.norm(Fx)
        else:                                    # Armijo
            lam = 1.0
            for _ in range(40):
                x_t = x + lam * s
                F_t = F_kelley(x_t)
                n_t = np.linalg.norm(F_t)
                if n_t <= (1.0 - 1e-4 * lam) * nF:
                    break
                lam *= 0.5
            else:
                return False, k
            x, Fx, nF = x_t, F_t, n_t
    return nF <= tol, max_iter


def parte_C():
    titulo("PARTE C · Mapa de cuencas con presupuesto de 12 iteraciones")

    print("""
  Barremos 51×51 = 2601 puntos iniciales y pintamos cada uno según cuántas
  iteraciones necesitó. El presupuesto es de 12 iteraciones: así se mide lo
  que de verdad importa en la práctica, no "¿converge alguna vez?" sino
  "¿converge con los recursos que tengo?".

  Dos aclaraciones de método, para que el experimento sea limpio:

    · La ventana excluye la franja |x2| < 0.15. Ahí J es singular (Parte D) y
      lo que se mediría es el efecto de la singularidad, no el de globalizar.
    · En dimensión 2 y con problemas suaves, Newton puro YA es bastante
      robusto. La diferencia aquí es real pero moderada. El abismo entre
      globalizar y no globalizar aparece en los sistemas grandes que vienen
      de discretizar una EDP: eso se ve en el ejercicio 3.
""")

    g1 = np.linspace(-1.5, 3.0, 51)
    g2 = np.linspace(0.15, 3.0, 51)
    mapas = {}
    for estrategia in ("puro", "armijo"):
        M = np.zeros((len(g2), len(g1)))
        for i, b in enumerate(g2):
            for j, a in enumerate(g1):
                ok, it = newton_2x2((a, b), estrategia)
                M[i, j] = it if ok else np.nan
        mapas[estrategia] = M

    filas = []
    for estrategia, M in mapas.items():
        exito = np.count_nonzero(~np.isnan(M))
        filas.append(["Newton puro" if estrategia == "puro" else "Newton + Armijo",
                      f"{exito} / {M.size}",
                      f"{100*exito/M.size:.1f} %",
                      f"{np.nanmedian(M):.0f}",
                      f"{np.nanmax(M):.0f}"])
    tabla(filas, ["estrategia", "puntos que convergen", "porcentaje",
                  "iters. medianas", "iters. máximas"])

    fig, ejes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    for ax, (estrategia, M) in zip(ejes, mapas.items()):
        nombre = "Newton puro" if estrategia == "puro" else "Newton + Armijo"
        im = ax.imshow(M, origin="lower", cmap="viridis_r",
                       extent=[g1[0], g1[-1], g2[0], g2[-1]],
                       aspect="auto", vmin=3, vmax=12)
        ax.set_facecolor("#7b241c")          # rojo oscuro = agotó el presupuesto
        pct = 100 * np.count_nonzero(~np.isnan(M)) / M.size
        ax.set_title(f"{nombre} — converge desde {pct:.1f} %")
        ax.plot(1, 1, "*", color="gold", ms=16, mec="black")
        ax.set_xlabel("x₁ inicial")
        ax.set_ylabel("x₂ inicial")
        fig.colorbar(im, ax=ax, label="iteraciones usadas")
    fig.suptitle("Cuencas con presupuesto de 12 iteraciones "
                 "(rojo oscuro = no llegó)", fontsize=13, fontweight="bold")
    fig.savefig(f"{SALIDA}/fig1c_cuencas.png", bbox_inches="tight")
    plt.close(fig)
    print(f"\n  → {SALIDA}/fig1c_cuencas.png")


# =============================================================================
# PARTE D — La letra chica: globalizar no es magia
# =============================================================================
def parte_D():
    titulo("PARTE D · Dónde la globalización TAMBIÉN falla")

    print("""
  El Jacobiano de este sistema es

        J = [ 2x1    2x2  ]
            [ e^(x1−1)  3x2² ]

  y sobre la recta x2 = 0 la segunda columna se anula: J es SINGULAR.
  La búsqueda de línea solo garantiza que ‖F‖ baje en cada paso; si el
  descenso arrastra la iteración hacia esa recta, el paso de Newton se hace
  enorme y mal condicionado, λ se va a cero y el método se atasca sin haber
  encontrado ninguna raíz. La región de confianza aguanta más porque acota
  el tamaño del paso, pero tampoco es una garantía.
""")

    filas = []
    for x0 in [(2.0, 0.5), (4.0, 0.0), (3.0, -1.0), (2.5, -2.0)]:
        fila = [str(x0)]
        for g, kw in [("ninguna", {}), ("linea", {}), ("region", dict(delta0=1.0))]:
            h = nk_lib.newton_krylov(F_kelley, np.array(x0), globalizacion=g,
                                     forzado=1e-12, tol=1e-11, max_iter=40, **kw)
            if h.convergio:
                fila.append(f"OK ({h.n_newton} it)")
            else:
                fila.append("falla")
        filas.append(fila)
    tabla(filas, ["x0", "Newton puro", "Armijo", "Región (dogleg)"])

    print("""
  Conclusión honesta para la exposición: la globalización garantiza que la
  función de mérito f(x) = ½‖F(x)‖² NO AUMENTE, y con eso convergencia a un
  punto estacionario de f. Un punto estacionario de f no tiene por qué ser
  una raíz de F: puede ser un mínimo local de f, o un punto donde J es
  singular. Eso es exactamente lo que dice el teorema, ni más ni menos.
""")


if __name__ == "__main__":
    import os
    os.makedirs(SALIDA, exist_ok=True)
    print(__doc__)
    parte_A()
    parte_B()
    parte_C()
    parte_D()
    titulo("FIN DEL EJERCICIO 1")
