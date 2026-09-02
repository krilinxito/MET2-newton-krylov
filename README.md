# Estrategias para la convergencia global — Métodos de Newton-Krylov

Paquete completo de exposición para **Métodos Numéricos II (DAT-252)**
UMSA · Facultad de Ciencias Puras y Naturales · Carrera de Informática
Docente: M.Sc. Carlos Mullisaca Choque
Expositores: Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

---

## Qué hay acá

| Entregable | Archivo | Cómo se usa |
|---|---|---|
| **Presentación** | `presentacion.html` | Doble clic. Un solo archivo, funciona sin internet. `→` avanza, `f` pantalla completa, `o` vista general, `?` ayuda. |
| **Informe** | `Informe_NewtonKrylov_DAT252.docx` | 21 páginas, formato de la plantilla institucional. |
| **Manual y guía** | `guia_de_estudio.md` | Dos partes: un **manual desde cero** que explica el tema sin dar nada por sabido (§1–§11), y el **guion de exposición** con minutaje, qué decir por diapositiva y preguntas probables con respuesta (§12–§18). |
| **Demostraciones en vivo** | `ejercicios_exposicion/` | Tres programas que se corren durante la exposición. |
| **Ejercicios para la clase** | `ejercicios_clase/` | Tres programas autónomos y comentados para repartir a los compañeros. |

La plantilla original `Gnombres_Dat252.docx` no se modifica: se usa solo como
molde de formato.

---

## Requisitos

```bash
pip install numpy scipy matplotlib
```

Python 3.9 o superior. Nada más: no hace falta `python-docx`, ni LaTeX, ni
conexión a internet.

---

## Correr todo

```bash
# 1. Demostraciones de la exposición (generan las figuras)
cd ejercicios_exposicion
python3 ej1_newton_vs_globalizado.py        # ~20 s  · por qué hace falta globalizar
python3 ej2_bratu1d_newton_krylov.py        # ~20 s  · el método completo sobre Bratu
python3 ej3_comparativa_globalizacion.py    # ~40 s  · las cuatro estrategias
cd ..

# 2. Ejercicios para los compañeros
cd ejercicios_clase
python3 clase1_armijo.py                    # ~5 s
python3 clase2_matrix_free.py               # ~25 s
python3 clase3_bratu_forcing.py             # ~35 s
cd ..

# 3. Regenerar la presentación y el informe con las figuras nuevas
python3 presentacion/build_presentacion.py
python3 informe/generar_informe.py
```

El orden importa: la presentación y el informe empotran las figuras que
producen los ejercicios.

---

## El tema, en una página

Se quiere resolver `F(x) = 0` con `F: Rⁿ → Rⁿ` y `n` grande. El método de
Newton tiene tres problemas prácticos:

| Problema | Solución | Dónde se ve |
|---|---|---|
| Converge solo **cerca** de la raíz | globalizar | ejercicios 1 y 3 |
| Formar y factorizar `J` es imposible | Krylov matriz-libre | ejercicio 2 |
| Resolver `J s = −F` exacto es un desperdicio | Newton inexacto | ejercicios 2 y 3 |

**La desigualdad que sostiene todo:** si `J s = −F + r` con `‖r‖ ≤ η‖F‖` y
`η < 1`, entonces

```
∇f(x)ᵀ s = Fᵀ J s = −‖F‖² + Fᵀr ≤ −(1 − η)‖F‖² < 0     con  f = ½‖F‖²
```

El paso de Newton inexacto siempre es dirección de descenso para la función de
mérito. Por eso la búsqueda de línea siempre encuentra un λ que sirve, y por eso
se exige η < 1.

---

## Resultados principales

Todos medidos con el código de este repositorio.

**Ecuación de Bratu 1D** (N = 250, λ = 3) — el término de forzado decide el costo:

| η | Iteraciones de Newton | Productos J·v |
|---|---|---|
| 10⁻¹ fijo | 7 | 9 854 |
| 10⁻¹² fijo | **4** | **24 804** |
| Eisenstat-Walker 2 | 8 | **9 707** |

El que menos iteraciones usa es el que más trabajo consume. Con precondicionador
laplaciano, Eisenstat-Walker baja de 9 707 a **22** productos J·v.

**Ecuación de Burgers 1D** (ν = 0.01, N = 200, 7 puntos iniciales) — la
globalización decide si se resuelve o no:

| Estrategia | Converge desde | J·v (mediana) |
|---|---|---|
| Sin globalizar | 1 de 7 | 10 631 |
| Búsqueda de línea (Armijo) | 5 de 7 | 3 114 |
| Región de confianza (dogleg) | 3 de 7 | 3 174 |
| Continuación pseudo-transitoria | **7 de 7** | **863** |

**Y el límite honesto:** sobre la función de Freudenstein-Roth, Newton *puro*
alcanza la raíz desde el 83 % de los puntos iniciales y Newton *con Armijo* solo
desde el 37 %, porque la búsqueda de línea queda atrapada en un mínimo local de
la función de mérito. La globalización garantiza que `‖F‖` no aumente, no que se
encuentre la raíz.

---

## Estructura

```
MET2/
├── presentacion.html                  ← deck autónomo (24 diapositivas)
├── Informe_NewtonKrylov_DAT252.docx   ← informe final
├── guia_de_estudio.md                 ← manual desde cero + guion de exposición
├── Gnombres_Dat252.docx               ← plantilla institucional (no se toca)
│
├── ejercicios_exposicion/
│   ├── nk_lib.py                      ← núcleo: J·v, forzado, Armijo, dogleg, Steihaug, Ψtc
│   ├── ej1_newton_vs_globalizado.py
│   ├── ej2_bratu1d_newton_krylov.py
│   ├── ej3_comparativa_globalizacion.py
│   └── figuras/
│
├── ejercicios_clase/
│   ├── README.md                      ← enunciados y preguntas de análisis
│   ├── clase1_armijo.py
│   ├── clase2_matrix_free.py
│   └── clase3_bratu_forcing.py
│
├── presentacion/
│   ├── plantilla.html                 ← fuente del deck (con marcadores %%FIG:...%%)
│   └── build_presentacion.py          ← empotra las figuras en base64
│
└── informe/
    ├── docx_min.py                    ← constructor OOXML mínimo (sin python-docx)
    ├── generar_informe.py             ← contenido del informe
    └── figuras/                       ← ecuaciones renderizadas y logo
```

---

## Referencias

- Kelley, C. T. (2003). *Solving Nonlinear Equations with Newton's Method*. SIAM.
- Dembo, Eisenstat & Steihaug (1982). «Inexact Newton Methods». *SIAM J. Numer. Anal.* 19(2), 400–408.
- Eisenstat & Walker (1996). «Choosing the forcing terms in an inexact Newton method». *SIAM J. Sci. Comput.* 17(1), 16–32.
- Knoll & Keyes (2004). «Jacobian-free Newton-Krylov methods». *J. Comput. Phys.* 193, 357–397.
- Pawlowski, Shadid, Simonis & Walker (2006). «Globalization techniques for Newton-Krylov methods». *SIAM Review* 48(4), 700–721.
