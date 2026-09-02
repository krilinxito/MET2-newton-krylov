# Ejercicios de clase — Métodos de Newton-Krylov y convergencia global

**Materia:** Métodos Numéricos II (DAT-252) — UMSA, Carrera de Informática
**Docente:** M.Sc. Carlos Mullisaca Choque
**Expositores:** Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

---

## Antes de empezar

```bash
pip install numpy scipy matplotlib
```

Los tres programas son **autónomos**: cada uno se puede copiar suelto y correr
sin nada más. No importan módulos del proyecto ni se necesita internet.

```bash
python3 clase1_armijo.py          # ~5 s
python3 clase2_matrix_free.py     # ~25 s
python3 clase3_bratu_forcing.py   # ~35 s
```

Cada uno imprime tablas por consola, guarda sus figuras en la misma carpeta y
termina con un bloque **PREGUNTAS PARA ANALIZAR**. Las preguntas son el
ejercicio; el código está completo justamente para que el tiempo se gaste en
entender, no en depurar.

---

## El hilo conductor

Queremos resolver `F(x) = 0` donde `F: Rⁿ → Rⁿ` y `n` es grande (viene de
discretizar una EDP). El método de Newton dice:

```
resolver   J(x_k) s_k = −F(x_k)
actualizar x_{k+1} = x_k + s_k
```

y tiene tres problemas prácticos, uno por ejercicio:

| Problema | Solución | Ejercicio |
|---|---|---|
| Newton solo converge **cerca** de la raíz | globalizar: búsqueda de línea, región de confianza, Ψtc | **1** |
| Formar y factorizar `J` es imposible si `n` es grande | Krylov matriz-libre: `J·v ≈ [F(x+εv) − F(x)]/ε` | **2** |
| Resolver `J s = −F` exactamente es un desperdicio | Newton inexacto: `‖Js + F‖ ≤ η‖F‖` con η adaptativo | **2 y 3** |

Los tres juntos dan el método **Newton-Krylov globalizado**, que es lo que usan
PETSc (SNES), SUNDIALS (KINSOL) y `scipy.optimize.newton_krylov`.

---

## Ejercicio 1 — `clase1_armijo.py`

### Búsqueda de línea de Armijo: qué garantiza y qué no

Resuelve dos sistemas 2×2 con Newton puro y con Newton + Armijo, y dibuja las
cuencas de convergencia de ambos.

- **Problema A** — circunferencia contra exponencial empinada. La búsqueda de
  línea convierte desbordamientos en raíces encontradas.
- **Problema B** — función de Freudenstein y Roth. Aquí pasa lo contrario:
  Newton puro encuentra la raíz desde muchos más puntos que Newton con Armijo.
  **No es un error del programa.** Es la lección del ejercicio.

**Qué hay que entender:** la búsqueda de línea garantiza que `‖F‖` no aumente
nunca. No garantiza llegar a una raíz. Si la función de mérito `f = ½‖F‖²`
tiene un mínimo local que no es raíz —y la de Freudenstein-Roth lo tiene, en
`x₂ ≈ −0.8968` con `f ≈ 49`— una sucesión que solo sabe bajar cae ahí y se
queda.

**Salidas:** `clase1_cuencas.png`

---

## Ejercicio 2 — `clase2_matrix_free.py`

### El Jacobiano que nunca se construye

Problema modelo: `−ν u'' + u³ = f(x)`, con `u(0) = u(1) = 0` y solución exacta
conocida `u(x) = sin(πx)`, para poder medir el error de verdad y no solo el
residuo.

Tres experimentos:

1. **El paso ε.** Barrido de `ε` en escala logarítmica. Se ve el mínimo donde
   se cruzan el error de truncamiento (`~ε`) y el de cancelación
   (`~eps_maq/ε`). También se ve que el óptimo *medido* no cae exactamente en
   `√eps_maq`, y por qué eso no importa en la práctica.
2. **El costo de formar J.** Cuántas evaluaciones de `F` cuesta armar la
   matriz contra cuántas cuesta un producto `J·v`.
3. **El oversolving.** Se resuelve el mismo problema con η fijo desde 0.5 hasta
   1e-12 y se cuenta el trabajo real. El η más exigente usa *menos* iteraciones
   de Newton y *mucho más* trabajo total, sin ganar ni una cifra de precisión
   en la respuesta final.

**Salidas:** `clase2_epsilon.png`, `clase2_oversolving.png`

---

## Ejercicio 3 — `clase3_bratu_forcing.py`

### Newton-Krylov completo sobre la ecuación de Bratu

`u'' + λ e^u = 0`, `u(0) = u(1) = 0`. Modelo clásico de ignición térmica: tiene
solución solo para `λ ≤ λ* ≈ 3.5138`, y en `λ*` el Jacobiano se vuelve singular.

Tres experimentos:

1. **Cinco maneras de elegir η**, con y sin precondicionador. La ganancia del
   precondicionador laplaciano es de dos órdenes de magnitud, salvo para
   `η = 1e-8`, que no mejora: hay que entender por qué.
2. **La trayectoria de η** que elige Eisenstat-Walker, iteración por iteración.
3. **Acercándose a λ\***: cómo crece `cond(J)` y qué se ve cuando el problema
   simplemente deja de tener solución.

**Salidas:** `clase3_forzado.png`, `clase3_punto_retorno.png`

---

## Una sola tabla para llevarse

| Concepto | Fórmula | Qué controla |
|---|---|---|
| Paso de Newton | `J s = −F` | la dirección |
| Newton inexacto | `‖J s + F‖ ≤ η‖F‖` | cuánto trabajo se gasta en el sistema lineal |
| Producto matriz-libre | `J·v ≈ [F(x+εv) − F(x)]/ε`, `ε ≈ √eps_maq` | no formar nunca `J` |
| Función de mérito | `f(x) = ½‖F(x)‖²` | qué significa "ir mejorando" |
| Condición de Armijo | `‖F(x+λs)‖ ≤ (1 − αλ)‖F(x)‖`, `α = 1e-4` | cuánto del paso se acepta |
| Eisenstat-Walker 2 | `η_k = 0.9 (‖F_k‖/‖F_{k−1}‖)^1.618` | η automático |

**La desigualdad que justifica todo:** si `J s = −F + r` con `‖r‖ ≤ η‖F‖` y
`η < 1`, entonces

```
∇f(x)ᵀ s  =  Fᵀ J s  =  −‖F‖² + Fᵀr  ≤  −(1 − η)‖F‖²  <  0
```

es decir, **el paso de Newton inexacto es siempre dirección de descenso para la
función de mérito, siempre que η < 1.** Por eso la búsqueda de línea siempre
encuentra algún λ que funciona, y por eso se exige η < 1 y no cualquier cosa.

---

## Referencias

- C. T. Kelley, *Solving Nonlinear Equations with Newton's Method*, SIAM, 2003.
- Dembo, Eisenstat & Steihaug, "Inexact Newton Methods", *SIAM J. Numer. Anal.*
  19(2), 400–408, 1982.
- Eisenstat & Walker, "Choosing the forcing terms in an inexact Newton method",
  *SIAM J. Sci. Comput.* 17(1), 16–32, 1996.
- Knoll & Keyes, "Jacobian-free Newton-Krylov methods: a survey of approaches
  and applications", *J. Comput. Phys.* 193, 357–397, 2004.
