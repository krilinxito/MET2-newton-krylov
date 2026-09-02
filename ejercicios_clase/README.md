# Notebooks — Métodos de Newton-Krylov aplicados

**Métodos Numéricos II (DAT-252) — UMSA, Carrera de Informática**
Docente: M.Sc. Carlos Mullisaca Choque
Expositores: Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

Tres problemas reales resueltos con métodos de Newton-Krylov. **No hay nada que
resolver**: se abre el notebook, se corren todas las celdas y se mira qué pasa.

---

## Cómo correrlos

**Opción 1 — en su computadora**

```bash
pip install numpy scipy matplotlib jupyterlab
jupyter lab
```

**Opción 2 — en Google Colab**
Subir el `.ipynb` a [colab.research.google.com](https://colab.research.google.com).
Ahí numpy, scipy y matplotlib ya vienen instalados.

En cualquiera de las dos: **Ejecutar todas las celdas** (`Run All`), de arriba abajo.
Los tres notebooks son independientes entre sí y no necesitan ningún otro archivo.

---

## Los tres notebooks

### `01_circuito_con_diodos.ipynb` · 2 s

El punto de operación de un circuito con diodos. Las leyes de Kirchhoff son lineales,
pero la corriente del diodo va como `exp(V/Vt)` y eso vuelve no lineal al sistema.

Newton con el paso completo manda 5 voltios al diodo, `exp(5/0.0259)` desborda y el
programa se muere en la primera iteración. Recortando el paso hasta que el residuo
baje, el circuito se resuelve en 8 iteraciones. **Es exactamente lo que hace SPICE.**

### `02_placa_que_irradia.ipynb` · 3 s

Temperatura de una placa de acero de 20 × 20 cm con un componente caliente en el
centro, que disipa por conducción y por radiación. El término de Stefan-Boltzmann va
como `T⁴`.

Malla de 50 × 50 = **2 500 incógnitas**, así que el Jacobiano tendría 6 250 000
entradas. No se construye: el método solo lo multiplica por vectores, y cada producto
sale de una evaluación más de la función. Se ve también qué cambia al precondicionar.

### `03_ignicion_termica.ipynb` · 15 s

Una lámina de material que genera calor por una reacción cuya velocidad crece
exponencialmente con su propia temperatura. Existe una **potencia crítica**
`λ* ≈ 3.5138`: por encima de ella no hay estado estacionario, y eso es la ignición.

Se sube la potencia de a poco usando cada solución como punto de partida de la
siguiente, se dibuja la rama de soluciones hasta el límite, y se ve qué pasa al
pasarse.

---

## De qué va el tema

Resolver `F(x) = 0` con `F: Rⁿ → Rⁿ` y `n` grande. El método de Newton resuelve, en
cada paso, un sistema **lineal**:

```
    J(x_k) · s_k = −F(x_k)          y luego    x_{k+1} = x_k + s_k
```

y tiene tres problemas que estos notebooks van mostrando:

| Problema | Solución | Notebook |
|---|---|---|
| Converge solo **cerca** de la raíz | recortar el paso (globalizar) | **01** |
| Formar y guardar `J` es imposible si `n` es grande | `J·v ≈ [F(x+εv) − F(x)]/ε` | **02** |
| Resolver `J s = −F` exacto es un desperdicio | resolverlo solo hasta `η‖F‖` | **02**, **03** |

La desigualdad que hace que todo encaje: si el paso `s` cumple `J s = −F + r` con
`‖r‖ ≤ η‖F‖` y `η < 1`, entonces para `f(x) = ½‖F(x)‖²`

```
    ∇f(x)ᵀ s  =  −‖F‖² + Fᵀr  ≤  −(1 − η)‖F‖²  <  0
```

es decir, **el paso de Newton siempre apunta hacia donde el residuo baja**, aunque
esté mal calculado. Por eso siempre existe un paso corto que mejora las cosas.
