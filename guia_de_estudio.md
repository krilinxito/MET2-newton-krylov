# Manual y guía de exposición
## Estrategias para la convergencia global — Métodos de Newton-Krylov

**Materia:** Métodos Numéricos II (DAT-252) — UMSA, Carrera de Informática
**Docente:** M.Sc. Carlos Mullisaca Choque
**Expositores:** Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

---

## 0. Cómo usar este manual

### 0.1. Dos maneras de leer esto

Este archivo tiene dos partes que sirven para cosas distintas.

**PARTE I (secciones 1 a 11) — el manual.** Explica el tema desde cero. No supone
que usted sepa qué es un Jacobiano, ni qué es una norma, ni qué es un método de
Krylov. Cada cosa se define antes de usarla, y casi todo viene con un ejemplo
numérico hecho a mano.

**PARTE II (secciones 12 a 18) — el guion.** Es la chuleta para el día de la
exposición: cuánto dura cada bloque, qué decir en cada diapositiva, qué escribir en
la pizarra y qué responder si preguntan.

Según su situación:

| Si usted… | Lea |
|---|---|
| Expone el viernes y ya entiende el tema | PARTE II, y §0.3 para refrescar |
| Expone el viernes y **no** entiende bien el tema | §0.3 → §3 → §6 → PARTE II |
| Quiere entenderlo de verdad | Todo, en orden, corriendo el código de §8 |
| Solo tiene 10 minutos | §0.3 |

A lo largo de la PARTE II hay marcas como **(→ §6.5)** que remiten a la sección de
la PARTE I donde eso se explica en serio.

### 0.2. Qué se da por sabido y qué no

**Se da por sabido** (cosas de Cálculo I–II y Álgebra Lineal básica):

- Qué es una derivada y qué significa geométricamente (la pendiente de la tangente).
- Qué es una derivada parcial `∂f/∂x`.
- Qué es un vector y qué es una matriz, y cómo se multiplica una matriz por un
  vector.
- Qué significa resolver un sistema lineal `A x = b`.

**No se da por sabido** — todo esto se explica acá:

- Qué es un sistema de ecuaciones **no** lineales y por qué es mucho más difícil que
  uno lineal (§1).
- Qué es una **norma** y por qué se usa (§2.1).
- Qué es el **Jacobiano** y cómo se arma (§2.2).
- Qué es el **gradiente** de una función escalar y qué significa que un vector
  «apunte cuesta abajo» (§2.3).
- Qué es el **épsilon de máquina** y por qué las computadoras no calculan exacto
  (§2.4).
- Qué es el **número de condición** (§2.5).
- Qué diferencia hay entre un método **directo** y uno **iterativo** (§2.6).
- Qué significa **orden de convergencia** lineal, superlineal o cuadrático (§2.7).
- Qué es un **subespacio** y qué es `span` (§5.1).

Si algo de la primera lista tampoco le suena, cualquier libro de Cálculo sirve; nada
de lo que sigue depende de manipulaciones algebraicas difíciles.

### 0.3. El tema entero, en una página

Queremos resolver un sistema de ecuaciones no lineales:

```
    F(x) = 0        donde x tiene n componentes y F devuelve n números
```

El **método de Newton** es la herramienta clásica. Funciona así: en cada paso
reemplaza `F` por su aproximación lineal, resuelve esa aproximación, y se mueve ahí.

```
    resolver   J(x_k) · s_k = −F(x_k)          ← un sistema LINEAL
    avanzar    x_{k+1} = x_k + s_k
```

Cuando funciona, es rapidísimo: **duplica las cifras correctas en cada paso**.
Pero tiene tres problemas graves en cuanto el problema es grande o el punto de
partida no es bueno:

```
┌─ PROBLEMA 1 ─────────────────────────────┐   ┌─ SOLUCIÓN ──────────────────┐
│ Resolver J s = −F exacto cuesta O(n³).   │ → │ NEWTON INEXACTO             │
│ Con n = 10⁶ eso es inviable.             │   │ Resolverlo solo "más o      │
│                                          │   │ menos": ‖Js + F‖ ≤ η‖F‖     │
└──────────────────────────────────────────┘   └─────────────────────────────┘

┌─ PROBLEMA 2 ─────────────────────────────┐   ┌─ SOLUCIÓN ──────────────────┐
│ Ni siquiera se puede FORMAR la matriz J. │ → │ KRYLOV MATRIZ-LIBRE         │
│ Con n = 10⁶ son 10¹² números: 8 TB.      │   │ J·v ≈ [F(x+εv) − F(x)]/ε    │
│                                          │   │ La matriz nunca existe.     │
└──────────────────────────────────────────┘   └─────────────────────────────┘

┌─ PROBLEMA 3 ─────────────────────────────┐   ┌─ SOLUCIÓN ──────────────────┐
│ Newton solo converge CERCA de la raíz.   │ → │ GLOBALIZACIÓN               │
│ Desde lejos puede diverger.              │   │ Exigir que ‖F‖ nunca suba.  │
│                                          │   │ ← EL TEMA DE HOY            │
└──────────────────────────────────────────┘   └─────────────────────────────┘
```

Y las tres soluciones encajan gracias a **una sola desigualdad**. Si el paso `s`
resuelve el sistema lineal solo aproximadamente —es decir `J s = −F + r` con un
error `r` que cumple `‖r‖ ≤ η‖F‖`— entonces, midiendo el progreso con la función
`f(x) = ½‖F(x)‖²`:

```
    ∇f(x)ᵀ s  =  −‖F‖² + Fᵀr  ≤  −(1 − η)‖F‖²  <  0     siempre que  η < 1
```

En castellano: **el paso de Newton, aunque esté mal calculado, siempre apunta hacia
donde `‖F‖` baja, con tal de que el error del sistema lineal sea menor que el propio
residuo.** Por eso siempre existe un paso corto que mejora las cosas, y por eso la
globalización funciona.

Todo lo demás de este manual es desarrollar esas cuatro cajas.

---

# PARTE I — El tema desde cero

## 1. El problema

### 1.1. De una ecuación a un sistema

Usted ya resolvió ecuaciones como `x² − 2 = 0`. Hay **una** incógnita y **una**
ecuación, y la respuesta es un número: `x = √2`.

Ahora imagine que hay **dos** incógnitas y **dos** ecuaciones, y que las ecuaciones
no son rectas:

```
    x₁² + x₂² − 2       = 0
    e^(x₁−1) + x₂³ − 2  = 0
```

La respuesta ya no es un número sino un **par** de números. En este caso la
respuesta es `x₁ = 1`, `x₂ = 1` (compruébelo: `1 + 1 − 2 = 0` y `e⁰ + 1 − 2 = 0`).

La palabra clave es **no lineal**. Un sistema es *lineal* si cada incógnita aparece
sola, multiplicada por un número, y sumada — como `3x₁ + 5x₂ = 7`. En cuanto aparece
un `x²`, un `e^x`, un `sen x` o un producto `x₁·x₂`, el sistema es **no lineal**.

La diferencia es enorme y conviene tenerla clara desde el principio:

| | Sistema lineal | Sistema no lineal |
|---|---|---|
| Forma | `A x = b` | `F(x) = 0` |
| ¿Cuántas soluciones? | 0, 1 o infinitas — se sabe de antemano | puede haber 0, 1, 7 o infinitas, y **no se sabe** |
| ¿Cómo se resuelve? | eliminación de Gauss, en un número fijo de pasos | **iterando**: se adivina y se corrige, sin garantía de llegar |
| ¿Cuánto tarda? | se puede calcular exactamente | depende del punto de partida |

Ese «depende del punto de partida» es el origen de todo este tema.

### 1.2. La notación, desarmada

Va a ver escrito esto:

```
    F : Rⁿ → Rⁿ ,      F(x) = 0
```

Símbolo por símbolo:

- **`Rⁿ`** es el conjunto de las listas de `n` números reales. `R²` son los pares
  como `(1.5, −3)`; `R¹⁰⁰⁰` son las listas de mil números.
- **`x`** (en negrita, o simplemente entendido como vector) es la incógnita: **no es
  un número, es una lista de `n` números** `x = (x₁, x₂, …, xₙ)`.
- **`F`** es una función que **toma** una lista de `n` números y **devuelve** otra
  lista de `n` números. La flecha `Rⁿ → Rⁿ` dice exactamente eso: entran `n`, salen
  `n`.
- **`F(x) = 0`** significa que las `n` componentes de la salida valen cero **todas a
  la vez**. No es una ecuación: son `n` ecuaciones simultáneas.

Que entren `n` y salgan `n` no es casualidad: si hay más ecuaciones que incógnitas
el sistema normalmente no tiene solución, y si hay menos, tiene infinitas. El caso
«cuadrado» es el que tiene una solución aislada, que es lo que buscamos.

Al valor `F(x)` se le llama **residuo**: mide cuánto le falta a `x` para ser
solución. Si `F(x)` da cero, `x` es la respuesta. Si da algo grande, `x` está lejos.

### 1.3. Un ejemplo con números

Tomemos el sistema de arriba y probemos con `x = (2, 0.5)`:

```
    F₁(2, 0.5) = 2² + 0.5² − 2       = 4 + 0.25 − 2      = +2.25
    F₂(2, 0.5) = e^(2−1) + 0.5³ − 2  = 2.7183 + 0.125 − 2 = +0.8433
```

Así que `F(2, 0.5) = (2.25, 0.8433)`. No es cero, o sea que `(2, 0.5)` **no** es
solución. ¿Cuán lejos está? Para responder eso hace falta una sola manera de medir
un vector, y eso es la norma (→ §2.1).

Este sistema concreto va a reaparecer muchas veces en el manual, así que vale la
pena bautizarlo: le llamaremos **el sistema 2×2**. Su solución exacta es `(1, 1)`.

### 1.4. De dónde salen estos sistemas en la vida real

Casi nadie se encuentra un sistema no lineal escrito así. Se lo encuentra al final
de otra cosa. El caso más común es **discretizar una ecuación en derivadas
parciales**.

Un ejemplo real, la **ecuación de Bratu**, que modela la ignición de un material
combustible. Dice que la temperatura `u(x)` a lo largo de una barra de longitud 1
cumple:

```
    u''(x) + λ·e^(u(x)) = 0 ,      u(0) = 0 ,  u(1) = 0
```

El término `u''` es la difusión del calor y `λ·e^u` es el calor que genera la
reacción química (crece exponencialmente con la temperatura: por eso hay ignición).
Esto es una ecuación diferencial: la incógnita no es un número, es una **función**
`u(x)`, y hay infinitos valores de `x`.

**Discretizar** significa dejar de buscar la función entera y buscar solo su valor
en unos cuantos puntos igualmente espaciados:

```
   x=0                                                            x=1
    |----o-----o-----o-----o-----o-----o-----o-----o-----o-----|
    u=0  u₁    u₂    u₃    u₄    u₅    u₆    u₇    u₈    u₉   u=0
         └──────────── n incógnitas ────────────────┘
              separación entre puntos:  h = 1/(n+1)
```

La derivada segunda en el punto `i` se aproxima con los vecinos (esto se ve en
Métodos Numéricos I y sale de un desarrollo de Taylor):

```
    u''(xᵢ)  ≈  (u_{i−1} − 2uᵢ + u_{i+1}) / h²
```

Sustituyendo en la ecuación, para cada punto `i = 1, …, n` queda:

```
    (u_{i−1} − 2uᵢ + u_{i+1}) / h²  +  λ·e^(uᵢ)  =  0
```

Y eso es exactamente un sistema de `n` ecuaciones **no lineales** (por el `e^u`) con
`n` incógnitas `u₁, …, uₙ`. Es decir: **`F(u) = 0`**.

La misma receta —discretizar y obtener un sistema no lineal— aparece en:

- Flujo de fluidos (Navier-Stokes estacionario).
- Deformación de estructuras con material no lineal.
- El paso implícito de un integrador temporal para problemas rígidos.
- Condiciones de optimalidad de un problema de optimización.
- Equilibrio de un circuito eléctrico con diodos o transistores.

### 1.5. Por qué `n` se vuelve enorme

En el dibujo de arriba puse 9 puntos. Con 9 puntos la solución sale grosera. Para
que sirva hacen falta cientos, y en dos o tres dimensiones el número explota:

| Problema | Puntos por lado | `n` (incógnitas) | Entradas de la matriz `n²` | Memoria si se guarda densa |
|---|---|---|---|---|
| Barra 1D | 1 000 | 10³ | 10⁶ | 8 MB |
| Placa 2D | 1 000 | 10⁶ | 10¹² | **8 TB** |
| Cubo 3D | 100 | 10⁶ | 10¹² | **8 TB** |
| Cubo 3D | 1 000 | 10⁹ | 10¹⁸ | inimaginable |

Ocho terabytes es más memoria de la que tiene cualquier computadora que usted vaya a
usar, y eso es solo para *guardar* la matriz: invertirla sería mucho peor.

**Este es el hecho que gobierna todo el tema.** Cualquier método que necesite
escribir la matriz completa está descartado de entrada. Ténganlo presente: no es que
los métodos de Newton-Krylov sean «más elegantes», es que son los únicos posibles.

---

## 2. Las herramientas

Esta sección es un cajón de herramientas. Cada una viene con un «para qué la vamos a
usar», así que si alguna ya la conoce, salte y vuelva cuando aparezca.

### 2.1. La norma de un vector

**El problema que resuelve:** `F(x)` es una lista de `n` números. Para decir «este
punto está más cerca de la solución que aquel» hace falta convertir esa lista en
**un solo número** que mida su tamaño.

**Qué es.** La norma euclídea (o norma 2) de un vector es la raíz de la suma de los
cuadrados. Es la longitud de la flecha, el teorema de Pitágoras en `n` dimensiones:

```
    ‖v‖ = √( v₁² + v₂² + … + vₙ² )
```

Se escribe con **dobles barras** `‖v‖` para distinguirla del valor absoluto `|a|` de
un número suelto.

**Ejemplo.** Para el residuo que calculamos en §1.3:

```
    ‖F(2, 0.5)‖ = √( 2.25² + 0.8433² ) = √( 5.0625 + 0.7112 ) = √5.7737 = 2.4028
```

**Propiedades que vamos a usar** (las tres son intuitivas si piensa en longitudes):

1. `‖v‖ ≥ 0`, y vale 0 **solo** si `v` es el vector de puros ceros.
2. `‖a·v‖ = |a|·‖v‖` — estirar el vector estira su longitud.
3. **Desigualdad de Cauchy-Schwarz:** `|uᵀv| ≤ ‖u‖·‖v‖`.
   (Aquí `uᵀv = u₁v₁ + u₂v₂ + … + uₙvₙ` es el producto escalar.) Dice que el producto
   escalar nunca es más grande que el producto de las longitudes. Esta la vamos a
   necesitar en §6.5, que es la sección más importante del manual.

**Para qué la usamos.** «Convergimos» significa `‖F(x)‖ ≤ tolerancia`. Todos los
criterios de parada y todas las gráficas de este trabajo miden `‖F‖`.

> ⚠ **Cuidado con una trampa.** Que `‖F(x)‖` sea chico **no** garantiza que `x` esté
> cerca de la solución. Lo que acota el error es `‖J⁻¹‖·‖F‖`. En §9 hay un caso real
> de este proyecto donde `‖F‖ = 6·10⁻¹⁰` y aun así la solución difiere de la exacta
> en 0.3.

### 2.2. El Jacobiano

**El problema que resuelve:** para aproximar `F` por algo lineal hace falta su
derivada. Pero `F` tiene `n` entradas y `n` salidas: ¿qué es «la derivada» de eso?

**Qué es.** El **Jacobiano** es la matriz de todas las derivadas parciales. La fila
`i` corresponde a la ecuación `Fᵢ`, y la columna `j` a la incógnita `xⱼ`:

```
            ⎡ ∂F₁/∂x₁   ∂F₁/∂x₂   ⋯   ∂F₁/∂xₙ ⎤
    J(x) =  ⎢ ∂F₂/∂x₁   ∂F₂/∂x₂   ⋯   ∂F₂/∂xₙ ⎥
            ⎢    ⋮          ⋮      ⋱      ⋮    ⎥
            ⎣ ∂Fₙ/∂x₁   ∂Fₙ/∂x₂   ⋯   ∂Fₙ/∂xₙ ⎦
```

En la regla mnemotécnica de siempre: **fila = ecuación, columna = incógnita.**

**Ejemplo, hecho a mano.** Para nuestro sistema 2×2:

```
    F₁(x₁,x₂) = x₁² + x₂² − 2            F₂(x₁,x₂) = e^(x₁−1) + x₂³ − 2

    ∂F₁/∂x₁ = 2x₁      ∂F₁/∂x₂ = 2x₂
    ∂F₂/∂x₁ = e^(x₁−1)  ∂F₂/∂x₂ = 3x₂²

                  ⎡  2x₁      2x₂  ⎤
    entonces  J = ⎢                ⎥
                  ⎣ e^(x₁−1)  3x₂² ⎦
```

Evaluado en `x = (2, 0.5)`:

```
         ⎡ 4.0000   1.0000 ⎤
    J =  ⎢                 ⎥        determinante = 4·0.75 − 1·2.7183 = 0.2817
         ⎣ 2.7183   0.7500 ⎦
```

Ese determinante tan chico (0.28 frente a entradas de tamaño 4) es una señal de que
la matriz está **casi** en el borde de ser no invertible. Guárdelo: en §3.4 vamos a
ver que eso es exactamente lo que hace estallar al método de Newton desde ese punto.

**Para qué lo usamos.** El Jacobiano es lo que convierte el problema no lineal en una
sucesión de problemas lineales. Y también es el objeto que **no podemos permitirnos
construir** cuando `n` es grande (§1.5), lo que da origen al Problema 2.

> **Analogía.** En una dimensión, la derivada `f'(x)` dice «si muevo `x` un poquito,
> `f` se mueve `f'(x)` veces ese poquito». El Jacobiano dice lo mismo pero para
> muchas variables a la vez: si muevo el vector `x` en la dirección `v`, la salida
> `F` se mueve aproximadamente `J·v`. Esa frase es literalmente lo que vamos a
> explotar en §5.4.

### 2.3. El gradiente y qué significa «cuesta abajo»

**El problema que resuelve:** queremos una noción de «voy mejorando». Para eso
haremos que un solo número mida el progreso, y necesitamos saber en qué dirección
ese número baja.

**Qué es.** Si `f` es una función que toma un vector y devuelve **un número**
(no `n` números — **uno**), su **gradiente** es el vector de sus derivadas parciales:

```
    ∇f(x) = ( ∂f/∂x₁ ,  ∂f/∂x₂ , … , ∂f/∂xₙ )
```

El símbolo `∇` se lee «nabla» o simplemente «gradiente de».

**La propiedad clave.** El gradiente apunta en la dirección en la que `f` **más
crece**. Por lo tanto `−∇f` apunta hacia donde `f` **más baja**.

**Direcciones de descenso.** Supongamos que estamos en `x` y queremos movernos en la
dirección `s`. ¿Sube o baja `f`? La respuesta está en el producto escalar:

```
    ∇f(x)ᵀ s  <  0     ⟹   moverse un poco en la dirección s HACE BAJAR a f
    ∇f(x)ᵀ s  >  0     ⟹   moverse en la dirección s hace SUBIR a f
    ∇f(x)ᵀ s  =  0     ⟹   s es perpendicular al gradiente: f no cambia (de momento)
```

Cuando `∇f ᵀs < 0` se dice que **`s` es una dirección de descenso**.

```
        curvas de nivel de f            ∇f apunta hacia arriba de la colina
              ___________
            /   ______   \              s ┐
           |   /  ___  \  |               └──→  si el ángulo entre s y ∇f
           |  |  ( ⊙ ) |  |    ← mínimo         es MAYOR que 90°, entonces
           |   \  ‾‾‾  /  |                     ∇fᵀs < 0  y  s baja.
            \   ‾‾‾‾‾‾   /
              ‾‾‾‾‾‾‾‾‾‾‾
```

Ojo con el «un poco»: la garantía es **local**. `∇fᵀs < 0` asegura que existe algún
paso pequeño en la dirección `s` que baja, pero no dice cuán grande puede ser ese
paso. Buscar ese tamaño es precisamente lo que hace la búsqueda de línea (→ §6.6).

**Para qué lo usamos.** Toda la sección 6 —el tema de la exposición— se apoya en
demostrar que el paso de Newton es una dirección de descenso.

### 2.4. Punto flotante y el épsilon de máquina

**El problema que resuelve:** la computadora no guarda números reales; guarda
aproximaciones. Eso tiene consecuencias que en este tema son centrales, no
anecdóticas.

**Qué es.** Un número `double` usa 64 bits: unos para el signo, unos para el
exponente y 52 para las cifras significativas. Eso da unas **16 cifras decimales**.
El **épsilon de máquina** `εₘ` es el número más chico tal que `1 + εₘ ≠ 1` en la
computadora:

```
    εₘ ≈ 2.22 × 10⁻¹⁶
```

**La consecuencia peligrosa: la cancelación.** Cuando se restan dos números **muy
parecidos**, las cifras que coinciden se anulan y solo sobreviven las últimas, que
son justo las menos confiables:

```
    a = 1.2345678901234567
    b = 1.2345678901234512
    ---------------------------
    a − b = 0.0000000000000055     ← quedaron 2 cifras confiables de las 16
```

Si además ese resultado se divide entre algo diminuto, el error relativo se
multiplica. Esto va a ser exactamente lo que pase en §5.5.

**Para qué lo usamos.** Aparece dos veces y en las dos es decisivo:

1. En §5.5, para elegir el paso `ε` de la diferencia finita.
2. Como **piso de precisión**: si el Jacobiano se calcula con diferencias finitas, se
   conoce con unas 8 cifras, no 16. Por eso en este trabajo nunca se pide un residuo
   final mucho menor que `‖F(x₀)‖ · √εₘ`. No es pereza: es lo que hay.

### 2.5. El número de condición

**El problema que resuelve:** hay matrices con las que trabajar es fácil y otras con
las que todo sale mal. Necesitamos una medida de eso.

**Qué es.** El **número de condición** `cond(A)` mide cuánto amplifica la matriz los
errores. Si usted resuelve `A x = b` y `b` tiene un error relativo del 0.001 %, la
solución `x` puede tener un error relativo de hasta `cond(A) × 0.001 %`.

```
    cond(A) = 1        matriz perfecta (una rotación, por ejemplo)
    cond(A) = 10³      se pierden ~3 cifras de las 16: aceptable
    cond(A) = 10⁸      se pierden ~8 cifras: la mitad de la precisión
    cond(A) = 10¹⁶     se pierde todo: la respuesta es basura
```

**Dato que importa acá.** La matriz que sale de discretizar una derivada segunda
tiene `cond ≈ 1/h² = O(n²)`. O sea: **cuanto más fina es la malla, peor
condicionado está el problema**. Con `n = 250` eso ya son unas 25 000 unidades. Es la
razón de ser del precondicionamiento (→ §5.6).

**Un detalle que casi nadie menciona y que en este trabajo importó.** Si en vez de
resolver con `A` uno resuelve con `AᵀA` (las llamadas *ecuaciones normales*), el
número de condición **se eleva al cuadrado**. De `10⁸` se pasa a `10¹⁶`, es decir a
basura. Eso fue exactamente lo que arruinó el método de Steihaug-CG en nuestras
pruebas (→ §6.7).

### 2.6. Métodos directos e iterativos

Para resolver un sistema lineal `A x = b` hay dos familias:

**Directos** (eliminación de Gauss, factorización LU). Hacen un número fijo de
operaciones y devuelven la respuesta exacta (salvo redondeo). Cuesta
`≈ n³/3` operaciones y hay que guardar la matriz entera.

```
    n = 1 000      →  3×10⁸ operaciones     → menos de un segundo
    n = 100 000    →  3×10¹⁴ operaciones    → varios días
    n = 1 000 000  →  3×10¹⁷ operaciones    → décadas
```

**Iterativos** (Jacobi, Gauss-Seidel, y los de Krylov que veremos). Empiezan con una
aproximación y la van mejorando. Se pueden **cortar cuando uno quiera**: si con 20
iteraciones ya alcanza la precisión que hace falta, se paran ahí.

**Para qué nos importa.** Los dos rasgos de los iterativos son justo lo que
necesitamos: se cortan a voluntad (Problema 1) y —los de Krylov— solo necesitan
multiplicar por la matriz, nunca verla (Problema 2).

### 2.7. Orden de convergencia

**El problema que resuelve:** dos métodos pueden «converger» y sin embargo uno tardar
5 pasos y el otro 5 000. Hace falta clasificar la velocidad.

Sea `e_k = ‖x_k − x*‖` el error en el paso `k`. Se dice que la convergencia es:

```
    LINEAL         si    e_{k+1} ≈ C · e_k        con 0 < C < 1
    SUPERLINEAL    si    e_{k+1} / e_k  →  0
    CUADRÁTICA     si    e_{k+1} ≈ C · e_k²
```

La diferencia entre lineal y cuadrática es abismal. Vea las **cifras correctas** que
gana cada paso, empezando ambos con 1 cifra correcta:

| Paso | Lineal (C = 0.1) | Cuadrática |
|---|---|---|
| 0 | 1 cifra | 1 cifra |
| 1 | 2 cifras | 2 cifras |
| 2 | 3 cifras | 4 cifras |
| 3 | 4 cifras | 8 cifras |
| 4 | 5 cifras | **16 cifras — límite de la máquina** |
| 10 | 11 cifras | (hace rato terminó) |

La convergencia cuadrática **duplica las cifras correctas en cada paso**. Por eso
Newton es el método de referencia: cuando arranca cerca, llega a precisión de máquina
en cuatro o cinco iteraciones.

Lo verá con números reales en §3.2.

---

## 3. El método de Newton

### 3.1. En una dimensión: la idea de la tangente

Queremos resolver `f(x) = 0` con una sola incógnita. La idea de Newton es:

> Si no sé resolver `f`, la reemplazo por la recta que mejor la imita cerca del punto
> donde estoy — la tangente — y resuelvo **esa**, que sí sé.

```
       f(x)
        │        ⟋ tangente en x_k
        │      ⟋
        │    ⟋           ● (x_k , f(x_k))
        │  ⟋           ⟋
    ────┼─⟋──────────●──────────────────── x
        │ x_{k+1}   x_k
        │
        └─ la tangente cruza el eje en x_{k+1}: ese es el nuevo candidato
```

La recta tangente en `x_k` tiene ecuación `y = f(x_k) + f'(x_k)·(x − x_k)`.
Igualando a cero y despejando `x`:

```
    0 = f(x_k) + f'(x_k)·(x − x_k)
    f'(x_k)·(x − x_k) = −f(x_k)
    x − x_k = −f(x_k)/f'(x_k)
    ────────────────────────────────
    x_{k+1} = x_k − f(x_k)/f'(x_k)          ← la fórmula de Newton
```

Fíjese en la estructura, porque es la que se generaliza:

```
    (derivada) · (paso)  =  −(residuo)
    luego:  nuevo punto = punto viejo + paso
```

### 3.2. Ejemplo a mano: calcular √2

Resolvamos `f(x) = x² − 2 = 0`, que tiene `f'(x) = 2x`. Empezamos en `x₀ = 1`:

```
    x₁ = 1 − (1² − 2)/(2·1)         = 1 + 1/2               = 1.5
    x₂ = 1.5 − (1.5² − 2)/(2·1.5)   = 1.5 − 0.25/3          = 1.416666666666667
    x₃ = 1.416666… − 0.00694…/2.833… = 1.414215686274510
```

Y la tabla completa, con el error contra el valor verdadero `√2 = 1.41421356237309…`:

| k | x_k | f(x_k) | error respecto de √2 | cifras correctas |
|---|---|---|---|---|
| 0 | 1.000000000000000 | −1.0·10⁰ | 4.1·10⁻¹ | 0 |
| 1 | 1.500000000000000 | +2.5·10⁻¹ | 8.6·10⁻² | 1 |
| 2 | 1.416666666666667 | +6.9·10⁻³ | 2.5·10⁻³ | 2 |
| 3 | 1.414215686274510 | +6.0·10⁻⁶ | 2.1·10⁻⁶ | 5 |
| 4 | 1.414213562374690 | +4.5·10⁻¹² | 1.6·10⁻¹² | 11 |
| 5 | 1.414213562373095 | +4.4·10⁻¹⁶ | 0 | **16 (exacto)** |

Mire la columna del error: `10⁻¹ → 10⁻² → 10⁻³ → 10⁻⁶ → 10⁻¹² → 0`. Los exponentes
se van **duplicando**. Eso es la convergencia cuadrática de §2.7, vista con números.
En cinco pasos llegó al límite de la aritmética de la computadora.

Con esto en la mano se entiende por qué Newton es el método de referencia. Ahora
viene la letra chica.

### 3.3. En n dimensiones: dónde aparece el sistema lineal

Con muchas incógnitas la idea es idéntica, pero la «recta tangente» se convierte en
una **aproximación lineal** hecha con el Jacobiano. El desarrollo de Taylor de
primer orden dice:

```
    F(x + s)  ≈  F(x) + J(x)·s
```

(léalo así: si me muevo `s` desde `x`, el residuo cambia aproximadamente `J·s`.)

Queremos que el nuevo residuo sea cero. Igualamos y despejamos:

```
    F(x) + J(x)·s = 0
    ────────────────────────────
    J(x_k) · s_k = −F(x_k)          ← ESTO ES UN SISTEMA LINEAL
    x_{k+1} = x_k + s_k
```

Compare con §3.1: donde antes había una división `f/f'`, ahora hay que **resolver un
sistema lineal**, porque no se puede «dividir» por una matriz.

Y acá está el punto donde nace todo el tema:

> **En cada iteración de Newton hay que resolver un sistema lineal de `n × n`.**
> Si `n` es un millón, esa frase esconde un problema gigantesco. Los Problemas 1 y 2
> de §0.3 son exactamente ese problema.

Al vector `s_k` se le llama **paso de Newton** (o *dirección de Newton*).

### 3.4. Ejemplo a mano: el sistema 2×2 desde un punto malo

Retomemos el sistema de §1.3 desde `x₀ = (2, 0.5)`. Ya calculamos en §1.3 y §2.2:

```
    F(x₀) = (2.250000, 0.843282)        ‖F(x₀)‖ = 2.402837

            ⎡ 4.0000   1.0000 ⎤
    J(x₀) = ⎢                 ⎥         det = 0.2817
            ⎣ 2.7183   0.7500 ⎦
```

Resolvemos `J·s = −F`, es decir

```
    4.0000·s₁ + 1.0000·s₂ = −2.250000
    2.7183·s₁ + 0.7500·s₂ = −0.843282
```

La solución de ese sistema 2×2 es:

```
    s = (−2.996676 , +9.736705)
```

**Mire el tamaño del paso.** Estamos parados en `(2, 0.5)` y Newton nos manda a:

```
    x₁ = x₀ + s = (2 − 2.9967 ,  0.5 + 9.7367) = (−0.996676 , 10.236705)
```

¡Nos manda a `x₂ ≈ 10`, cuando la solución está en `(1, 1)`! ¿Por qué? Porque
`det(J) = 0.2817` es chico: la matriz está casi al borde de no ser invertible, y al
«dividir» entre algo casi cero el paso se dispara.

Y las cosas empeoran. Siguiendo con Newton puro desde ahí:

| k | x_k | ‖F(x_k)‖ |
|---|---|---|
| 0 | (2.000000, 0.500000) | 2.4 |
| 1 | (−0.996676, 10.236705) | 1.08·10³ |
| 2 | (16.007064, 6.823056) | 3.29·10⁶ |
| 3 | … | desborda |

En 16 iteraciones el residuo se hace `inf` y el programa se rinde. Este es
**exactamente** el resultado del ejercicio 1 (`ej1_newton_vs_globalizado.py`, tabla
de la Parte B): *Newton puro — NO converge — divergió (‖F‖ desbordó)*.

**Y ahora la solución, adelantada.** ¿Qué pasa si en vez de dar el paso completo
damos solo una fracción `λ` de él? Probemos, en el primer paso:

```
   partimos de  ‖F(x₀)‖ = 2.4028  y queremos BAJAR de ahí

    λ = 1      → x = (−0.9967, 10.2367)   ‖F‖ = 1075.86    peor
    λ = 1/2    → x = ( 0.5017,  5.3684)   ‖F‖ =  155.69    peor
    λ = 1/4    → x = ( 1.2508,  2.9342)   ‖F‖ =   25.87    peor
    λ = 1/8    → x = ( 1.6254,  1.7171)   ‖F‖ =    6.10    peor
    λ = 1/16   → x = ( 1.8127,  1.1085)   ‖F‖ =    2.9894  peor
    λ = 1/32   → x = ( 1.9064,  0.8043)   ‖F‖ =    2.4888  peor
    λ = 1/64   → x = ( 1.9532,  0.6521)   ‖F‖ =    2.4037  peor por poquito
    λ = 1/128  → x = ( 1.9766,  0.5761)   ‖F‖ =    2.3935  ← ¡por fin baja!
```

Dos cosas para notar en esa escalera. La primera es lo **violento** que era el paso
completo: multiplicaba el residuo por 450. La segunda es que en `λ = 1/64` faltó muy
poco (2.4037 contra 2.4028): la mejora aparece de golpe, no de a poco.

Con `λ = 1/128` el residuo por fin baja. Eso es **búsqueda de línea con retroceso**,
y es la primera de las tres estrategias del tema (→ §6.6). Con ella, el mismo
problema desde el mismo punto converge en **7 iteraciones** a `(1, 1)` con
`‖F‖ = 9.9·10⁻¹²`.

### 3.5. El teorema de convergencia, traducido frase por frase

El enunciado formal dice:

> Sea `F` continuamente diferenciable en un entorno de `x*`, con `F(x*) = 0`, `J(x*)`
> no singular y `J` Lipschitz continua cerca de `x*`. **Entonces existe `δ > 0` tal
> que** para todo `x₀` con `‖x₀ − x*‖ < δ`, la sucesión de Newton converge a `x*` y
> satisface `‖x_{k+1} − x*‖ ≤ C‖x_k − x*‖²`.

Desarmémoslo:

| Trozo | Qué significa en castellano |
|---|---|
| «`F` continuamente diferenciable» | `F` es suave: se puede derivar y la derivada no da saltos. |
| «`J(x*)` no singular» | En la solución, el Jacobiano es invertible (determinante ≠ 0). Si fuera singular, el sistema lineal no tendría solución única. |
| «`J` Lipschitz continua» | El Jacobiano no cambia demasiado brusco de un punto a otro. Es una condición técnica de suavidad. |
| «**existe `δ > 0` tal que**» | ⚠ **Aquí está todo el asunto.** Ver abajo. |
| «`‖x₀ − x*‖ < δ`» | El punto inicial tiene que estar a distancia menor que `δ` de la solución. |
| «`≤ C‖x_k − x*‖²`» | Convergencia cuadrática: las cifras se duplican (§2.7). |

**Las cuatro palabras que originan el tema de la exposición:** *«existe δ > 0 tal
que»*.

El teorema afirma que **existe** un radio de convergencia. Pero:

- **No dice cuánto vale `δ`.**
- **No da ninguna forma de calcularlo ni de estimarlo.**
- **En problemas reales `δ` puede ser ridículamente pequeño.**
- Y **fuera de esa bola el teorema no promete absolutamente nada** — y «nada»
  incluye que la sucesión diverja, como acabamos de ver en §3.4.

O sea: Newton es un método **local**. Es rapidísimo *si ya estaba cerca*. El problema
es que uno normalmente no sabe dónde está la solución (por eso la busca), así que no
tiene manera de saber si su punto inicial está dentro de la bola o no.

**A eso apunta la palabra «global» del título de la exposición:** conseguir que el
método funcione desde puntos iniciales *arbitrarios*, no solo desde puntos que ya
estaban cerca.

### 3.6. El caso más simple donde Newton falla: arctan

No hace falta un sistema grande para romper a Newton. Alcanza con una función de una
variable, suave, monótona y sin nada raro: `f(x) = arctan(x)`, cuya raíz es `x* = 0`.

El paso de Newton es `x_{k+1} = x_k − arctan(x_k)·(1 + x_k²)`. Desde `x₀ = 2`:

| k | x_k | arctan(x_k) |
|---|---|---|
| 0 | +2.0000 | +1.1071 |
| 1 | −3.5357 | −1.2952 |
| 2 | +13.9510 | +1.4992 |
| 3 | −279.3441 | −1.5672 |
| 4 | +122 017.00 | +1.5708 |
| 5 | −2.34·10¹⁰ | −1.5708 |
| 6 | +8.59·10²⁰ | +1.5708 |

Salta de un lado al otro, cada vez más lejos. ¿Por qué? Porque `arctan` **se
aplana** cuando `x` crece: la tangente es casi horizontal y cruza el eje muy lejos.

```
     arctan(x)
        │      ⋯⋯⋯⋯⋯⋯⋯⋯⋯  ← se aplana: la tangente es casi horizontal
   +1.57│    ⋯
        │  ⟋●  x₀ = 2
   ─────┼─⟋───────────────── x
       ⟋│
     ⟋  │
   ●    │  x₁ = −3.54    ← la tangente cruzó el eje mucho más lejos
```

El umbral exacto es `|x₀| ≈ 1.3917`: por debajo converge, por encima diverge. Nadie
que resuelva un problema real sabe de antemano dónde está su umbral.

**Y el arreglo es el mismo de antes.** Con búsqueda de línea, desde el mismo `x₀ = 2`:

```
    x:  2  →  −0.768  →  0.273  →  −0.0134  →  1.6·10⁻⁶  →  −2.75·10⁻¹⁸
    λ:     0.5       1.0        1.0        1.0         1.0
```

Solo el primer paso necesitó recortarse a la mitad. A partir de ahí ya estaba dentro
de la bola de convergencia y aceptó pasos completos, recuperando la velocidad
cuadrática. Converge en 5 iteraciones.

**Esa es la idea entera de la globalización en una frase:** recortar los primeros
pasos lo justo para entrar en la zona donde Newton funciona, y después dejarlo
correr a toda velocidad.

### 3.7. Resumen: los tres problemas de Newton

Ya podemos nombrar con precisión lo que hay que arreglar:

| # | Problema | Aparece en | Se arregla en |
|---|---|---|---|
| 1 | Resolver `J s = −F` exacto cuesta `O(n³)` | §3.3, §2.6 | §4 — Newton inexacto |
| 2 | Ni siquiera se puede formar `J` | §1.5, §2.2 | §5 — Krylov matriz-libre |
| 3 | Solo converge cerca de la raíz | §3.4, §3.5, §3.6 | §6 — Globalización |

Los tres tienen solución, y las tres soluciones se combinan sin estorbarse. El
resultado es el método **Newton-Krylov globalizado**.

---

## 4. Problema 1 — el sistema lineal cuesta demasiado

### 4.1. Cuánto cuesta resolverlo exacto

En cada iteración de Newton hay que resolver `J s = −F`. Con eliminación de Gauss o
factorización LU eso son `≈ n³/3` operaciones (§2.6). Y hay que hacerlo **en cada
iteración**, porque `J` cambia con `x`.

```
    n = 1 000        →  3×10⁸   operaciones   →  ~0.3 s      × 10 iteraciones = 3 s
    n = 100 000      →  3×10¹⁴  operaciones   →  ~4 días     × 10 = imposible
    n = 1 000 000    →  3×10¹⁷  operaciones   →  ~10 años    × 10 = ni hablar
```

Hay que buscar otra cosa.

### 4.2. La idea: resolver «más o menos»

Acá viene una observación que parece trivial y no lo es.

Recuerde de §3.3 que `J s = −F` **no es el problema que queremos resolver**. Es una
*aproximación* al problema: viene de reemplazar `F` por su recta tangente. Esa
aproximación es buena si estamos cerca de la solución, y **es una mentira si estamos
lejos**.

Entonces:

> Si estamos lejos de la raíz, el modelo lineal ni siquiera describe bien a `F`.
> ¿Para qué gastar horas resolviendo con dieciséis cifras un modelo que solo vale
> una o dos?

La propuesta es **resolver el sistema lineal solo aproximadamente**, y de paso usar
un método iterativo (§2.6), que se puede cortar cuando uno quiera.

### 4.3. El criterio de Newton inexacto

Si resolvemos el sistema solo a medias, el paso `s` que obtenemos no cumple
`J s = −F` exactamente. Le sobra un poquito, que llamamos **residuo lineal**:

```
    r  =  J s + F           (si el paso fuera exacto, r sería cero)
```

Necesitamos una regla que diga cuán grande podemos dejar ese sobrante. La regla es:

```
    ‖J s + F‖  ≤  η · ‖F‖
        ↑           ↑   ↑
        │           │   └─ el residuo NO lineal actual: cuán lejos estamos
        │           └───── el "término de forzado", un número entre 0 y 1
        └───────────────── lo que le sobró al sistema lineal
```

En castellano: **«el error que dejo en el sistema lineal tiene que ser una fracción
`η` de lo mal que estoy ahora»**.

Y eso es astuto, porque la exigencia se **autoajusta**:

- Al principio, cuando `‖F‖` es grande, la cota `η‖F‖` es grande: se permite un paso
  muy chapucero, que sale barato.
- Cerca del final, cuando `‖F‖` es diminuto, la cota se vuelve diminuta: se exige un
  paso muy preciso, pero para entonces cada resolución cuesta poco porque el iterado
  ya está cerca.

Al número `η ∈ [0, 1)` se le llama **término de forzado** (*forcing term* en inglés).
Con `η = 0` se recupera Newton exacto.

### 4.4. El teorema de Dembo, Eisenstat y Steihaug (1982)

La pregunta obvia es: *si resuelvo mal a propósito, ¿no arruino la convergencia
cuadrática que era toda la gracia de Newton?* La respuesta es el teorema fundamental
del tema:

| Si se elige… | se obtiene… |
|---|---|
| `η_k ≤ η < 1` constante | convergencia **lineal** |
| `η_k → 0` | convergencia **superlineal** |
| `η_k = O(‖F_k‖)` | convergencia **cuadrática** (¡Newton completo!) |

Léalo despacio, porque dice algo notable:

> **Se puede ser perezoso sin pagar absolutamente nada en velocidad de
> convergencia**, siempre que la pereza se apriete al ritmo correcto.

La tercera fila es la clave: si `η` se hace proporcional a `‖F‖` —o sea, si a medida
que nos acercamos exigimos proporcionalmente más precisión— se recupera la
convergencia cuadrática íntegra, **sin haber resuelto exactamente ni un solo sistema
lineal en todo el proceso**.

### 4.5. Oversolving: resolver de más sale caro

Si un `η` chico da mejor convergencia, la tentación es poner `η = 10⁻¹²` y olvidarse.
Es un error, y tiene nombre: **oversolving** (resolver de más).

Lo medimos. El ejercicio 2 (`ej2_bratu1d_newton_krylov.py`) resuelve la ecuación de
Bratu con `n = 250`, cambiando **únicamente** el valor de `η`:

| η | Iteraciones de Newton | Productos J·v | Tiempo |
|---|---|---|---|
| 10⁻¹ fijo | 7 | 9 854 | 1.8 s |
| 10⁻³ fijo | 4 | 21 332 | 3.9 s |
| **10⁻¹² fijo** | **4** | **24 804** | **4.4 s** |
| Eisenstat-Walker 2 | 8 | **9 707** | 1.8 s |

Mire la fila resaltada. Resolver con doce cifras usa **la mitad** de iteraciones de
Newton… y **2.6 veces más trabajo total**.

**Por qué pasa esto.** Cada iteración de Newton es barata o cara según cuántas
iteraciones internas necesite el solver lineal. Con `η` chico, cada paso de Newton
obliga al solver a trabajar muchísimo. Se cambian «pocas iteraciones caras» por
«muchas iteraciones baratas», y las baratas ganan.

**La lección metodológica, que vale para cualquier trabajo numérico:**

> Contar iteraciones externas **engaña**. Lo que se paga son las evaluaciones de la
> función. En este proyecto todo se mide en **productos `J·v`**, porque —como veremos
> en §5.4— cada producto `J·v` cuesta exactamente una evaluación de `F`.

> **Analogía.** Es como corregir un examen. Si el alumno va por la pregunta 1 de 20,
> no tiene sentido corregir con lupa cada coma: se corrige por encima y se sigue.
> Recién al final, cuando ya casi está la nota, vale la pena mirar en detalle.

### 4.6. Eisenstat-Walker: elegir η automáticamente

Queda una pregunta práctica: ¿quién elige `η` en cada paso? Eisenstat y Walker (1996)
propusieron hacerlo automático mirando **cuánto se equivocó el modelo lineal en el
paso anterior**. La versión que usamos (su *Choice 2*) es:

```
                ⎛  ‖F_k‖   ⎞^α
    η_k  =  γ · ⎜ ───────  ⎟          con  γ = 0.9   y   α = (1+√5)/2 ≈ 1.618
                ⎝ ‖F_{k−1}‖ ⎠
```

Desarmada:

- **`‖F_k‖ / ‖F_{k−1}‖`** es la razón entre el residuo de ahora y el de antes. Si vale
  0.01, el residuo bajó cien veces: el modelo lineal está funcionando muy bien.
  Si vale 0.9, apenas bajó: el modelo no sirve de mucho.
- **Elevar a `α ≈ 1.618`** exagera esa señal: si la razón es chica, `η` se hace todavía
  más chico.
- **`γ = 0.9`** es un factor de seguridad, para no apretar del todo.

En una frase: **«si el paso anterior salió muy bien, apretá; si salió regular,
aflojá».**

Además lleva una **salvaguarda**: si `γ·η_{k−1}^α > 0.1`, se toma
`η_k = max(η_k, γ·η_{k−1}^α)`. Impide que `η` se desplome por un único paso
afortunado que después no se sostenga.

**Cómo se ve funcionando.** Estos son los valores reales que eligió en el ejercicio 2:

| k | η elegido | ‖F_k‖ |
|---|---|---|
| 0 | 9.0·10⁻¹ | 7.5·10⁻⁴ |
| 1 | 7.6·10⁻¹ | 1.7·10⁻⁴ |
| 2 | 5.8·10⁻¹ | 1.2·10⁻⁴ |
| 3 | 3.7·10⁻¹ | 1.4·10⁻⁵ |
| 4 | 1.8·10⁻¹ | 1.6·10⁻⁶ |
| 5 | 9.2·10⁻³ | 5.4·10⁻⁸ |

Empieza en 0.9 —o sea, casi sin resolver nada— y termina en 0.0092. **Nadie se lo
dijo:** lo dedujo solo, mirando cómo iba bajando el residuo. Y con eso fue la opción
más barata de las cinco que probamos.

---

## 5. Problema 2 — la matriz ni siquiera cabe

### 5.1. Antes de empezar: subespacio y `span`

Dos palabras que van a aparecer y conviene tener claras.

Un **subespacio** de `Rⁿ` es un conjunto de vectores cerrado bajo sumas y
multiplicación por números: si `u` y `v` están, entonces `3u − 2v` también está.
Geométricamente, en `R³`, los subespacios son: el origen, las rectas por el origen,
los planos por el origen, y todo `R³`.

El **`span`** (o *envolvente lineal*) de un conjunto de vectores es el subespacio de
todas sus combinaciones:

```
    span{ v₁, v₂ }  =  { a·v₁ + b·v₂   para todos los números a, b }
```

Si `v₁` y `v₂` no son paralelos, `span{v₁, v₂}` es el plano que los contiene.

### 5.2. El subespacio de Krylov

Los métodos de Krylov resuelven `A s = b` buscando la mejor solución posible dentro
de un subespacio que se va agrandando:

```
    𝒦₁ = span{ b }
    𝒦₂ = span{ b, Ab }
    𝒦₃ = span{ b, Ab, A²b }
      ⋮
    𝒦ₘ = span{ b, Ab, A²b, …, A^(m−1)b }
```

Cada iteración agrega una potencia más de `A` aplicada a `b`, y busca la mejor
respuesta ahí adentro. Cuanto más grande el subespacio, mejor la respuesta.

**¿Por qué esas potencias y no otros vectores?** Porque `b`, `Ab`, `A²b`, … contienen
exactamente la información que la matriz «revela» sobre `b`. De hecho hay un teorema
(Cayley-Hamilton) que garantiza que la solución exacta de `A s = b` **está** en `𝒦ₙ`.
En la práctica se llega mucho antes: con unas decenas de iteraciones ya alcanza.

**Y acá está la observación que cambia todo:**

> Para construir `b, Ab, A²b, …` lo único que se necesita es **multiplicar por `A`**.
> Nunca hace falta ver la matriz entrada por entrada, ni factorizarla, ni siquiera
> almacenarla.
>
> **Basta con tener una función que, dado un vector `v`, devuelva `A·v`.**

Guarde esa frase. Es la bisagra de toda la sección.

### 5.3. GMRES, y por qué no gradientes conjugados

De los métodos de Krylov, el que usamos es **GMRES** (*Generalized Minimal
RESidual*). En cuatro frases:

1. Construye una base ortonormal del subespacio `𝒦ₘ` (procedimiento de Arnoldi).
2. Dentro de ese subespacio busca el `s` que **minimiza `‖A s − b‖`**.
3. Cada iteración agranda el subespacio en uno y mejora esa mejor respuesta.
4. Como guarda todos los vectores de la base, la memoria crece con `m`; por eso se
   usa con **reinicio** (se reinicia cada 30 o 40 iteraciones).

Comparado con las alternativas:

| Método | Requisito sobre `A` | Memoria |
|---|---|---|
| **GMRES** | ninguno | crece con `m` (por eso el reinicio) |
| BiCGSTAB | ninguno | fija, pero menos estable |
| TFQMR | ninguno | fija |
| Gradientes conjugados (CG) | **simétrica y definida positiva** | fija, el más barato |

**Por qué no CG, que sería más barato.** Porque CG exige que la matriz sea simétrica
definida positiva, y el Jacobiano de un problema real **no lo es**. En cuanto hay
convección, transporte o acoplamiento asimétrico entre variables, `J ≠ Jᵀ`. En
nuestros problemas de difusión-reacción concretos sí resulta simétrica —y lo
aprovechamos en §6.7— pero no se puede suponer en general.

Es una pregunta que el docente puede hacer perfectamente. La respuesta corta:
**«porque J no es simétrica»**.

### 5.4. El truco central: `J·v` sin construir `J`

Recapitulemos. GMRES solo necesita poder calcular `J·v`. Y `J·v` resulta ser algo
que ya sabemos aproximar.

Vuelva a la definición de derivada de toda la vida, en una dimensión:

```
    f'(x)  =  lím       f(x + ε) − f(x)
              ε→0       ───────────────
                              ε
```

Si no tomamos el límite y nos quedamos con un `ε` pequeño pero finito, tenemos una
aproximación. En varias variables, la versión direccional de eso es exactamente
`J·v`:

```
    ┌──────────────────────────────────────────────┐
    │                     F(x + ε·v) − F(x)        │
    │      J(x)·v   ≈    ─────────────────────     │
    │                              ε               │
    └──────────────────────────────────────────────┘
```

Léalo con calma, porque es el corazón del método:

- El lado izquierdo requiere, en principio, **conocer la matriz `J` entera**.
- El lado derecho requiere **una sola evaluación de `F`** (ya tenemos `F(x)` guardado
  del paso anterior, así que solo hay que calcular `F(x + εv)`).

**El Jacobiano nunca se construye. Ni una sola entrada.** A esto se le llama
*Jacobian-Free Newton-Krylov*, **JFNK**.

Los números de nuestro ejercicio 2, con `n = 250`:

```
    Entradas de J si la formáramos ............... 62 500
    Evaluaciones de F para formarla por columnas ..   250
    Evaluaciones de F para UN producto J·v .......      1
    Productos J·v por paso de Newton, con precondicionador ...  ~4
    Productos J·v por paso de Newton, sin precondicionador  ... ~1200
```

Con `n = 10⁶` la primera fila es 10¹² y la segunda 10⁶, mientras que la tercera sigue
siendo 1. Ahí se ve por qué esto es lo único viable.

**El precio a pagar.** La fórmula es una *aproximación*, no la derivada exacta. El
Jacobiano queda conocido con unas 8 cifras significativas en vez de 16. Consecuencias
concretas:

- Cerca de la raíz, la convergencia deja de ser exactamente cuadrática.
- No se puede pedir un residuo final mucho menor que `‖F(x₀)‖·√εₘ`.
- Ciertos algoritmos que amplifican el error dejan de funcionar (lo veremos en §6.7).

### 5.5. Elegir el paso ε: dos errores que tiran para lados opuestos

¿Cuánto vale `ε`? La tentación es «lo más chico posible, así se parece más al
límite». Es un error, y entenderlo es de lo más instructivo del tema.

Hay **dos** fuentes de error compitiendo:

```
  ┌─ ERROR DE TRUNCAMIENTO ───────────────┐   ┌─ ERROR DE CANCELACIÓN ────────────┐
  │ La fórmula es la derivada MÁS un      │   │ F(x+εv) y F(x) se parecen cada    │
  │ término de orden ε (Taylor).          │   │ vez más; al restarlos se pierden  │
  │                                       │   │ cifras (§2.4) y después se divide │
  │      ≈ (ε/2)·‖F″‖·‖v‖²                │   │ entre ε, que amplifica.           │
  │                                       │   │      ≈ εₘ·‖F‖ / (ε·‖v‖)           │
  │  ↓ BAJA cuando ε baja                 │   │  ↑ SUBE cuando ε baja             │
  └───────────────────────────────────────┘   └───────────────────────────────────┘
```

El error total es la suma, y se minimiza donde los dos se cruzan. Igualándolos y
despejando:

```
                       ┌──────────────┐
                       │  εₘ · ‖F‖    │
    ε_óptimo  ≈   √    │ ───────────  │        y si ‖F‖ ≈ ‖F″‖:   ε ≈ √εₘ ≈ 1.5·10⁻⁸
                       │    ‖F″‖      │
                       └──────────────┘
```

De ahí sale la receta famosa `ε ≈ √εₘ ≈ 1.5 × 10⁻⁸`.

**Medido de verdad.** Esta es la tabla que produce `clase2_matrix_free.py`,
comparando el `J·v` aproximado contra el Jacobiano analítico exacto:

| ε | error relativo | quién manda |
|---|---|---|
| 10⁻¹⁶ | 9.5·10⁻¹ (¡95 %!) | cancelación |
| 10⁻¹⁴ | 2.9·10⁻² | cancelación |
| 10⁻¹² | 2.8·10⁻⁴ | cancelación |
| 10⁻¹⁰ | 3.0·10⁻⁶ | cancelación |
| 10⁻⁸ | 2.7·10⁻⁸ | cancelación |
| 10⁻⁶ | 2.9·10⁻¹⁰ | ← mínimo |
| 10⁻⁴ | 5.7·10⁻⁹ | truncamiento |
| 10⁻² | 5.7·10⁻⁷ | truncamiento |

Dibujada, la curva tiene forma de **V**:

```
  error
   1e+0 │●
        │ ●                                            ● ← truncamiento
   1e-4 │   ●                                    ●
        │     ●                            ●
   1e-8 │       ●                    ●
        │          ●           ●
   1e-10│              ●  ●  ●          ← el mínimo
        └──┴────┴────┴────┴────┴────┴────┴────┴───── ε
         1e-16 1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2
         └── cancelación ──┘        └─ truncamiento ─┘
```

**Un detalle honesto que salió midiendo.** El mínimo medido cayó en `2.5·10⁻⁶`, dos
órdenes por encima de `√εₘ = 1.5·10⁻⁸`. No es un error del experimento: en ese
problema concreto el operador está dominado por su parte lineal, así que `‖F″‖` es
chica y —según la fórmula de arriba— el óptimo se corre hacia arriba.

¿Importa? **No.** Usar `√εₘ` da un error de `1.9·10⁻⁸` en vez del óptimo `1.8·10⁻¹⁰`:
cien veces peor y perfectamente irrelevante, porque Newton inexacto ya tolera un
residuo lineal de `η‖F‖` con `η ~ 10⁻²` (§4.3). **Un Jacobiano con 7 cifras correctas
le sobra.** Por eso `√εₘ` se usa siempre: no es óptimo, pero nunca es catastrófico y
no exige conocer `‖F″‖`.

### 5.6. Precondicionamiento: donde se gana el 90 % del rendimiento

**El problema.** GMRES converge rápido cuando los autovalores de la matriz están
**agrupados**. Si están muy dispersos —o sea, si el número de condición es grande
(§2.5)— tarda muchísimo. Y ya sabemos que el Jacobiano de una EDP tiene
`cond ≈ O(n²)`: cuanto más fina la malla, peor.

**La idea.** En vez de resolver `J s = −F`, resolver un sistema equivalente pero
mejor portado:

```
    M⁻¹ J s  =  −M⁻¹ F
```

Si `M` se parece a `J`, entonces `M⁻¹J` se parece a la identidad, cuyos autovalores
están todos en 1: agrupadísimos. GMRES lo resuelve en poquísimas iteraciones.

**La tensión.** `M` tiene que cumplir dos cosas a la vez, y son contradictorias:

```
    parecerse a J   ←────── tensión ──────→   ser BARATO de invertir
    (si M = J, perfecto…)                      (…pero invertir J era el problema)
```

Opciones habituales: la parte lineal del operador · una factorización incompleta
(ILU) · multigrid · descomposición de dominios · el Jacobiano **congelado** de alguna
iteración anterior, reutilizado varias veces.

**Lo que ganó en nuestro caso.** En Bratu, el Jacobiano es
`J = L − h²λ·diag(e^u)`, donde `L` es el laplaciano discreto. El término no lineal es
una perturbación diagonal pequeña, así que `L` solo ya es una aproximación excelente
— y `L` es **tridiagonal**, o sea que se factoriza una sola vez en `O(n)`.

| η | J·v sin precondicionador | J·v con precondicionador | ganancia |
|---|---|---|---|
| 10⁻¹ fijo | 9 854 | 19 | **519×** |
| 10⁻³ fijo | 21 332 | 24 | **889×** |
| 10⁻¹² fijo | 24 804 | 24 603 | 1.0× |
| Eisenstat-Walker 2 | 9 707 | 22 | **441×** |

Dos comentarios sobre esa tabla, y los dos son buen material de exposición:

1. **La ganancia es enorme porque el precondicionador es muy bueno para *este*
   problema.** No hay precondicionador universal: elegirlo exige conocer el operador.
   En un problema dominado por convección, `L` solo no alcanzaría.

2. **La fila de `η = 10⁻¹²` no mejora nada.** Y tiene explicación: le estamos pidiendo
   a GMRES un residuo relativo de `10⁻¹²`, que está por debajo de lo que la aritmética
   de doble precisión puede entregar sobre este operador. No puede alcanzarlo nunca,
   agota su presupuesto de iteraciones en cada paso, y el precondicionador no lo
   salva. **Pedir más precisión de la que existe no acelera: quema trabajo.**

---

## 6. Problema 3 — la convergencia global

Esta es la sección del tema de la exposición. Todo lo anterior era preparación.

### 6.1. Primero: «global» no significa lo que parece

Aclaremos esto antes que nada, porque es la confusión más común y da mala impresión
en una exposición:

| Frase | Qué significa **acá** |
|---|---|
| «convergencia **global**» | converger **desde cualquier punto inicial**, no solo desde uno cercano |
| «mínimo **global**» | el punto más bajo de toda la función |

**Son cosas distintas y no relacionadas.** La convergencia global de la que hablamos
*no* garantiza encontrar el mínimo global de nada. De hecho, en §6.9 vamos a ver que
puede quedarse atrapada en un mínimo *local*.

### 6.2. La función de mérito

**El problema.** El método de Newton no tiene ninguna noción interna de «estar
mejorando». Calcula un paso y lo acepta, punto. Si el resultado es peor que el punto
de partida —como en §3.4— no se entera.

Hay que darle una **medida de progreso**. Pero `F(x)` es un vector de `n` números:
no se pueden comparar dos vectores directamente. Necesitamos convertirlo en un
número. Se llama **función de mérito** y la elección estándar es:

```
    f(x)  =  ½ ‖F(x)‖²  =  ½ · ( F₁(x)² + F₂(x)² + … + Fₙ(x)² )
```

**Por qué esta y no otra:**

- Es **escalar**: un solo número, comparable.
- Es **no negativa** siempre.
- Vale **cero exactamente en las raíces de `F`**, y solo ahí. Así que minimizar `f` y
  resolver `F(x) = 0` son la misma búsqueda.
- El **cuadrado** la hace diferenciable en todos lados (la norma sola no lo es en el
  origen, igual que `|x|` no es derivable en 0).
- El **½** está solo para que al derivar se cancele el 2 y la fórmula quede limpia.
  No tiene ningún significado.

Con esta traducción, «resolver un sistema no lineal» se convierte en «minimizar una
función escalar», y toda la maquinaria de la optimización queda disponible.

### 6.3. De dónde sale `∇f = JᵀF`

Vamos a necesitar el gradiente de `f`. Sale con la regla de la cadena, y conviene
verlo en detalle porque después lo usamos en la demostración clave.

Escribamos `f` como suma:

```
    f(x) = ½ · Σⱼ  Fⱼ(x)²                        (suma sobre j = 1 … n)
```

Derivemos respecto de una incógnita cualquiera `xᵢ`. Regla de la cadena sobre cada
término del sumatorio:

```
    ∂f/∂xᵢ  =  ½ · Σⱼ  2 · Fⱼ(x) · ∂Fⱼ/∂xᵢ  =  Σⱼ  Fⱼ(x) · ∂Fⱼ/∂xᵢ
```

(acá se ve para qué servía el ½: canceló el 2.)

Ahora observe qué es `∂Fⱼ/∂xᵢ`: es la entrada de la **fila j, columna i** del
Jacobiano (§2.2). Y la suma `Σⱼ Fⱼ · J_{ji}` es exactamente la componente `i` del
producto `Jᵀ F` — el traspuesto aparece porque estamos recorriendo la columna `i`
bajando por las filas `j`. Entonces:

```
    ┌────────────────────────┐
    │   ∇f(x)  =  J(x)ᵀ F(x) │
    └────────────────────────┘
```

Guárdela: la usamos en §6.5 y en §6.7.

### 6.4. Qué es una dirección de descenso, en este contexto

De §2.3: la dirección `s` hace bajar a `f` si `∇f(x)ᵀ s < 0`.

Aplicado a nuestro caso, con `∇f = JᵀF`:

```
    ∇f(x)ᵀ s  =  (JᵀF)ᵀ s  =  Fᵀ J s
```

(usando que `(AB)ᵀ = BᵀAᵀ` y que `(Fᵀ)ᵀ = F`.)

Entonces la pregunta «¿el paso de Newton hace bajar el residuo?» se convierte en
«¿es `Fᵀ J s` negativo?». Eso es lo que vamos a calcular ahora.

### 6.5. La desigualdad que sostiene todo ⭐

Esta es la sección más importante del manual. Si algo hay que entender de verdad, es
esto.

**Punto de partida.** Tenemos un paso `s` de Newton **inexacto** (§4.3), o sea que no
cumple `J s = −F` exacto sino que le sobra un residuo `r`:

```
    J s  =  −F + r        con     ‖r‖ ≤ η‖F‖
```

**La pregunta.** ¿Es `s` una dirección de descenso para `f = ½‖F‖²`?

**La cuenta**, paso por paso, justificando cada igualdad:

```
    ∇f(x)ᵀ s
      = Fᵀ J s                       ← por §6.3 y §6.4
      = Fᵀ (−F + r)                  ← porque J s = −F + r, que es lo que sabemos
      = −FᵀF + Fᵀr                   ← distribuyendo el producto escalar
      = −‖F‖² + Fᵀr                  ← porque FᵀF = ‖F‖² (definición de norma, §2.1)
      ≤ −‖F‖² + ‖F‖·‖r‖              ← por Cauchy-Schwarz: Fᵀr ≤ |Fᵀr| ≤ ‖F‖‖r‖
      ≤ −‖F‖² + ‖F‖·(η‖F‖)           ← porque ‖r‖ ≤ η‖F‖, que es lo que exigimos
      = −‖F‖² + η‖F‖²
      = −(1 − η)‖F‖²
```

**El resultado:**

```
    ┌──────────────────────────────────────┐
    │    ∇f(x)ᵀ s   ≤   −(1 − η)·‖F‖²      │
    └──────────────────────────────────────┘
```

**Qué dice, en castellano.** Si `η < 1` y `F ≠ 0`, entonces `(1 − η) > 0` y `‖F‖² > 0`,
así que el lado derecho es **estrictamente negativo**. Por lo tanto `∇fᵀs < 0`:

> **El paso de Newton inexacto SIEMPRE apunta cuesta abajo para la función de mérito,
> por mal que esté calculado, con tal de que η < 1.**

Y como es dirección de descenso, por §2.3 **siempre existe algún `λ > 0` que hace
bajar `f`**. La búsqueda de línea nunca se queda sin opciones: siempre hay un paso
corto que sirve. Solo hay que encontrarlo.

**Y de yapa, explica por qué se exige `η < 1`.** Hasta ahora podía parecer una cota de
conveniencia. Ahora se ve que no:

- Si `η < 1`, la cota es negativa → hay descenso garantizado.
- Si `η ≥ 1`, la cota vale `≤ 0` y **no garantiza nada**: el paso podría ser de
  subida.

**`η < 1` es exactamente la condición que hace compatibles el Newton inexacto (§4) y
la globalización (§6).** Es la costura entre las dos mitades del tema, y es la
pregunta más probable del docente.

### 6.6. Estrategia A — búsqueda de línea con retroceso (Armijo)

**La idea.** Confiar en la *dirección* de Newton pero negociar la *distancia*:

```
    x_{k+1}  =  x_k  +  λ_k · s_k          con  λ ∈ (0, 1]
```

Ya la vimos funcionar a mano en §3.4 y en §3.6. Falta la regla formal de cuándo
aceptar un `λ`.

**La condición de Armijo**, escrita en términos del residuo:

```
    ‖F(x + λs)‖   ≤   (1 − α·λ) · ‖F(x)‖         con  α = 10⁻⁴
```

Desarmada:

- El lado izquierdo es el residuo **después** de dar el paso.
- El lado derecho es el residuo de **ahora**, multiplicado por un número apenas menor
  que 1.
- O sea: **exigimos que el residuo baje, aunque sea un poquito.**

**¿Por qué `α` tan chico?** Porque si exigiéramos mucha mejora casi ningún `λ` la
cumpliría y el método se arrastraría. `α = 10⁻⁴` pide una mejora casi simbólica: lo
justo para que el teorema de convergencia funcione (impide que el método se estanque
bajando cantidades cada vez más pequeñas) sin frenar nada. Se deja en `10⁻⁴` y no se
toca.

Medido: en el ejercicio de clase, con `α = 10⁻⁴` el problema converge en 9
iteraciones; con `α = 0.9` no converge en 40.

**El algoritmo (backtracking):**

```
    λ = 1
    mientras  ‖F(x + λs)‖ > (1 − αλ)‖F(x)‖ :
        λ = λ / 2                 (o mejor: interpolación cuadrática)
        si λ < λ_mín : abandonar
    aceptar x + λs
```

Se prueba `λ = 1` **primero**, y eso es importante: cerca de la solución el paso
completo cumple la condición, se acepta, y **se recupera la convergencia cuadrática
íntegra** (§2.7). La búsqueda de línea solo actúa cuando hace falta.

En la implementación real no se parte a la mitad a ciegas sino que se ajusta una
parábola a la información disponible y se salta a su mínimo, acotado al intervalo
seguro `[0.1λ, 0.5λ]`. Es más rápido y no cambia la idea.

**Ventajas y límites:**

| A favor | En contra |
|---|---|
| Trivial de implementar (~15 líneas) | Solo se mueve **sobre la recta de Newton**: si esa dirección es mala, no hay nada que hacer |
| Cada retroceso cuesta **una evaluación de F**, no un sistema lineal nuevo | Si `J` está mal condicionada, `s` es enorme y `λ` se va a cero |
| Es el valor por defecto de `scipy.optimize.newton_krylov` | Baja monótonamente… incluso hacia un mínimo local que no es raíz (§6.9) |

### 6.7. Estrategia B — región de confianza

**La idea.** Invertir la pregunta. En vez de «¿cuánto de este paso acepto?»,
preguntar **«¿hasta dónde le creo al modelo lineal?»**.

Se declara un radio `Δ` dentro del cual uno confía en la aproximación lineal, y se
busca el mejor paso **dentro de ese radio**:

```
    minimizar   m(s) = ½‖F + J s‖²          ← el modelo lineal
    sujeto a    ‖s‖ ≤ Δ                     ← pero sin alejarse más de Δ
```

**Cómo se ajusta `Δ`.** Se compara la mejora que el modelo **predijo** con la que
**realmente** ocurrió:

```
              reducción real          f(x) − f(x + s)
    ρ  =  ──────────────────────  =  ──────────────────
            reducción predicha         f(x) − m(s)

    ρ ≈ 1   → el modelo acertó       → agrandar Δ (confiar más)
    ρ < ¼   → el modelo mintió       → achicar Δ y rechazar el paso
```

**La ventaja real sobre la búsqueda de línea.** Como `Δ` acota explícitamente el
tamaño del paso, **la región de confianza sigue funcionando aunque `J` sea singular**.
En ese caso el paso de Newton es infinito o disparatado, y la búsqueda de línea se
rinde con `λ → 0`; la región de confianza simplemente lo recorta al borde y sigue.

**Cómo se resuelve el subproblema.** Hay dos maneras, y la elección **no** es
indiferente en el contexto matriz-libre. Este fue un hallazgo real de este trabajo:

*Opción 1 — Steihaug-CG.* Aplica gradientes conjugados a las **ecuaciones normales**
`JᵀJ s = −JᵀF`. Problema: como vimos en §2.5, eso **eleva al cuadrado el número de
condición**. Con un Jacobiano conocido solo con 8 cifras (§5.4), el resultado es
inutilizable: **en nuestras pruebas sobre Bratu se estancó alrededor de `10⁻⁶` y dejó
de progresar**. No es viable acá.

*Opción 2 — Dogleg de Powell.* Interpola entre dos direcciones que ya tenemos:

```
              s_C = punto de Cauchy                  s_N = paso de Newton
              (mínimo del modelo a lo largo          (el que ya calculó GMRES)
               del máximo descenso)

                    ⎧ s_N                       si ‖s_N‖ ≤ Δ   (cabe: tomarlo entero)
              s  =  ⎨ Δ·s_C/‖s_C‖               si ‖s_C‖ ≥ Δ   (ni el Cauchy cabe)
                    ⎩ s_C + τ(s_N − s_C)        en otro caso   (el "codo del perro")

         x ●───────→● s_C
                     ╲                     ← el camino tiene forma de pata doblada:
                      ╲                       de ahí el nombre "dogleg"
                       ╲
                        ● s_N
              ╭┈┈┈┈┈┈┈┈┈┈┈┈╮
              ╎  radio Δ   ╎   el paso se corta donde el camino cruza el borde
              ╰┈┈┈┈┈┈┈┈┈┈┈┈╯
```

El punto de Cauchy se calcula con `s_C = −(‖g‖²/‖Jg‖²)·g` donde `g = JᵀF = ∇f`
(§6.3). Cuesta **dos productos matriz-vector extra** y es lo que se usa en la
práctica (Pawlowski et al., 2006). Es la opción que implementamos.

> **Una limitación honesta del enfoque matriz-libre.** El dogleg necesita `g = JᵀF`,
> o sea multiplicar por `Jᵀ`. Y la fórmula de diferencias finitas de §5.4 da `J·v`,
> **no** `Jᵀ·w`: no hay manera de obtener el traspuesto con una evaluación de `F`.
> Se sale del paso aprovechando que `J` es simétrica cuando lo es —cierto de forma
> exacta en operadores de difusión-reacción con diferencias centradas— o formando `J`
> si `n` es chico. Es una limitación real, no un detalle de implementación.

### 6.8. Estrategia C — continuación pseudo-transitoria (Ψtc)

**La idea, que es distinta y más física.** El estado que buscamos es un estado
**estacionario**: el punto donde el sistema se queda quieto. Un estado estacionario
es el límite de un transitorio. Entonces: en vez de saltar directo a la respuesta,
**simulemos la evolución en el tiempo hasta que se estabilice**.

Se introduce un tiempo artificial `t` y se integra:

```
    dx/dt  =  −F(x) ,        x(0) = x₀
```

Cuando `dx/dt = 0`, o sea cuando el sistema se detiene, es porque `F(x) = 0`: nuestra
solución.

Integrando con **Euler implícito** con paso de tiempo `δ`, cada paso queda:

```
    (I + δ·J) · s  =  −δ · F
```

**Los dos extremos son muy ilustrativos:**

```
    δ → 0     ⟹   s ≈ −δ·F        máximo descenso muy amortiguado.
                                   Robustísimo, lentísimo.

    δ → ∞     ⟹   J·s = −F        ¡Newton puro!
                                   Rápido, frágil.
```

O sea que `δ` es una perilla que va **continuamente** de «lento y seguro» a «rápido y
peligroso». Lo que hace Ψtc es girar esa perilla automáticamente con la **regla SER**
(*Switched Evolution Relaxation*):

```
    δ_{k+1}  =  δ_k · ‖F_k‖ / ‖F_{k+1}‖
```

Si el residuo bajó (`‖F_{k+1}‖ < ‖F_k‖`), la fracción es mayor que 1 y `δ` **crece**:
se acelera. Si no bajó, `δ` se achica. El método empieza amortiguado y **termina
siendo Newton puro**, recuperando la convergencia cuadrática al final.

> **Analogía.** En vez de saltar directo al estado final, simula la película desde el
> principio. Al comienzo va en cámara lenta; cuando ve que va bien, acelera; y al
> final ya está corriendo a toda velocidad.

**Dos advertencias, ambas descubiertas midiendo en este trabajo:**

1. **El signo importa.** El método exige que el flujo `dx/dt = −F(x)` sea
   **estable**, lo que obliga a elegir el signo del residuo de modo que `F′` quede
   definida positiva. Con el signo contrario, Ψtc diverge sistemáticamente. Por eso
   en el código el residuo de Bratu se escribe con un menos delante.

2. **`δ₀` es un parámetro de verdad**, no un detalle. Medido sobre Burgers:

   | δ₀ | Resultado | Productos J·v |
   |---|---|---|
   | 10⁻² | no converge | 2 499 |
   | 1 | no converge | 2 772 |
   | 70.7 = 1/‖F₀‖ (por defecto) | converge | 863 |
   | 10⁴ | converge | 200 |
   | 10⁶ | converge | 315 |
   | 10⁹ | no converge | 15 586 |

   Falla por defecto (sobre-amortigua y se arrastra) **y** por exceso (equivale a
   Newton puro y se pierde la robustez). El valor por defecto razonable es
   `δ₀ = 1/‖F(x₀)‖`.

**Cuándo brilla.** Problemas dominados por convección y estados estacionarios de
EDPs: dinámica de fluidos, transporte, combustión. Es el estándar de facto para
Navier-Stokes estacionario. Y en nuestro experimento con Burgers fue, por lejos, la
mejor (→ §9.3).

### 6.9. Qué NO garantiza la globalización

Acá está el límite del tema, y decirlo bien da más crédito que vender el método como
infalible.

**El enunciado exacto del teorema es:**

> Si `F` es continuamente diferenciable y los iterados permanecen en un conjunto
> acotado, la globalización garantiza que la sucesión converge a un **punto
> estacionario de `f = ½‖F‖²`**.

Un **punto estacionario** de `f` es un punto donde `∇f = 0`. Y `∇f = JᵀF` (§6.3), así
que `∇f = 0` puede ocurrir de dos maneras:

```
    ① F = 0                    ← ¡bien! es una raíz, es lo que buscábamos
    ② F ≠ 0  pero  JᵀF = 0     ← mal: J es singular ahí. No es raíz de nada.
```

El caso ② es real y se manifiesta de dos formas:

**Modo de fallo 1 — mínimo local de la función de mérito.** La función de
Freudenstein y Roth (un clásico de la literatura) tiene un mínimo local en
`x₂ ≈ −0.8968` con `f ≈ 49 ≠ 0`. Una sucesión que solo sabe bajar cae en ese valle y
se queda ahí para siempre, sin haber resuelto nada.

**Modo de fallo 2 — Jacobiano singular.** En nuestro sistema 2×2, el Jacobiano es

```
         ⎡  2x₁      2x₂  ⎤
    J =  ⎢                ⎥      →   sobre la recta x₂ = 0 la segunda columna
         ⎣ e^(x₁−1)  3x₂² ⎦          se anula entera:  J es SINGULAR.
```

Si el descenso arrastra la iteración hacia esa recta, el paso de Newton se vuelve
enorme y mal condicionado, `λ` tiende a cero y la búsqueda de línea se rinde.

**Y la consecuencia incómoda: a veces Newton puro gana.** Medido sobre
Freudenstein-Roth, con 1681 puntos iniciales:

```
    Newton puro     → llega a la raíz desde  82.9 %  de los puntos
    Newton + Armijo → llega a la raíz desde  36.6 %  de los puntos
```

¿Cómo puede ser? Porque los pasos salvajes de Newton puro **saltan por encima** del
valle donde el método prudente queda atrapado. La globalización impide empeorar, y a
veces empeorar temporalmente era justo lo que hacía falta.

**Esto no invalida el método.** En el mismo experimento, sobre un problema con
exponencial empinada, la búsqueda de línea gana (100 % contra 98 %), y en el problema
grande de EDP la diferencia es abrumadora a favor de globalizar (→ §9.3). Lo que hay
que retener es qué se promete exactamente:

> **La globalización promete que `‖F‖` nunca aumenta y que se converge a un punto
> estacionario de `f`. No promete encontrar la raíz.**
>
> Saber exactamente qué se promete es lo que permite diagnosticar cuando falla, en
> vez de culpar al programa.

### 6.10. Cuál usar

| | Búsqueda de línea | Región de confianza | Ψtc |
|---|---|---|---|
| **Qué ajusta** | longitud `λ` del paso | radio `Δ` de confianza | paso de tiempo `δ` |
| **Dirección** | siempre la de Newton | puede cambiarla (dogleg) | interpola Newton ↔ descenso |
| **Costo extra por paso** | 1 evaluación de `F` por retroceso | 2 productos `J·v` + posible rechazo | ninguno |
| **Si `J` es singular** | `λ → 0`, se atasca | `Δ` lo acota, continúa | `I + δJ` es regular |
| **Implementación** | ~15 líneas | ~60 líneas | ~20 líneas |
| **Parámetros que tocar** | `α` (nunca se toca) | `Δ₀`, umbrales de `ρ` | **`δ₀` (sí importa)** |
| **Cuándo usarla** | por defecto, siempre | `J` mal condicionada | estados estacionarios de EDPs |

Y no son excluyentes: en producción se combinan. Lo habitual es hacer **continuación
en un parámetro físico** (ir subiendo `λ`, o el número de Reynolds, de a poquito,
usando cada solución como punto inicial de la siguiente) para generar un buen `x₀`, y
después aplicar Newton-Krylov con búsqueda de línea desde ahí.

---

## 7. El método completo, armado

Ya tenemos todas las piezas. Así queda el método entero:

```python
# ══ NEWTON-KRYLOV GLOBALIZADO ══════════════════════════════════════════
x = x0
F = evaluar(x)                       # el residuo en el punto actual

mientras ‖F‖ > tolerancia:

    # ─── 1. ¿CUÁN BIEN resolver el sistema lineal?  ────────  → §4.6
    η = 0.9 · (‖F‖ / ‖F_anterior‖)^1.618          # Eisenstat-Walker
    η = min(0.9, η)                               # nunca más flojo que 0.9

    # ─── 2. El Jacobiano como FUNCIÓN, no como matriz  ─────  → §5.4
    Jv = lambda v:  (evaluar(x + ε·v) − F) / ε     # una evaluación de F

    # ─── 3. Resolver SOLO HASTA η  ─────────────────────────  → §4.3, §5.3
    s = gmres(Jv, −F, tolerancia_relativa=η, precondicionador=M)
    #   ahora se cumple  ‖J·s + F‖ ≤ η‖F‖   pero NO  J·s = −F

    # ─── 4. GLOBALIZAR: ¿cuánto de ese paso acepto?  ───────  → §6  ← EL TEMA
    λ = buscar_paso(x, s, F)         # Armijo, o región de confianza, o Ψtc

    # ─── 5. Avanzar  ───────────────────────────────────────
    x = x + λ·s
    F_anterior, F = F, evaluar(x)
# ═══════════════════════════════════════════════════════════════════════
```

Tres observaciones sobre este pseudocódigo:

1. **Lo único que el usuario tiene que escribir es la función `F(x)`.** Todo lo demás
   es maquinaria genérica. Opcionalmente puede aportar un precondicionador, que es
   donde se gana o se pierde el 90 % del rendimiento (§5.6).

2. **Las líneas 1 a 3 hacen que el método sea *posible*** en problemas grandes.
   **La línea 4 lo hace *confiable*.** Sin ella, todo lo anterior sirve de poco: se
   puede resolver eficientísimamente un sistema lineal que lo lleva a uno a ninguna
   parte.

3. **Nunca aparece la matriz `J`.** Ni en la línea 2, ni en la 3, ni en ningún lado.

---

## 8. El código, por dentro

Todo el método está implementado desde cero en este proyecto, con NumPy y SciPy como
únicas dependencias. Esta sección lo recorre.

### 8.1. Instalar y correr

```bash
pip install numpy scipy matplotlib
```

```bash
# Demostraciones de la exposición (generan las figuras)
cd ejercicios_exposicion
python3 ej1_newton_vs_globalizado.py        # ~20 s
python3 ej2_bratu1d_newton_krylov.py        # ~20 s
python3 ej3_comparativa_globalizacion.py    # ~40 s
cd ..

# Ejercicios para los compañeros
cd ejercicios_clase
python3 clase1_armijo.py                    # ~5 s
python3 clase2_matrix_free.py               # ~25 s
python3 clase3_bratu_forcing.py             # ~35 s
```

No hace falta `python-docx`, ni LaTeX, ni conexión a internet.

### 8.2. `nk_lib.py`, función por función

Este archivo es el núcleo compartido. Cada función corresponde a una sección de este
manual.

---

**`class ContadorF`** — el que mide el costo honestamente.

```python
class ContadorF:
    def __init__(self, F):
        self.F = F
        self.n_F = 0        # evaluaciones totales de F
        self.n_Jv = 0       # de esas, cuántas vinieron de un producto J·v

    def __call__(self, x):
        self.n_F += 1
        return self.F(np.asarray(x, dtype=float))
```

Envuelve la función `F` del usuario y cuenta cuántas veces se la llama. Existe porque
—como argumentamos en §4.5— contar iteraciones de Newton engaña: lo que se paga son
evaluaciones de `F`. Todas las tablas de costo de este trabajo salen de este contador.

---

**`jv_diferencias_finitas(F, x, Fx, v)`** — el corazón del método (**→ §5.4, §5.5**).

```python
    xs = float(np.dot(x, v)) / norma_v      # proyección de x sobre v
    eps = SQRT_EPS
    if xs != 0.0:
        eps = eps * max(abs(xs), 1.0) * math.copysign(1.0, xs)
    eps = eps / norma_v

    return (F(x + eps * v) - Fx) / eps
```

Es la fórmula `[F(x+εv) − F(x)]/ε` con la elección de `ε` de Kelley (2003). Tres
detalles que valen la pena:

- **`max(abs(xs), 1.0)`** escala `ε` con la magnitud de `x`. Si `x` tuviera
  componentes del orden de `10⁶`, un `ε` fijo de `10⁻⁸` sería relativamente absurdo.
- **`/ norma_v`** normaliza por el tamaño de `v`, para que el desplazamiento real sea
  el que queremos.
- **`copysign`** preserva el signo, para no restar cantidades que se cancelen (§2.4).

Cuesta **una** evaluación de `F` (ya teníamos `Fx` guardado).

---

**`operador_jacobiano(F, x, Fx, contador)`** — el envoltorio para SciPy.

Devuelve un `LinearOperator` de SciPy, que es un objeto que *se comporta* como una
matriz para quien lo multiplica, pero por dentro solo tiene la función de arriba.
Esto es lo que permite pasárselo a `gmres` sin que GMRES se entere de que la matriz
no existe.

---

**`forzado_eisenstat_walker(...)`** — la elección automática de `η` (**→ §4.6**).

```python
    elif tipo == "ew2":
        eta = GAMMA_EW * (normaF / normaF_prev) ** ALPHA_EW

    # Salvaguarda de Eisenstat-Walker: no bajar η más rápido de lo sensato.
    if eta_prev is not None:
        piso = GAMMA_EW * eta_prev ** ALPHA_EW
        if piso > 0.1:
            eta = max(eta, piso)
```

(`GAMMA_EW = 0.9` y `ALPHA_EW = (1+√5)/2 ≈ 1.618` son constantes del módulo.)

Implementa las dos variantes (*Choice 1* y *Choice 2*) más la salvaguarda descrita en
§4.6 y un tope final para no pedir más precisión de la que la tolerancia justifica.

---

**`linea_armijo(F, x, s, normaF, alpha=1e-4, ...)`** — la estrategia A (**→ §6.6**).

```python
    while not (n_t <= (1.0 - alpha * lam) * normaF):
        ...
        lam_nuevo = min(max(lam_nuevo, sigma0 * lam), sigma1 * lam)
```

Dos cosas para notar:

- La condición está escrita como `not (n_t <= cota)` y **no** como `n_t > cota`.
  Parecen lo mismo pero no lo son: si el paso completo desborda, `n_t` vale `nan`, y
  en punto flotante `nan > cota` es `False` — o sea que la segunda versión
  **aceptaría un paso inválido**. Esto era un bug real que apareció al probar el
  código y se corrigió.
- `sigma0 = 0.1`, `sigma1 = 0.5` son las salvaguardas del retroceso: el nuevo `λ`
  siempre queda entre el 10 % y el 50 % del anterior.

---

**`paso_dogleg(Jadj, Fx, s_newton, delta)`** — la estrategia B (**→ §6.7**).

Recibe el paso de Newton que GMRES ya calculó y decide cuánto de él cabe en la región
de confianza, interpolando con el punto de Cauchy si hace falta. Tres casos, tal como
en el dibujo de §6.7.

---

**`steihaug_cg(Jop, Fx, delta, ...)`** — la alternativa que **no** usamos.

Está implementada y documentada a propósito: es la manera «de libro» de resolver el
subproblema de la región de confianza, y sirve para mostrar por qué no funciona en el
contexto matriz-libre (§6.7). Es útil en problemas chicos donde el Jacobiano se puede
formar.

---

**`newton_krylov(F, x0, globalizacion=..., forzado=..., ...)`** — el bucle principal.

```python
def newton_krylov(F, x0, globalizacion="linea", forzado="ew2",
                  tol=1e-10, max_iter=60, precond=None,
                  gmres_restart=30, gmres_maxiter=200,
                  delta0=1.0, delta_max=1e3, ptc_delta0="auto",
                  jacobiano_simetrico=False, max_evals_F=None,
                  guardar_trayectoria=False, verbose=False):
```

Los dos parámetros que importan para el tema:

| Parámetro | Valores | Qué hace |
|---|---|---|
| `globalizacion` | `"ninguna"` | Newton inexacto puro, paso completo. Sirve para **mostrar que falla**. |
| | `"linea"` | búsqueda de línea de Armijo (§6.6) |
| | `"region"` | región de confianza con dogleg (§6.7) |
| | `"ptc"` | continuación pseudo-transitoria (§6.8) |
| `forzado` | un número | `η` fijo en todas las iteraciones |
| | `"ew1"` / `"ew2"` | Eisenstat-Walker adaptativo (§4.6) |

Devuelve un objeto `Historial` con todo lo que hace falta para graficar y explicar:

```python
    residuales    # ‖F‖ en cada iteración
    etas          # el η que se eligió en cada paso
    lambdas       # el λ que aceptó la búsqueda de línea
    krylov_iters  # iteraciones de GMRES por paso
    n_F, n_Jv     # los contadores de trabajo
    convergio     # True / False
    motivo        # por qué se detuvo
```

### 8.3. Los tres ejercicios de la exposición

| Archivo | Qué demuestra |
|---|---|
| `ej1_newton_vs_globalizado.py` | **Por qué hace falta globalizar.** arctan en 1D, el sistema 2×2, mapas de cuencas de convergencia, y los casos donde la globalización también falla. |
| `ej2_bratu1d_newton_krylov.py` | **El método completo.** Verifica `J·v` contra el Jacobiano analítico, barre `ε`, compara cinco maneras de elegir `η`, mide el precondicionador y contrasta contra `scipy.optimize.newton_krylov`. |
| `ej3_comparativa_globalizacion.py` | **Las cuatro estrategias comparadas** sobre la ecuación de Burgers, desde 7 puntos iniciales distintos. |

### 8.4. Cómo leer la salida

Cuando corra los programas, esto es lo que hay que mirar:

**En las tablas.** Casi todas tienen una columna de **iteraciones de Newton** y otra
de **productos J·v**. La segunda es la que importa (§4.5). Si las dos apuntan en
direcciones opuestas, esa es exactamente la lección del ejercicio.

**En las gráficas de residuo.** El eje vertical es `‖F(x_k)‖` en escala logarítmica.
Dos cosas para mirar:

- **¿La curva sube alguna vez?** Si sube, no había globalización. Si es monótona
  decreciente, la globalización está haciendo su trabajo — y esa monotonía *es*
  literalmente la garantía del teorema (§6.9).
- **¿Cae en picada al final?** Eso es la convergencia cuadrática entrando en acción:
  el iterado ya está dentro de la bola del teorema (§3.5).

**En el eje horizontal.** Cuando dice «productos J·v acumulados» en vez de
«iteración», la gráfica está midiendo trabajo real. Dos métodos que se ven parecidos
por iteración pueden estar separados por un factor 10 por trabajo.

**En los mapas de cuencas.** Cada píxel es un punto inicial distinto; el color dice
cuántas iteraciones hicieron falta, y el rojo oscuro dice «no llegó a ninguna raíz».
Cuanta más superficie de color, más robusto el método.

### 8.5. Experimentos que puede hacer usted

El código está pensado para que se lo pueda toquetear. Algunas cosas que vale la pena
probar, en orden de dificultad:

1. **En `ej1`**, cambie `x0` en la Parte B. Busque un punto desde el cual Newton puro
   converja y Armijo no (existen: vea §6.9).
2. **En `ej2`**, cambie `LAMBDA = 3.0` a `3.5`. Está más cerca del punto de retorno
   `λ* ≈ 3.5138`, donde el problema deja de tener solución. ¿Crece el trabajo?
   ¿Cuál `η` aguanta mejor?
3. **En `ej2`**, cambie `N = 250` a `N = 1000`. Sin precondicionador el trabajo
   debería crecer mucho (porque `cond = O(n²)`, §2.5); con precondicionador, casi
   nada. Eso es lo que significa que un método «escale».
4. **En `ej3`**, baje `NU = 0.01` a `0.005`. La capa límite se afila y el problema se
   endurece. ¿Cuál estrategia aguanta más?
5. **En `clase3`**, comente la línea de la salvaguarda de Eisenstat-Walker y vea qué
   pasa con el total de productos `J·v`.

---

## 9. Los resultados, explicados

Todas las cifras de esta sección salen de correr el código de §8. Acá no solo se
listan: se explica **por qué** son las que son.

### 9.1. Por qué hace falta globalizar (ejercicio 1)

Desde `x₀ = (2.0, 0.5)` en el sistema 2×2:

| Estrategia | ¿Converge? | Iteraciones | ‖F‖ final | Evaluaciones de F |
|---|---|---|---|---|
| Newton puro | **no** | — | `inf` (desbordó) | 138 |
| Búsqueda de línea (Armijo) | sí | 7 | 9.9·10⁻¹² | 108 |
| Región de confianza (dogleg) | sí | 6 | 0 (exacto) | 799 |

**Por qué falla Newton puro:** ya lo vimos en §3.4. El determinante del Jacobiano en
`x₀` vale 0.28, casi cero, así que el paso se dispara a `x₂ ≈ 10` y a partir de ahí
el término `e^(x₁−1)` explota.

**Por qué la región de confianza gasta 799 evaluaciones y Armijo solo 108:** porque
el dogleg necesita `Jᵀ` (§6.7), que en dimensión 2 se obtiene formando el Jacobiano
columna a columna. En un problema chico eso es aceptable; en uno grande sería
impensable.

Y un matiz importante, que salió midiendo y no aparece en las presentaciones
habituales del tema:

> **En sistemas pequeños y suaves, Newton puro ya es bastante robusto.** Con un
> presupuesto de 12 iteraciones y 2601 puntos iniciales, Newton puro converge desde
> el **86.6 %** y Newton con Armijo desde el **95.3 %**. La diferencia es real pero
> moderada.

El abismo aparece cuando el sistema viene de discretizar una EDP. Eso es §9.3.

### 9.2. El costo de resolver de más y el valor de precondicionar (ejercicio 2)

Ya vimos las dos tablas en §4.5 y §5.6. Lo que hay que saber decir sobre ellas:

**Sobre el oversolving.** El punto no es que `η = 10⁻¹²` sea «malo»: converge
perfectamente, y en menos iteraciones de Newton que nadie. El punto es que
**cambiar pocas iteraciones caras por muchas baratas conviene**, y que la métrica
equivocada (iteraciones) da la respuesta al revés.

**Sobre el precondicionador.** La ganancia de 519× parece de otro mundo, y hay una
razón concreta: el laplaciano **es** casi todo el Jacobiano en Bratu, porque el
término no lineal `−h²λ·e^u` es una perturbación diagonal pequeña. Es un
precondicionador casi perfecto para *este* problema. En uno dominado por convección
no lo sería.

**Sobre la validación.** El ejercicio contrasta el resultado contra
`scipy.optimize.newton_krylov`, que es la implementación de referencia:

```
    ‖F‖ final, nuestra implementación : 2.8·10⁻¹¹
    ‖F‖ final, SciPy                  : 1.4·10⁻⁰⁹
    diferencia relativa entre soluciones: 8.6·10⁻⁰⁸
```

Coinciden en el orden de la precisión alcanzable con un Jacobiano por diferencias
finitas, que es exactamente lo que predice §5.4. Esto es lo que da confianza en que
la implementación propia está bien.

### 9.3. El resultado central: las cuatro estrategias sobre Burgers (ejercicio 3)

**El problema.** La ecuación de Burgers estacionaria:

```
    −ν u″ + u u′ = 0 ,     u(0) = 1 ,  u(1) = −1 ,   ν = 0.01
```

Su solución tiene una **capa límite interna**: una transición casi vertical en
`x = ½` de anchura del orden de `2ν = 0.02`. Es el prototipo de lo que aparece al
resolver flujo de fluidos estacionario.

```
     u(x)
    +1 ├──────────────────┐
       │                  │      ← la transición ocurre en un espesor 2ν
     0 ┼──────────────────┼──────
       │                  │
    −1 │                  └──────────────────
       └──────────────────────────────────── x
       0                 0.5                1
```

Un detalle metodológico que conviene mencionar en la exposición: el **número de
Péclet de celda** vale `h/ν = 0.5 < 2`, o sea que las diferencias centradas **no**
producen oscilaciones espurias. La dificultad del problema es genuinamente la no
linealidad, no un esquema numérico inestable. Si lo fuera, ninguna estrategia de
globalización lo arreglaría.

**El experimento.** Siete puntos iniciales `u₀ = a·(1−2x)` con `a` de −3 a +3.
Todas las corridas comparten el paso de Newton, GMRES, el precondicionador y el
término de forzado. **Lo único que cambia es la globalización.**

| Estrategia | Converge desde | Tasa | J·v (mediana) |
|---|---|---|---|
| Newton inexacto, sin globalizar | 1 de 7 | 14 % | 10 631 |
| Búsqueda de línea (Armijo) | 5 de 7 | 71 % | 3 114 |
| Región de confianza (dogleg) | 3 de 7 | 43 % | 3 174 |
| **Continuación pseudo-transitoria** | **7 de 7** | **100 %** | **863** |

Y el detalle caso por caso:

```
     a  │ sin glob. │ línea    │ región   │ Ψtc
    ────┼───────────┼──────────┼──────────┼──────────
    −3  │    ··     │ OK  88it │    ··    │ OK 108it
    −2  │    ··     │    ··    │    ··    │ OK  90it
    −1  │ OK  52it  │ OK  25it │    ··    │ OK  81it
    +0  │    ··     │    ··    │ OK  54it │ OK 141it
    +1  │    ··     │ OK  18it │ OK  12it │ OK  55it
    +2  │    ··     │ OK 149it │ OK  93it │ OK  86it
    +3  │    ··     │ OK 119it │    ··    │ OK 130it
```

**Tres lecturas de esta tabla:**

1. **El contraste es categórico.** Sin globalizar se resuelve uno de siete casos; con
   Ψtc, los siete.

2. **Ψtc es además la más barata**, por un factor de doce respecto a no globalizar.
   Esto suele sorprender, y la razón es que las corridas sin globalización consumen
   enormes cantidades de trabajo **divergiendo** antes de agotar su presupuesto. La
   robustez no se paga: acá se cobra.

3. **La búsqueda de línea le gana a la región de confianza**, contra lo que uno
   esperaría. En este problema la dirección de Newton es buena; lo que hacía falta era
   recortarla, no cambiarla.

### 9.4. El teorema, visto en una gráfica

Para `u₀ ≡ 0`, un punto inicial perfectamente razonable:

| Estrategia | ‖F‖ inicial | ‖F‖ **máximo** | ‖F‖ final | ¿Empeora? |
|---|---|---|---|---|
| Sin globalizar | 1.4·10⁻² | **3.9·10¹** | 2.8·10⁰ | **SÍ, sube** |
| Búsqueda de línea | 1.4·10⁻² | 1.4·10⁻² | 1.3·10⁻⁵ | nunca |
| Región de confianza | 1.4·10⁻² | 1.4·10⁻² | 6.1·10⁻¹⁰ | nunca |
| Ψtc | 1.4·10⁻² | 1.4·10⁻² | 5.2·10⁻¹¹ | nunca |

Mire la columna del máximo. Sin globalizar, el residuo llega a ser **2 700 veces
peor** que el punto de partida. Con cualquiera de las tres estrategias, la columna
«máximo» es idéntica a la «inicial»: **nunca subió, ni una vez.**

Eso no es suerte. Es literalmente lo que impone la condición de Armijo
`‖F(x+λs)‖ ≤ (1−αλ)‖F(x)‖` y su equivalente `ρ > 0` en la región de confianza. **La
monotonía es la garantía del teorema, escrita en el código.**

### 9.5. Un residuo chico no es una solución precisa

Este resultado apareció de casualidad al comparar contra la solución analítica de
Burgers, y es de lo más instructivo del trabajo:

```
    ‖F(u_numérica)‖              = 6.1·10⁻¹⁰      ← convergió perfectamente
    ‖F(u_exacta muestreada)‖     = 2.6·10⁻⁰⁵      ← ¡la solución exacta NO es
                                                    solución del sistema discreto!
    error máximo entre ambas     = 3.0·10⁻⁰¹      ← una diferencia enorme
    la capa numérica cruza cero en x = 0.4925 (la exacta, en 0.5000)
```

**¿Qué pasó?** No es un error del solver. La **posición** de la capa límite está
exponencialmente mal determinada: el Jacobiano tiene un autovalor del orden de
`e^(−1/ν)`, o sea que desplazar la capa entera casi no cambia el residuo. Bajar `‖F‖`
a `10⁻¹⁰` no fija la posición mejor que eso.

**La lección**, que ya adelantamos en §2.1:

> Lo que acota el error de la solución no es `‖F‖`, es `‖J⁻¹‖·‖F‖`. Si `J` está mal
> condicionada, `‖J⁻¹‖` es enorme y un residuo diminuto no garantiza nada.

---

## 10. Glosario

| Término | Qué es | Fórmula | Dónde |
|---|---|---|---|
| **Sistema no lineal** | `n` ecuaciones con `n` incógnitas donde alguna incógnita no aparece «sola» | `F(x) = 0` | §1.1 |
| **Residuo** | Cuánto le falta a `x` para ser solución | `F(x)` | §1.2 |
| **Discretizar** | Cambiar una función incógnita por sus valores en unos pocos puntos | — | §1.4 |
| **Norma** | Longitud de un vector; convierte `n` números en uno | `‖v‖ = √(Σvᵢ²)` | §2.1 |
| **Cauchy-Schwarz** | El producto escalar no supera al producto de las longitudes | `∣uᵀv∣ ≤ ‖u‖‖v‖` | §2.1 |
| **Jacobiano** | Matriz de todas las derivadas parciales de `F` | `Jᵢⱼ = ∂Fᵢ/∂xⱼ` | §2.2 |
| **Gradiente** | Vector de derivadas parciales de una función **escalar** | `∇f` | §2.3 |
| **Dirección de descenso** | Dirección en la que `f` baja, al menos un poco | `∇fᵀs < 0` | §2.3 |
| **Épsilon de máquina** | El menor `ε` con `1 + ε ≠ 1` en la computadora | `εₘ ≈ 2.2·10⁻¹⁶` | §2.4 |
| **Cancelación** | Pérdida de cifras al restar números parecidos | — | §2.4 |
| **Número de condición** | Cuánto amplifica una matriz los errores | `cond(A)` | §2.5 |
| **Método directo** | Da la respuesta exacta en un número fijo de pasos (LU, Gauss) | `O(n³)` | §2.6 |
| **Método iterativo** | Mejora una aproximación; se puede cortar cuando uno quiera | — | §2.6 |
| **Orden de convergencia** | Cuán rápido se achica el error por paso | lineal / superlineal / cuadrático | §2.7 |
| **Paso de Newton** | El vector `s` que resuelve el sistema lineal de Newton | `J s = −F` | §3.3 |
| **Convergencia local** | Funciona solo si el punto inicial ya estaba cerca | — | §3.5 |
| **Newton inexacto** | Resolver el sistema lineal solo aproximadamente | `‖Js + F‖ ≤ η‖F‖` | §4.3 |
| **Término de forzado** | El `η` de arriba: cuánta pereza se permite | `η ∈ [0,1)` | §4.3 |
| **Oversolving** | Resolver el sistema lineal con más precisión de la que el modelo merece | — | §4.5 |
| **Eisenstat-Walker** | Regla que elige `η` sola mirando cómo bajó el residuo | `η_k = 0.9(‖F_k‖/‖F_{k−1}‖)^1.618` | §4.6 |
| **`span`** | Todas las combinaciones lineales de unos vectores | — | §5.1 |
| **Subespacio de Krylov** | El subespacio generado por las potencias de `A` sobre `b` | `𝒦ₘ = span{b, Ab, …}` | §5.2 |
| **GMRES** | Método de Krylov que minimiza el residuo; no exige simetría | — | §5.3 |
| **JFNK / matriz-libre** | Newton-Krylov con `J·v` por diferencias finitas; `J` nunca se forma | `J·v ≈ [F(x+εv)−F(x)]/ε` | §5.4 |
| **Precondicionador** | Matriz `M` fácil de invertir que agrupa el espectro de `M⁻¹J` | `M⁻¹Js = −M⁻¹F` | §5.6 |
| **Convergencia global** | Converger desde cualquier `x₀`. **No** es «mínimo global» | — | §6.1 |
| **Función de mérito** | El escalar que mide el progreso | `f = ½‖F‖²` | §6.2 |
| **Condición de Armijo** | Regla que decide si un paso `λ` se acepta | `‖F(x+λs)‖ ≤ (1−αλ)‖F(x)‖` | §6.6 |
| **Backtracking** | Probar `λ = 1`, después `½`, `¼`… hasta que se cumpla Armijo | — | §6.6 |
| **Región de confianza** | Radio `Δ` dentro del cual uno le cree al modelo lineal | `min m(s)` con `‖s‖ ≤ Δ` | §6.7 |
| **Punto de Cauchy** | Mínimo del modelo a lo largo del máximo descenso | `s_C = −(‖g‖²/‖Jg‖²)g` | §6.7 |
| **Dogleg** | Camino quebrado entre el punto de Cauchy y el paso de Newton | — | §6.7 |
| **Steihaug-CG** | CG truncado sobre las ecuaciones normales. Inviable matriz-libre | `JᵀJs = −JᵀF` | §6.7 |
| **Ψtc** | Continuación pseudo-transitoria: simular el transitorio | `(I + δJ)s = −δF` | §6.8 |
| **Regla SER** | La regla que alarga `δ` cuando el residuo baja | `δ_{k+1} = δ_k‖F_k‖/‖F_{k+1}‖` | §6.8 |
| **Punto estacionario** | Donde `∇f = 0`. Puede no ser raíz de `F` | `JᵀF = 0` | §6.9 |
| **Punto de retorno** | Valor del parámetro donde dos ramas de soluciones se juntan y `J` se vuelve singular | `λ* ≈ 3.5138` en Bratu | §8.5 |

---

## 11. Preguntas de quien recién aprende

Estas son distintas de las que va a hacer el docente (esas están en §16). Son las que
uno se hace estudiando.

**¿Por qué no uso simplemente `scipy.optimize.newton_krylov` y me olvido?**
En un trabajo real, eso es exactamente lo que hay que hacer, y así lo dice el informe.
Implementarlo desde cero sirve para dos cosas: entenderlo, y poder diagnosticar
cuando la biblioteca falle. Además SciPy no trae continuación pseudo-transitoria, que
en nuestro experimento fue la única estrategia que resolvió los siete casos.

**¿Por qué no calculo el Jacobiano a mano y lo programo?**
Se puede, y si el problema es chico o el Jacobiano es fácil, conviene: es más preciso
que las diferencias finitas y no tiene el piso de 8 cifras de §5.4. El problema es
que en un código de simulación real `F` puede tener miles de líneas y derivar eso a
mano es impracticable y propenso a errores. La diferenciación automática es una
alternativa moderna, pero no siempre está disponible.

**¿Qué pasa si `F` no es diferenciable?**
Entonces Newton no aplica: todo el método se apoya en la aproximación lineal.
Hay métodos para ese caso (semismooth Newton, métodos sin derivadas tipo Nelder-Mead
o DFO), pero son otro tema. Si `F` es diferenciable a trozos y uno cae justo en un
quiebre, en la práctica se sale del paso perturbando un poco el punto.

**¿Cómo sé si convergí de verdad?**
Buena pregunta, y §9.5 muestra que no es obvia. Criterios habituales:
`‖F(x)‖ ≤ tol_abs + tol_rel·‖F(x₀)‖` (residuo), y `‖x_{k+1} − x_k‖` pequeño
(el iterado dejó de moverse). Si los dos se cumplen y `J` no está patológicamente mal
condicionada, se puede confiar. Si `‖F‖` es chico pero el iterado sigue moviéndose,
sospeche.

**¿Esto sirve para optimización, o solo para resolver `F(x) = 0`?**
Sirve, y de hecho la relación es íntima. Minimizar `g(x)` equivale a resolver
`∇g(x) = 0`, que es un sistema no lineal con `F = ∇g` y `J = ∇²g` (la Hessiana).
Los métodos de Newton-Krylov aplicados a eso se llaman *Newton truncado* y son
estándar en optimización a gran escala. La búsqueda de línea y la región de confianza
que vimos vienen justamente de ahí.

**Si el paso de Newton siempre es de descenso (§6.5), ¿por qué a veces la búsqueda de
línea falla?**
Porque «es de descenso» significa que existe **algún** `λ > 0` que baja `f`, pero ese
`λ` puede ser absurdamente pequeño. Si `λ` cae por debajo de `10⁻¹⁰`, seguir no tiene
sentido numérico: se está haciendo un paso más chico que el ruido de redondeo. Ahí el
código se rinde. Eso pasa típicamente cuando `J` está casi singular.

**¿Por qué se mide en «productos J·v» y no en segundos?**
Porque los segundos dependen de la máquina, del lenguaje y de qué más esté corriendo.
Los productos `J·v` son una medida de trabajo **reproducible**, y en un método
matriz-libre son proporcionales al costo real, porque cada uno cuesta una evaluación
de `F`.

**¿Qué es un «operador» y por qué SciPy usa esa palabra?**
Un `LinearOperator` es un objeto que se comporta como matriz para quien lo multiplica
pero que por dentro puede ser cualquier cosa — en nuestro caso, una función que hace
una diferencia finita. Es la abstracción que permite que GMRES no se entere de que la
matriz no existe.

**Todo esto es de los años 80 y 90. ¿Sigue vigente?**
Sí. Es lo que hay dentro de PETSc y SUNDIALS, que son las bibliotecas que mueven
buena parte de la simulación científica del mundo hoy. Lo que cambió desde entonces
son los precondicionadores (multigrid algebraico, descomposición de dominios) y el
paralelismo, no el esqueleto del método.

**¿Cuál de las tres estrategias uso si tengo que elegir una y ya?**
Búsqueda de línea. Es la más simple, la más barata por paso, la que menos parámetros
tiene y el valor por defecto de las bibliotecas serias. Si le falla, y su problema es
un estado estacionario de una EDP, pruebe Ψtc.

---

# PARTE II — Guion de exposición

> Esta parte supone que ya leyó la PARTE I, o al menos §0.3. Es la chuleta para el
> día de la exposición. Cada bloque lleva una marca **(→ §X)** que remite a donde eso
> se explica en serio.

## 12. Minutaje y reparto

| Bloque | Diapositivas | Minutos | Quién | Manual |
|---|---|---|---|---|
| Apertura y planteo del problema | 1–3 | 2:30 | Maximiliano | §1 |
| Newton y su límite + demo arctan | 4–5 | 3:00 | Maximiliano | §3 |
| Newton inexacto y término de forzado | 6–8 | 4:00 | Maximiliano | §4 |
| Krylov y el Jacobiano sin Jacobiano | 9–12 | 4:30 | Iver | §5 |
| **Convergencia global: mérito y descenso** | 13–14 | 3:00 | Iver | **§6.2–§6.5** |
| Las tres estrategias | 15–18 | 4:00 | Iver | §6.6–§6.10 |
| Resultados y cierre | 19–24 | 4:00 | Maximiliano | §9 |
| **Total** | | **~22 min** | | |

**Si van atrasados**, se puede recortar en este orden:

1. La diapositiva 23 (errores comunes) se salta entera.
2. La 17 (Ψtc) se resume en dos frases: «integra el transitorio; empieza lento y
   acelera».
3. La 11 (pseudocódigo) se muestra 15 segundos sin recorrerla.

**Lo que nunca se salta:** la diapositiva **14**. Es el corazón del tema y es lo que
va a preguntar el docente.

## 13. Guion diapositiva por diapositiva

### 1 · Portada — 20 s
Presentarse, nombrar el tema, y la frase que orienta todo:
> «El método de Newton es rapidísimo cuando funciona. Hoy hablamos de qué se hace
> para que funcione siempre.»

### 2 · El problema — 50 s  **(→ §1)**
Plantear `F(x) = 0` y decir de dónde sale: discretizar una EDP no lineal, el paso
implícito de un integrador rígido, condiciones de optimalidad.
**Lo que hay que remarcar:** `n` no es 2 ni 10, es 10⁶. Una malla de 100×100×100 da
un millón de incógnitas y un Jacobiano de 8 terabytes.

### 3 · Agenda — 40 s  **(→ §3.7)**
Los tres dolores. Decirlo como hoja de ruta, no leer las cajas:
> «Newton tiene tres problemas. Los dos primeros se resuelven con Krylov; el tercero
> es el tema de hoy.»

### 4 · Newton y su teorema — 1:30  **(→ §3.3, §3.5)**
Escribir la iteración. Enunciar la convergencia cuadrática, y **detenerse en la letra
chica**: «existe δ > 0 tal que».
> «El teorema garantiza que existe un radio de convergencia, pero no dice cuánto vale
> ni da forma de calcularlo. Fuera de esa bola no promete nada, y "nada" incluye
> divergir.»

### 5 · Demo arctan — 1:30  **(→ §3.6)**
Es la diapositiva que convence. Señalar la figura izquierda con el dedo:
> «Miren: la tangente en x₀ = 2 cruza el eje en −3.5. La tangente ahí cruza en 13.9.
> Cada paso empeora. Y esto es arcotangente: una función suave, monótona, sin nada
> raro. El umbral es 1.3917.»

Después la derecha: mismo `x₀`, con búsqueda de línea converge en 5 iteraciones.

### 6 · Newton inexacto — 1:30  **(→ §4.2, §4.3, §4.4)**
La idea en una frase:
> «Si estamos lejos, el modelo lineal es una mentira. ¿Para qué resolver una mentira
> con dieciséis cifras?»

Escribir `‖J s + F‖ ≤ η‖F‖` y nombrar `η` como *término de forzado*.
Leer la tabla del teorema DES: `η` constante → lineal; `η → 0` → superlineal;
`η = O(‖F‖)` → cuadrática.

### 7 · Eisenstat-Walker — 1:15  **(→ §4.6)**
No hace falta que memoricen la fórmula; sí la **idea**:
> «Mira cuánto bajó el residuo en el paso anterior y pide para el siguiente
> exactamente esa precisión.»

Señalar la tabla de la derecha: `η` va de 0.9 a 0.012 sin que nadie se lo diga.

### 8 · Oversolving — 1:15  **(→ §4.5)**
**Diapositiva de impacto.** Mostrar la tabla y decir:
> «Fíjense en la fila resaltada. Resolver con doce cifras usa *la mitad* de
> iteraciones de Newton… y 2.6 veces más trabajo. Contar iteraciones externas
> engaña. Lo que se paga son evaluaciones de F.»

### 9 · Krylov — 1:15  **(→ §5.2, §5.3)**
Definir `𝒦ₘ(A,b) = span{b, Ab, …, A^(m−1)b}` y remarcar:
> «Construir esa base solo necesita multiplicar por A. Nunca hace falta la matriz.»

Justificar GMRES en vez de CG: **el Jacobiano no es simétrico** en cuanto hay
convección o transporte.

### 10 · J sin J — 1:15  **(→ §5.4, §5.5)**
Escribir `J·v ≈ [F(x+εv) − F(x)]/ε` y decir «una sola evaluación de F».
La elección de `ε`: dos errores opuestos, óptimo cerca de `√εₘ ≈ 1.5·10⁻⁸`.
Mencionar el precio: **el Jacobiano solo se conoce con ~8 cifras**.

### 11 · El algoritmo — 45 s  **(→ §7)**
Recorrer el pseudocódigo señalando dónde está cada dolor. Cerrar con:
> «Lo único que el usuario tiene que escribir es la función F.»

### 12 · Precondicionamiento — 1:15  **(→ §5.6)**
La ganancia de 519× es el número que se les va a quedar. Y **señalar la fila que no
mejora**:
> «η = 10⁻¹² no mejora ni con precondicionador. Le estamos pidiendo a GMRES una
> precisión que la doble precisión no puede dar. Pedir más precisión de la que existe
> no acelera: quema trabajo.»

### 13 · Función de mérito — 1:15  **(→ §6.2, §6.3)**
`f(x) = ½‖F(x)‖²`, `∇f = JᵀF`. Traducir raíces en minimización.
**Advertir de inmediato el precio:** toda raíz es mínimo de `f`, pero no todo mínimo
de `f` es raíz. Se retoma en la 22.

### 14 · La desigualdad — 1:45 ⭐ **LA DIAPOSITIVA IMPORTANTE**  **(→ §6.5)**
Hacer la cadena despacio, término por término:

```
∇f ᵀs = Fᵀ J s = Fᵀ(−F + r) = −‖F‖² + Fᵀr ≤ −‖F‖² + ‖F‖‖r‖ ≤ −(1−η)‖F‖²
```

Y rematar con las dos conclusiones:

1. El paso de Newton inexacto **siempre** apunta cuesta abajo si `η < 1`, así que
   siempre existe algún `λ` que sirve.
2. **Por eso se exige `η < 1`.** Si `η ≥ 1` la cota se vuelve `≤ 0` y no garantiza
   nada. Es la costura entre las dos mitades de la exposición.

### 15 · Búsqueda de línea — 1:15  **(→ §6.6)**
Condición de Armijo con `α = 10⁻⁴`. Por qué `α` tan chico: se pide una mejora casi
simbólica, lo justo para que el teorema funcione sin frenar el método.
Ventaja: cada retroceso cuesta **una evaluación de F**, no un sistema nuevo.

### 16 · Región de confianza — 1:15  **(→ §6.7)**
«Invierte la pregunta: no *cuánto acepto* sino *hasta dónde le creo al modelo*.»
Mencionar `ρ` y el ajuste de `Δ`. Y el detalle que demuestra que lo entendieron:
> «Steihaug-CG resuelve las ecuaciones normales, que elevan al cuadrado el
> condicionamiento. Con un Jacobiano por diferencias finitas se estanca en 10⁻⁶. Lo
> probamos. Por eso usamos dogleg.»

### 17 · Ψtc — 1:30  **(→ §6.8)**
La idea física: el estacionario es el límite de un transitorio, simulémoslo.
`(I + δJ)s = −δF`. Los dos extremos: `δ→0` es descenso amortiguado, `δ→∞` es Newton.
Regla SER. Y la honestidad sobre `δ₀`: la tabla muestra que falla por defecto y por
exceso.

### 18 · Comparativa — 45 s  **(→ §6.10)**
No leer la tabla entera. Señalar dos filas: **«si J es singular»** y **«parámetros a
ajustar»**. Cerrar diciendo que en producción se combinan.

### 19 · El problema de Burgers — 45 s  **(→ §9.3)**
Capa límite interna, anchura `≈ 2ν`. Y adelantarse a la objeción:
> «Péclet de celda 0.5, menor que 2: la discretización es estable. La dificultad es la
> no linealidad, no el esquema.»

### 20 · El resultado central — 1:15  **(→ §9.3)**
Dejar que la gráfica hable. Números: sin globalizar 1 de 7; Ψtc 7 de 7 **y doce veces
más barata**. Frase de cierre:
> «La robustez no se paga. Aquí se cobra.»

### 21 · El teorema en una gráfica — 1:00  **(→ §9.4)**
Señalar la curva roja subiendo. Es la demostración visual de que sin globalización no
hay ninguna garantía de mejora.

### 22 · Lo que no garantiza — 1:00  **(→ §6.9)**
Enunciado exacto: converge a un **punto estacionario de `f`**, que puede no ser raíz.
Los dos modos de fallo. Y el caso incómodo:
> «Sobre Freudenstein-Roth, Newton puro gana: 83 % contra 37 %. Sus pasos salvajes
> saltan el valle donde el método prudente se queda atrapado.»

Decirlo con naturalidad: mostrar que conocen los límites da más crédito que vender el
método como infalible.

### 23 · En la práctica — 45 s
`scipy.optimize.newton_krylov`, PETSc SNES, KINSOL. Recorrer rápido la tabla de
síntomas.

### 24 · Conclusiones — 45 s
Las cinco conclusiones, y cerrar con la frase de la caja:
> «La globalización no promete encontrar la raíz. Promete que ‖F‖ nunca aumenta.
> Saber exactamente qué se promete es lo que permite diagnosticar cuando falla.»

## 14. Las seis fórmulas que hay que poder escribir en la pizarra

```
1)  J(x_k) s_k = −F(x_k),      x_{k+1} = x_k + s_k              → §3.3

2)  ‖J s + F‖ ≤ η ‖F‖                             (Newton inexacto)   → §4.3

3)  J(x)·v ≈ [F(x + εv) − F(x)] / ε,  ε ≈ √εₘ     (matriz-libre)      → §5.4

4)  f(x) = ½‖F(x)‖²,     ∇f(x) = J(x)ᵀ F(x)       (función de mérito) → §6.2

5)  ∇f ᵀs ≤ −(1 − η)‖F‖²                          (descenso)          → §6.5

6)  ‖F(x + λs)‖ ≤ (1 − αλ)‖F(x)‖,   α = 10⁻⁴      (Armijo)            → §6.6
```

Extras si hay tiempo o preguntas:

```
η_k = 0.9 (‖F_k‖/‖F_{k−1}‖)^1.618                 (Eisenstat-Walker 2) → §4.6
(I + δJ) s = −δF,   δ_{k+1} = δ_k ‖F_k‖/‖F_{k+1}‖ (Ψtc con regla SER)  → §6.8
```

## 15. Analogías que funcionan

- **Término de forzado.** «Es como corregir un examen. Si el alumno va por la
  pregunta 1 de 20, no tiene sentido corregir con lupa: corriges por encima y sigues.
  Recién al final vale la pena mirar cada coma.» **(→ §4.5)**
- **Krylov matriz-libre.** «GMRES no necesita ver la matriz; solo necesita poder
  preguntarle "¿qué le haces a este vector?". Y esa pregunta se responde con una
  evaluación de F.» **(→ §5.2)**
- **Búsqueda de línea vs región de confianza.** «La búsqueda de línea confía en la
  dirección y negocia la distancia. La región de confianza desconfía del mapa entero
  y decide hasta dónde le cree.» **(→ §6.6, §6.7)**
- **Ψtc.** «En vez de saltar directo al estado final, simula la película desde el
  principio. Al comienzo va en cámara lenta; cuando ve que va bien, acelera.»
  **(→ §6.8)**

## 16. Preguntas probables del docente, con respuesta

**¿Por qué se exige η < 1 y no η ≤ 1, o cualquier otra cota?**  **(→ §6.5)**
Porque de la desigualdad `∇f ᵀs ≤ −(1 − η)‖F‖²` se ve que `η < 1` es exactamente lo
que hace que la cota sea estrictamente negativa. Con `η = 1` la cota es `≤ 0` y ya no
garantiza descenso: el paso podría ser de subida y la búsqueda de línea se quedaría
sin justificación teórica.

**¿Por qué GMRES y no gradientes conjugados, si CG es más barato?**  **(→ §5.3)**
Porque CG exige que la matriz sea simétrica definida positiva y el Jacobiano de un
problema real no lo es: en cuanto hay convección, transporte o acoplamiento
asimétrico entre variables, `J` deja de ser simétrica. En nuestros problemas de
difusión-reacción sí lo es —y lo aprovechamos para el dogleg— pero no se puede
suponer en general.

**¿Qué pasa si el Jacobiano es singular en la solución?**  **(→ §6.7, §6.9)**
Se pierde la convergencia cuadrática (pasa a lineal, con razón 1/2 para raíces
dobles) y la búsqueda de línea se degrada: el paso de Newton se dispara y `λ → 0`.
Ahí la región de confianza es preferible, porque `Δ` acota el paso independientemente
de lo mal condicionada que esté `J`, y Ψtc también, porque `I + δJ` sigue siendo
regular aunque `J` no lo sea.

**¿La globalización garantiza encontrar la raíz?**  **(→ §6.9)**
No. Garantiza que `‖F‖` no aumente y que la sucesión converja a un **punto
estacionario de `f = ½‖F‖²`**. Ese punto puede ser un mínimo local de `f` que no es
raíz de `F`. Lo medimos: sobre Freudenstein-Roth, Newton con Armijo llega a la raíz
desde el 37 % de los puntos y Newton puro desde el 83 %.

**¿Cuándo conviene región de confianza en vez de búsqueda de línea?**  **(→ §6.10)**
Cuando `J` está mal condicionada o es casi singular, y cuando la dirección de Newton
es poco fiable. La búsqueda de línea solo puede moverse sobre la recta de Newton; la
región de confianza puede elegir otra dirección. El precio son dos productos `J·v`
extra por paso y el posible rechazo del paso.

**¿Por qué ε ≈ √εₘ y no ε mucho más chico?**  **(→ §5.5)**
Porque hay dos errores en competencia: el truncamiento crece con `ε` y la cancelación
por redondeo crece al bajarlo. El óptimo está donde se cruzan, en
`√(εₘ‖F‖/‖F″‖)`. Con `ε = 10⁻¹⁶` el error de `J·v` llega a ser del 95 %: se restan
dos números casi iguales y se divide entre algo diminuto.

**¿Cómo se elige el precondicionador?**  **(→ §5.6)**
No hay receta universal: hay que conocer el operador. La regla es que `M` se parezca a
`J` y sea barato de invertir. En nuestro caso de Bratu el laplaciano solo ya sirve,
porque el término no lineal es una perturbación diagonal pequeña. En un problema
dominado por convección habría que incluir esa parte.

**¿Qué diferencia hay entre Ψtc y simplemente resolver el transitorio?**  **(→ §6.8)**
Que a Ψtc no le importa que el transitorio sea físicamente correcto: solo lo usa como
camino hacia el estacionario, y por eso puede alargar el paso de tiempo sin límite
(regla SER) hasta convertirse en Newton puro. Un integrador temporal de verdad tendría
que respetar la precisión temporal.

**¿Por qué el residuo puede ser 10⁻⁹ y la solución estar mal?**  **(→ §9.5)**
Porque el error está acotado por `‖J⁻¹‖·‖F‖`, no por `‖F‖`. Nos pasó en Burgers:
convergimos a `6×10⁻¹⁰` y la solución difiere de la analítica en 0.3, porque la
posición de la capa límite está exponencialmente mal determinada.

**¿Cuál es el costo real de una iteración?**  **(→ §4.5)**
Una evaluación de `F` por cada producto `J·v` de GMRES, más una por cada retroceso de
la búsqueda de línea. Por eso medimos en productos `J·v` y no en iteraciones de
Newton: en dos de nuestros tres experimentos, contar iteraciones llevaba a la
conclusión contraria.

## 17. Checklist antes de exponer

- [ ] Correr los tres ejercicios y confirmar que las figuras están actualizadas.
- [ ] Abrir `presentacion.html` en el navegador del aula **y probar sin internet**
      (las figuras van empotradas en base64, debe funcionar).
- [ ] Probar `f` (pantalla completa) y las flechas en ese navegador.
- [ ] Tener a mano `ej1_newton_vs_globalizado.py` por si piden ver el código.
- [ ] **Plan B si falla el proyector:** las seis fórmulas de §14 en la pizarra
      alcanzan para dar la charla entera.
- [ ] **Plan B si falla Python:** los números están en las tablas de las diapositivas
      8, 12, 17 y 20; no hace falta correr nada en vivo.
- [ ] Repartir `ejercicios_clase/` a los compañeros (los tres `.py` y el README).

## 18. Errores a evitar al exponer

- **No leer las diapositivas.** Están escritas para leerse después; en vivo hay que
  señalar y explicar.
- **No decir «convergencia global» como si significara «mínimo global».** No lo
  significa, y es la confusión más común del tema **(→ §6.1)**.
- **No prometer que la globalización siempre funciona.** La diapositiva 22 existe
  precisamente para eso, y da más crédito que ocultarlo **(→ §6.9)**.
- **No perderse en el pseudocódigo de la 11.** Es un mapa, no un tutorial: 45
  segundos.
- **Si preguntan algo que no saben,** decir qué sí saben y dónde se verificaría.
  «No lo medimos, pero se comprobaría corriendo el ejercicio 3 con ν más chico» es una
  respuesta perfectamente buena.

---

## 19. Fuentes

Ordenadas por dificultad, de menor a mayor. Si va a leer **una sola**, que sea la de
Kelley 2003: es corta, clarísima y está escrita para quien recién entra.

1. **Kelley, C. T. (2003).** *Solving Nonlinear Equations with Newton's Method*. SIAM,
   Fundamentals of Algorithms 1. — El más accesible. Unas 100 páginas, con código.
2. **Nocedal, J. & Wright, S. J. (2006).** *Numerical Optimization*, 2ª ed. Springer.
   — Los capítulos 3 (búsqueda de línea) y 4 (región de confianza) son la referencia
   estándar de la globalización.
3. **Knoll, D. A. & Keyes, D. E. (2004).** «Jacobian-free Newton-Krylov methods: a
   survey of approaches and applications». *J. Comput. Phys.* 193, 357–397. — El
   panorama completo del método matriz-libre y sus aplicaciones.
4. **Kelley, C. T. (1995).** *Iterative Methods for Linear and Nonlinear Equations*.
   SIAM. — Más técnico que el de 2003, cubre también los métodos de Krylov lineales.
5. **Saad, Y. (2003).** *Iterative Methods for Sparse Linear Systems*, 2ª ed. SIAM.
   — Si quiere entender GMRES y los precondicionadores de verdad.
6. **Dembo, R. S., Eisenstat, S. C. & Steihaug, T. (1982).** «Inexact Newton Methods».
   *SIAM J. Numer. Anal.* 19(2), 400–408. — El artículo original del teorema de §4.4.
7. **Eisenstat, S. C. & Walker, H. F. (1996).** «Choosing the forcing terms in an
   inexact Newton method». *SIAM J. Sci. Comput.* 17(1), 16–32. — El artículo original
   de §4.6.
8. **Pawlowski, R. P., Shadid, J. N., Simonis, J. P. & Walker, H. F. (2006).**
   «Globalization techniques for Newton-Krylov methods and applications to the fully
   coupled solution of the Navier-Stokes equations». *SIAM Review* 48(4), 700–721.
   — Comparación experimental de las estrategias de globalización, en el espíritu de
   nuestro ejercicio 3.
9. **Kelley, C. T. & Keyes, D. E. (1998).** «Convergence analysis of pseudo-transient
   continuation». *SIAM J. Numer. Anal.* 35(2), 508–523. — La teoría de Ψtc (§6.8).
