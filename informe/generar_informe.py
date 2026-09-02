# -*- coding: utf-8 -*-
"""
generar_informe.py — Construye Informe_NewtonKrylov_DAT252.docx

Materia : Métodos Numéricos II (DAT-252) — UMSA
Tema    : Estrategias para la convergencia global · Métodos de Newton-Krylov
Autores : Maximiliano Gómez Mallo · Iver Samuel Medina Balboa

El documento se construye SOBRE la plantilla institucional
Gnombres_Dat252.docx, que se abre como ZIP y se reescribe solo en
word/document.xml. Así hereda tal cual la tipografía, los márgenes y el tamaño
de página de la plantilla, sin depender de python-docx.

Las ecuaciones desplegadas se rinden con matplotlib (mathtext) a PNG
transparente y se insertan centradas: generar OMML a mano es frágil y esto se
ve mejor.

    python3 informe/generar_informe.py
"""

import io
import pathlib
import shutil
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from docx_min import Documento, EMU_POR_CM, parrafo, run, runs_con_marcas, salto_pagina

warnings.filterwarnings("ignore")

RAIZ = pathlib.Path(__file__).resolve().parent.parent
PLANTILLA = RAIZ / "Gnombres_Dat252.docx"
SALIDA = RAIZ / "Informe_NewtonKrylov_DAT252.docx"
FIG_EXPO = RAIZ / "ejercicios_exposicion" / "figuras"
FIG_CLASE = RAIZ / "ejercicios_clase"
CACHE = pathlib.Path(__file__).resolve().parent / "figuras"
LOGO_ORIGEN = pathlib.Path("/home/max1/ml/jogo/rl_combat/informe/figuras/Logo_Umsa.png")

AZUL = "1F4E79"
ANCHO_TEXTO_CM = 16.5          # 8.5" − 2" de márgenes


# =============================================================================
# Ecuaciones desplegadas
# =============================================================================
_n_eq = 0


def ecuacion(doc, tex, tam=13, escala=1.0):
    """Rinde $tex$ a PNG transparente y lo inserta centrado.

    Se rinde a 300 dpi y se coloca a su tamaño natural en pulgadas, de manera
    que `tam` puntos en la figura son `tam` puntos en la página.
    """
    global _n_eq
    _n_eq += 1
    destino = CACHE / f"eq{_n_eq:02d}.png"
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, f"${tex}$", fontsize=tam, color="#111111")
    fig.savefig(destino, dpi=300, bbox_inches="tight", pad_inches=0.03,
                transparent=True)
    plt.close(fig)

    datos = destino.read_bytes()
    px_w = int.from_bytes(datos[16:20], "big")
    ancho_emu = int(px_w / 300 * 914400 * escala)
    tope = int(ANCHO_TEXTO_CM * EMU_POR_CM)
    if ancho_emu > tope:
        ancho_emu = tope
    doc.imagen(destino, ancho_emu=ancho_emu, despues=180)


def figura(doc, nombre, pie, ancho_cm=15.0):
    for carpeta in (FIG_EXPO, FIG_CLASE):
        ruta = carpeta / nombre
        if ruta.exists():
            doc.imagen(ruta, ancho_cm=ancho_cm, despues=60)
            doc.pie_figura(pie)
            return
    raise FileNotFoundError(
        f"Falta {nombre}. Corré primero los ejercicios de ejercicios_exposicion/.")


# =============================================================================
# PORTADA
# =============================================================================
def portada(doc):
    for texto, sz in [("UNIVERSIDAD MAYOR DE SAN ANDRÉS", 36),
                      ("FACULTAD DE CIENCIAS PURAS Y NATURALES", 28),
                      ("CARRERA DE INFORMÁTICA", 28)]:
        doc.add(parrafo(run(texto, negrita=True, sz=sz), jc="center",
                        despues=0, linea=None))

    doc.espacio(200)
    if LOGO_ORIGEN.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        logo = CACHE / "Logo_Umsa.png"
        if not logo.exists():
            shutil.copy(LOGO_ORIGEN, logo)
        doc.imagen(logo, ancho_cm=2.9, despues=200)
    else:
        doc.espacio(400)

    doc.regla(AZUL, 12, antes=0, despues=200)
    doc.add(parrafo(run("Estrategias para la Convergencia Global",
                        negrita=True, sz=40), jc="center", despues=80,
                    linea=None))
    doc.add(parrafo(run("Métodos de Newton-Krylov", sz=28), jc="center",
                    despues=180, linea=None))
    doc.regla(AZUL, 12, antes=0, despues=320)

    doc.tabla([
        ["MATERIA:", "Métodos Numéricos II"],
        ["SIGLA:", "DAT-252"],
        ["DOCENTE:", "M.Sc. Carlos Mullisaca Choque"],
        ["ESTUDIANTES:", "Gómez Mallo, Maximiliano\nMedina Balboa, Iver Samuel"],
        ["FECHA DE ENTREGA:", "La Paz, Bolivia — 1 de Septiembre de 2026"],
    ], anchos=[2600, 6760], encabezado=False, sz=21,
        relleno_alt="F7FAFC")
    doc.add(salto_pagina())


# =============================================================================
# CUERPO
# =============================================================================
def introduccion(doc):
    doc.titulo_seccion("1. INTRODUCCIÓN")
    doc.texto(
        "Buena parte de los problemas de la ingeniería y de la física computacional "
        "terminan, después de discretizar, en el mismo lugar: un sistema de "
        "ecuaciones **no lineales** que hay que resolver.")
    ecuacion(doc, r"F(x) = 0, \qquad F : \mathbb{R}^n \rightarrow \mathbb{R}^n")
    doc.texto(
        "El método de Newton es la herramienta natural para atacarlo, y su "
        "convergencia cuadrática es una de las propiedades más citadas del análisis "
        "numérico. Sin embargo, esa convergencia es **local**: el teorema garantiza "
        "que existe un radio __δ__ alrededor de la solución dentro del cual el método "
        "funciona, pero no dice cuánto vale ese radio ni cómo estimarlo. Fuera de él, "
        "el método puede estancarse, oscilar o divergir.")
    doc.texto(
        "A eso se le suma un problema de escala. En los sistemas que interesan, __n__ "
        "no es 2 ni 10: es 10⁴, 10⁶ o más. Una malla de 100×100×100 nodos en tres "
        "dimensiones ya da un millón de incógnitas, y el Jacobiano correspondiente "
        "tendría 10¹² entradas — unos ocho terabytes si se almacenara denso. Formarlo "
        "es imposible, y factorizarlo, más todavía.")
    doc.texto(
        "Los **métodos de Newton-Krylov** resuelven las dos dificultades a la vez. "
        "Sustituyen la resolución exacta del sistema lineal de Newton por un método "
        "iterativo de subespacios de Krylov que solo necesita productos matriz-vector, "
        "y aproximan esos productos con diferencias finitas de la propia función "
        "residuo, sin construir nunca el Jacobiano. Sobre esa base, las **estrategias "
        "de convergencia global** —búsqueda de línea, región de confianza y "
        "continuación pseudo-transitoria— son las que convierten un método rápido "
        "pero frágil en un método utilizable desde puntos iniciales arbitrarios.")
    doc.texto(
        "Este informe expone la teoría de esas estrategias, las implementa desde cero "
        "en Python y las compara experimentalmente sobre tres problemas de dificultad "
        "creciente. Todas las cifras que aparecen en las secciones de resultados "
        "provienen de las corridas reales del código que se adjunta, no de la "
        "literatura.")


def marco_teorico(doc):
    doc.titulo_seccion("2. MARCO TEÓRICO")

    doc.titulo_sub("2.1. El método de Newton y su convergencia local")
    doc.texto(
        "El método de Newton linealiza el residuo alrededor del iterado actual y "
        "resuelve el modelo lineal resultante:")
    ecuacion(doc, r"J(x_k)\,s_k = -F(x_k), \qquad x_{k+1} = x_k + s_k")
    doc.texto(
        "donde __J(x) = ∂F/∂x__ es el Jacobiano. El teorema clásico de convergencia "
        "establece que si __F__ es continuamente diferenciable, __J(x*)__ es no "
        "singular y __J__ es Lipschitz continua en un entorno de la solución __x*__, "
        "entonces existe __δ > 0__ tal que, para todo punto inicial con "
        "‖x₀ − x*‖ < δ, la sucesión converge y satisface")
    ecuacion(doc, r"\|x_{k+1}-x^*\| \;\leq\; C\,\|x_k - x^*\|^2")
    doc.texto(
        "La convergencia cuadrática es real y espectacular: el número de cifras "
        "correctas se duplica en cada paso. Pero la hipótesis crítica está en las "
        "cuatro palabras **«existe δ > 0 tal que»**. El teorema no proporciona ninguna "
        "manera de calcular __δ__, y en problemas mal condicionados o fuertemente no "
        "lineales ese radio puede ser muy pequeño. Fuera de él el teorema no promete "
        "absolutamente nada, y «nada» incluye la divergencia.")

    doc.titulo_sub("2.2. Métodos de Newton inexactos y el término de forzado")
    doc.texto(
        "Cuando __xₖ__ está lejos de la solución, el modelo lineal __J s = −F__ es "
        "una aproximación pobre de __F__. Resolverlo con precisión de máquina es, por "
        "lo tanto, un desperdicio: se gasta trabajo en obtener con dieciséis cifras la "
        "solución de un modelo que solo vale una o dos. Dembo, Eisenstat y Steihaug "
        "(1982) formalizaron esta idea con los **métodos de Newton inexactos**, que "
        "aceptan cualquier paso __sₖ__ que satisfaga")
    ecuacion(doc, r"\|J(x_k)\,s_k + F(x_k)\| \;\leq\; \eta_k\,\|F(x_k)\|")
    doc.texto(
        "El parámetro __ηₖ ∈ [0,1)__ se denomina **término de forzado**. El teorema "
        "de Dembo-Eisenstat-Steihaug relaciona su elección con el orden de "
        "convergencia local, de acuerdo con la Tabla 1.")
    doc.tabla([
        ["Elección del término de forzado", "Orden de convergencia local"],
        ["ηₖ ≤ η < 1 constante", "lineal"],
        ["ηₖ → 0", "superlineal"],
        ["ηₖ = O(‖Fₖ‖)", "cuadrática (se recupera Newton)"],
    ], anchos=[5200, 4160])
    doc.pie_figura("**Tabla 1.** Teorema de Dembo, Eisenstat y Steihaug (1982): "
                   "el término de forzado controla el orden de convergencia.")
    doc.texto(
        "La lectura práctica es que **se puede ser perezoso sin pagar ningún precio en "
        "velocidad de convergencia**, siempre que la pereza se ajuste al ritmo "
        "correcto. Eisenstat y Walker (1996) propusieron dos reglas para hacerlo "
        "automáticamente. La segunda, que es la que usamos, mide cuánto se redujo el "
        "residuo en el paso anterior:")
    ecuacion(doc, r"\eta_k = \gamma \left( \frac{\|F_k\|}{\|F_{k-1}\|} \right)^{\alpha}, "
                  r"\qquad \gamma = 0.9, \quad \alpha = \frac{1+\sqrt{5}}{2} \approx 1.618")
    doc.texto(
        "Si el residuo bajó mucho, el modelo lineal está describiendo bien a __F__ y "
        "conviene apretar; si bajó poco, no vale la pena invertir en resolverlo. Se "
        "añade una salvaguarda, __ηₖ = max(ηₖ, γ ηₖ₋₁^α)__ cuando esta cota supera "
        "0.1, para impedir que __η__ se desplome por un paso afortunado aislado.")

    doc.titulo_sub("2.3. Subespacios de Krylov y el producto Jacobiano-vector")
    doc.texto(
        "Admitida la resolución inexacta, el sistema lineal puede resolverse con un "
        "método iterativo. Los métodos de Krylov buscan la solución de __A s = b__ "
        "dentro del subespacio")
    ecuacion(doc, r"\mathcal{K}_m(A,b) = \mathrm{span}\{\, b,\; Ab,\; A^2b,\; \ldots,\; A^{m-1}b \,\}")
    doc.texto(
        "y la observación decisiva es que **construir esa base solo requiere multiplicar "
        "por A**. Nunca hace falta la matriz entrada por entrada, ni factorizarla, ni "
        "almacenarla: basta con una función que, dado __v__, devuelva __A·v__. "
        "Se emplea **GMRES** con reinicio, porque el Jacobiano de un problema con "
        "convección o transporte no es simétrico y el método de gradientes conjugados "
        "no es aplicable.")
    doc.texto(
        "Y ese producto es una derivada direccional, que se aproxima con **una sola "
        "evaluación de F**:")
    ecuacion(doc, r"J(x)\,v \;\approx\; \frac{F(x+\varepsilon v) - F(x)}{\varepsilon}")
    doc.texto(
        "Esto es todo el método **Jacobian-Free Newton-Krylov** (JFNK). La elección de "
        "__ε__ equilibra dos errores opuestos: el de truncamiento, del orden de "
        "__(ε/2)‖F″‖‖v‖²__, que crece con __ε__; y el de cancelación por redondeo, del "
        "orden de __εₘ‖F‖/(ε‖v‖)__ —con __εₘ ≈ 2.2×10⁻¹⁶__ el épsilon de máquina—, "
        "que crece al disminuirlo. Igualándolos:")
    ecuacion(doc, r"\varepsilon_{\mathrm{opt}} \approx \sqrt{\frac{\epsilon_{maq}\,\|F\|}{\|F''\|}}"
                  r" \;\approx\; \sqrt{\epsilon_{maq}} \approx 1.5\times10^{-8}")
    doc.texto(
        "La aproximación final supone que __‖F‖__ y __‖F″‖__ son del mismo orden. "
        "Cuando no lo son —por ejemplo si el operador está dominado por su parte "
        "lineal— el óptimo se desplaza, como se comprueba en la Sección 5.2. En la "
        "práctica se usa __√εₘ__ de todos modos: no es óptimo, pero nunca es "
        "catastrófico y no exige conocer __‖F″‖__.")
    doc.texto(
        "El precio de este enfoque es que el Jacobiano solo se conoce con unas ocho "
        "cifras significativas. Cerca de la solución la convergencia deja de ser "
        "exactamente cuadrática, y no puede exigirse un residuo final mucho menor que "
        "__‖F(x₀)‖·√εₘ__.")

    doc.titulo_sub("2.4. Precondicionamiento")
    doc.texto(
        "GMRES converge rápido cuando los autovalores del operador están agrupados. El "
        "Jacobiano de una ecuación en derivadas parciales hereda del laplaciano "
        "discreto un número de condición del orden de __1/h² = O(n²)__, lo que hace "
        "que el número de iteraciones crezca con el refinamiento de la malla. El "
        "remedio es resolver el sistema precondicionado:")
    ecuacion(doc, r"M^{-1} J\,s = -M^{-1}F")
    doc.texto(
        "El precondicionador __M__ debe cumplir dos condiciones simultáneas y en "
        "tensión: parecerse a __J__ lo suficiente como para agrupar el espectro, y ser "
        "barato de invertir. Las opciones habituales son la parte lineal del operador, "
        "una factorización incompleta ILU, multigrid, descomposición de dominios, o el "
        "Jacobiano congelado de alguna iteración anterior. **No existe un "
        "precondicionador universal**: elegirlo exige conocer el operador.")


def globalizacion(doc):
    doc.add(salto_pagina())
    doc.titulo_seccion("3. ESTRATEGIAS PARA LA CONVERGENCIA GLOBAL")
    doc.texto(
        "Las dos secciones anteriores hacen que el método sea **posible** en problemas "
        "grandes. Esta sección es la que lo hace **confiable**.")

    doc.titulo_sub("3.1. La función de mérito y la dirección de descenso")
    doc.texto(
        "El método de Newton no tiene ninguna noción interna de «estar mejorando»: "
        "acepta el paso completo sin preguntarse si el resultado es mejor que el punto "
        "de partida. Hay que darle una medida de progreso, y la elección estándar es "
        "la **función de mérito**")
    ecuacion(doc, r"f(x) = \frac{1}{2}\|F(x)\|^2 = \frac{1}{2}F(x)^{T}F(x),"
                  r"\qquad \nabla f(x) = J(x)^{T}F(x)")
    doc.texto(
        "Es escalar, no negativa, y se anula exactamente en las raíces de __F__. Con "
        "ella, resolver __F(x) = 0__ se convierte en minimizar __f__, y toda la "
        "maquinaria de la optimización sin restricciones queda disponible.")
    doc.texto(
        "El resultado que sostiene todo lo que sigue responde a una pregunta: ¿es el "
        "paso de Newton inexacto una dirección de descenso para __f__? Si "
        "__J s = −F + r__ con __‖r‖ ≤ η‖F‖__, entonces")
    ecuacion(doc, r"\nabla f(x)^{T}s = F^{T}Js = F^{T}(-F+r) = -\|F\|^2 + F^{T}r"
                  r" \;\leq\; -\|F\|^2 + \|F\|\|r\| \;\leq\; -(1-\eta)\|F\|^2")
    doc.texto(
        "Por lo tanto, **si η < 1 y F ≠ 0, el producto ∇f ᵀs es estrictamente "
        "negativo**: el paso de Newton inexacto siempre apunta cuesta abajo, y en "
        "consecuencia siempre existe algún __λ > 0__ que reduce __f__. La búsqueda de "
        "línea nunca se queda sin opciones.")
    doc.texto(
        "Esto explica además por qué se exige __η < 1__ y no cualquier cota. Si "
        "__η ≥ 1__ la desigualdad anterior se vuelve __≤ 0__ y deja de garantizar nada: "
        "el paso podría ser de subida. **La condición η < 1 es exactamente lo que hace "
        "compatibles el Newton inexacto y la globalización.**")

    doc.titulo_sub("3.2. Búsqueda de línea con retroceso (Armijo)")
    doc.texto(
        "La estrategia más simple mantiene fija la dirección de Newton y negocia solo "
        "su longitud: __xₖ₊₁ = xₖ + λₖ sₖ__. El paso __λₖ__ se acepta cuando "
        "satisface la **condición de Armijo**, escrita aquí en términos de la norma del "
        "residuo:")
    ecuacion(doc, r"\|F(x+\lambda s)\| \;\leq\; (1-\alpha\lambda)\,\|F(x)\|,"
                  r"\qquad \alpha = 10^{-4}")
    doc.texto(
        "Se prueba primero __λ = 1__, de modo que cerca de la solución el paso completo "
        "se acepta y la convergencia cuadrática se conserva íntegra. Si no se cumple, "
        "se retrocede usando interpolación cuadrática sobre la información ya "
        "disponible, con salvaguardas __λ nuevo ∈ [0.1λ, 0.5λ]__ para que el retroceso "
        "no sea ni demasiado tímido ni demasiado brusco, y se abandona si __λ__ cae por "
        "debajo de un umbral.")
    doc.texto(
        "El valor __α = 10⁻⁴__ exige una mejora casi simbólica: lo justo para que el "
        "teorema de convergencia funcione, sin frenar al método. Su virtud principal es "
        "el costo: **cada retroceso cuesta una evaluación de F, no un sistema lineal "
        "nuevo**. Su limitación es que solo puede moverse sobre la recta de Newton; si "
        "esa dirección es mala, no hay nada que hacer.")

    doc.titulo_sub("3.3. Región de confianza: dogleg y Steihaug-CG")
    doc.texto(
        "La región de confianza invierte la pregunta. En lugar de «¿cuánto de este paso "
        "acepto?», pregunta «¿hasta dónde le creo al modelo lineal?», y resuelve un "
        "problema restringido:")
    ecuacion(doc, r"\min_{s}\; \frac{1}{2}\|F + Js\|^2 \quad \mathrm{sujeto\;a}\quad \|s\| \leq \Delta")
    doc.texto(
        "El paso se acepta o se rechaza según la razón entre la reducción real y la "
        "predicha por el modelo, y el radio __Δ__ se ajusta con ese mismo indicador:")
    ecuacion(doc, r"\rho = \frac{f(x)-f(x+s)}{f(x)-m(s)}; \qquad"
                  r"\rho < \frac{1}{4} \Rightarrow \Delta \downarrow, \qquad"
                  r"\rho > \frac{3}{4} \Rightarrow \Delta \uparrow")
    doc.texto(
        "Como __Δ__ acota explícitamente el tamaño del paso, la región de confianza "
        "**sigue funcionando aunque J sea singular**, que es precisamente donde la "
        "búsqueda de línea se rinde. Ahí está su ventaja real.")
    doc.texto(
        "El subproblema restringido admite dos tratamientos, y la elección entre ellos "
        "no es indiferente en el contexto matriz-libre. El **CG truncado de Steihaug** "
        "aplica gradientes conjugados a las ecuaciones normales __JᵀJ s = −JᵀF__, lo "
        "que **eleva al cuadrado el número de condición**. Con un Jacobiano conocido "
        "solo por diferencias finitas, con ruido relativo del orden de 10⁻⁸, esto lo "
        "vuelve inviable: en nuestras pruebas sobre Bratu el método se estancó "
        "alrededor de 10⁻⁶ y dejó de progresar. El **dogleg de Powell**, en cambio, "
        "interpola entre el punto de Cauchy y el paso de Newton que GMRES ya calculó:")
    ecuacion(doc, r"s_C = -\frac{\|g\|^2}{\|Jg\|^2}\,g, \qquad g = J^{T}F;"
                  r"\qquad s = s_C + \tau\,(s_N - s_C),\quad \|s\| = \Delta")
    doc.texto(
        "Cuesta dos productos matriz-vector adicionales y es lo que se usa en la "
        "práctica (Pawlowski et al., 2006). Requiere __Jᵀ__, que las diferencias "
        "finitas hacia adelante no proporcionan; se resuelve aprovechando la simetría "
        "de __J__ cuando existe —cierta de forma exacta en operadores de "
        "difusión-reacción discretizados con diferencias centradas— o formando el "
        "Jacobiano si la dimensión es pequeña. **Esta es una limitación real del "
        "enfoque matriz-libre, no un detalle de implementación.**")

    doc.titulo_sub("3.4. Continuación pseudo-transitoria")
    doc.texto(
        "La tercera estrategia parte de una idea distinta y más física: el estado "
        "estacionario que se busca es el límite de un transitorio, de modo que puede "
        "obtenerse simulándolo. Se integra")
    ecuacion(doc, r"\frac{dx}{dt} = -F(x), \qquad x(0) = x_0")
    doc.texto("con Euler implícito y paso __δ__, lo que da en cada paso el sistema")
    ecuacion(doc, r"(I + \delta J)\,s = -\delta F")
    doc.texto(
        "Los dos extremos son ilustrativos. Con __δ → 0__ se obtiene __s = −δF__, es "
        "decir máximo descenso fuertemente amortiguado: robustísimo y lentísimo. Con "
        "__δ → ∞__ se recupera __J s = −F__, o sea Newton puro: rápido y frágil. La "
        "**regla SER** (Switched Evolution Relaxation) interpola automáticamente entre "
        "ambos:")
    ecuacion(doc, r"\delta_{k+1} = \delta_k \cdot \frac{\|F_k\|}{\|F_{k+1}\|}")
    doc.texto(
        "El paso de tiempo se alarga solo cuando el residuo baja, de modo que el método "
        "empieza amortiguado y termina siendo Newton. Nótese que el sistema se escribe "
        "multiplicado por __δ__: la forma __(I/δ + J)s = −F__ es equivalente en "
        "aritmética exacta pero está peor condicionada en ambos extremos.")
    doc.texto(
        "Dos advertencias que conviene retener. La primera es que el método exige que "
        "el flujo __dx/dt = −F(x)__ sea **estable**, lo que obliga a elegir el signo del "
        "residuo de manera que __F′__ resulte definida positiva; con el signo contrario "
        "el método diverge sistemáticamente. La segunda es que __δ₀__ es un parámetro "
        "real, no un detalle: demasiado pequeño sobre-amortigua y el método se arrastra, "
        "demasiado grande equivale a Newton puro y se pierde la robustez. Se cuantifica "
        "en la Sección 5.3.")

    doc.titulo_sub("3.5. Comparación de las tres estrategias")
    doc.tabla([
        ["", "Búsqueda de línea", "Región de confianza", "Ψtc"],
        ["Qué ajusta", "longitud λ del paso", "radio Δ de confianza", "paso de tiempo δ"],
        ["Dirección", "siempre la de Newton", "puede cambiarla", "interpola Newton–descenso"],
        ["Costo extra por paso", "1 evaluación de F por retroceso", "2 productos J·v y posible rechazo", "ninguno"],
        ["Si J es singular", "λ → 0, se atasca", "Δ lo acota, continúa", "I + δJ es regular"],
        ["Parámetros a ajustar", "α (no se toca)", "Δ₀ y umbrales de ρ", "δ₀ (sí importa)"],
        ["Cuándo usarla", "por defecto", "J mal condicionada", "estados estacionarios de EDPs"],
    ], anchos=[2300, 2400, 2500, 2160], sz=17)
    doc.pie_figura("**Tabla 2.** Las tres estrategias de globalización comparadas. "
                   "También se combinan entre sí: es habitual usar continuación en un "
                   "parámetro físico para generar un buen punto inicial y aplicar "
                   "después Newton-Krylov con búsqueda de línea.")


def implementacion(doc):
    doc.add(salto_pagina())
    doc.titulo_seccion("4. IMPLEMENTACIÓN")

    doc.titulo_sub("4.1. Estructura del código")
    doc.texto(
        "Todo el método se implementó desde cero en Python, con NumPy y SciPy como "
        "únicas dependencias. El núcleo está en `ejercicios_exposicion/nk_lib.py` y "
        "contiene las piezas descritas en las secciones 2 y 3: el producto "
        "Jacobiano-vector por diferencias finitas, los términos de forzado de "
        "Eisenstat-Walker, la búsqueda de línea de Armijo con interpolación "
        "cuadrática, el CG truncado de Steihaug, el paso de dogleg y el bucle de "
        "Newton-Krylov que los combina. De SciPy se reutiliza únicamente `gmres` a "
        "través de un `LinearOperator`, que es la manera correcta de hacerlo: "
        "reimplementar GMRES no aportaría nada al tema del trabajo.")
    doc.tabla([
        ["Archivo", "Contenido"],
        ["`nk_lib.py`", "Núcleo compartido: J·v, forzado, Armijo, Steihaug, dogleg, Ψtc"],
        ["`ej1_newton_vs_globalizado.py`", "Por qué hace falta globalizar (1D y 2×2)"],
        ["`ej2_bratu1d_newton_krylov.py`", "El método completo sobre Bratu 1D"],
        ["`ej3_comparativa_globalizacion.py`", "Las cuatro estrategias sobre Burgers 1D"],
    ], anchos=[3400, 5960], sz=17)
    doc.pie_figura("**Tabla 3.** Organización del código de la exposición.")

    doc.titulo_sub("4.2. Medición del costo")
    doc.texto(
        "Una decisión metodológica atraviesa todos los experimentos: **el costo se mide "
        "en productos J·v, no en iteraciones de Newton**. La razón es que en un método "
        "matriz-libre cada producto J·v cuesta exactamente una evaluación de __F__, y "
        "en un problema real esa evaluación es una simulación completa. Contar "
        "iteraciones externas mide lo que se ve, no lo que se paga, y como se muestra "
        "en la Sección 5.2 puede llevar a conclusiones exactamente invertidas.")

    doc.titulo_sub("4.3. Los problemas de prueba")
    doc.tabla([
        ["Problema", "Ecuación", "Qué pone a prueba"],
        ["Arcotangente", "arctan(x) = 0", "la localidad de Newton en el caso más simple"],
        ["Sistema 2×2", "x₁²+x₂²−2 = 0 ; exp(x₁−1)+x₂³−2 = 0",
         "trayectorias, cuencas y la singularidad de J"],
        ["Bratu 1D", "u″ + λ·exp(u) = 0, u(0)=u(1)=0",
         "términos de forzado y precondicionamiento"],
        ["Burgers 1D", "−ν u″ + u u′ = 0, u(0)=1, u(1)=−1",
         "las cuatro estrategias de globalización"],
    ], anchos=[1900, 3660, 3800], sz=17)
    doc.pie_figura("**Tabla 4.** Problemas de prueba. La dificultad crece de arriba "
                   "hacia abajo, y con ella el margen entre globalizar y no hacerlo.")


def resultados(doc):
    doc.add(salto_pagina())
    doc.titulo_seccion("5. RESULTADOS Y DISCUSIÓN")

    # ---------------------------------------------------------------
    doc.titulo_sub("5.1. La necesidad de globalizar")
    doc.texto(
        "El primer experimento resuelve __arctan(x) = 0__, cuya raíz es __x* = 0__. La "
        "función es suave y monótona, pero se aplana lejos del origen: la tangente "
        "cruza el eje cada vez más lejos y el paso de Newton sobrepasa la raíz de "
        "manera creciente. El umbral es __|x₀| ≈ 1.3917__; por encima de él, Newton "
        "puro diverge.")
    figura(doc, "fig1a_arctan.png",
           "**Figura 1.** Izquierda: construcción geométrica del disparo. Cada tangente "
           "cruza el eje más lejos que la anterior. Derecha: desde el mismo x₀ = 2, "
           "Newton puro diverge y Newton con búsqueda de línea converge en 5 "
           "iteraciones.", ancho_cm=15.5)
    doc.texto(
        "El segundo experimento usa el sistema 2×2, cuya raíz es __x* = (1, 1)__. Desde "
        "__x₀ = (2.0, 0.5)__ las tres estrategias parten del mismo punto y usan el mismo "
        "paso de Newton; lo único que cambia es cuánto de ese paso se acepta.")
    doc.tabla([
        ["Estrategia", "¿Converge?", "Iteraciones", "‖F‖ final", "Evals. de F"],
        ["Newton puro", "no", "—", "desbordamiento", "138"],
        ["Búsqueda de línea (Armijo)", "sí", "7", "9.9×10⁻¹²", "108"],
        ["Región de confianza (dogleg)", "sí", "6", "0 (exacto)", "799"],
    ], anchos=[3200, 1500, 1600, 1700, 1360], sz=17)
    doc.pie_figura("**Tabla 5.** Sistema 2×2 desde x₀ = (2.0, 0.5). Newton puro "
                   "desborda; ambas globalizaciones alcanzan la raíz.")
    figura(doc, "fig1b_trayectorias.png",
           "**Figura 2.** Izquierda: trayectorias sobre las curvas de nivel de ‖F‖; las "
           "líneas blancas son F₁ = 0 y F₂ = 0, que se cortan en la raíz. Derecha: "
           "historia del residuo. La curva roja crece hasta desbordar.", ancho_cm=15.5)
    doc.texto(
        "Conviene registrar un matiz que el experimento reveló y que no aparece en las "
        "presentaciones habituales del tema: **en sistemas pequeños y suaves, Newton "
        "puro ya es bastante robusto**. Con un presupuesto de 12 iteraciones y sobre "
        "2601 puntos iniciales, Newton puro converge desde el 86.6 % y Newton con "
        "Armijo desde el 95.3 %. La diferencia es real pero moderada. El abismo entre "
        "globalizar y no hacerlo aparece cuando el sistema proviene de discretizar una "
        "EDP, y eso se documenta en la Sección 5.3.")

    # ---------------------------------------------------------------
    doc.add(salto_pagina())
    doc.titulo_sub("5.2. Término de forzado, oversolving y precondicionamiento")
    doc.texto(
        "El segundo experimento resuelve la ecuación de Bratu, __u″ + λ·exp(u) = 0__ con "
        "__u(0) = u(1) = 0__ y __λ = 3__, discretizada con __N = 250__ nodos interiores. "
        "Se comparan cinco maneras de elegir __η__, manteniendo idénticos el paso de "
        "Newton, GMRES y la globalización.")
    doc.tabla([
        ["Término de forzado", "Iters. de Newton", "Productos J·v", "Tiempo", "‖F‖ final"],
        ["η = 10⁻¹ fijo", "7", "9 854", "1.8 s", "2.2×10⁻¹⁰"],
        ["η = 10⁻³ fijo", "4", "21 332", "3.9 s", "9.5×10⁻¹²"],
        ["η = 10⁻¹² fijo", "4", "24 804", "4.4 s", "1.1×10⁻¹¹"],
        ["Eisenstat-Walker 1", "8", "10 495", "1.9 s", "5.0×10⁻¹⁰"],
        ["Eisenstat-Walker 2", "8", "9 707", "1.8 s", "5.1×10⁻¹⁰"],
    ], anchos=[2700, 1900, 1800, 1300, 1660], sz=17)
    doc.pie_figura("**Tabla 6.** Bratu 1D, N = 250, λ = 3, sin precondicionador. "
                   "La fila más cara en trabajo es la que menos iteraciones de Newton usa.")
    doc.texto(
        "El resultado es el fenómeno del **oversolving**, y es contraintuitivo hasta que "
        "se mira la columna correcta. Resolver el sistema lineal con doce cifras "
        "(η = 10⁻¹²) reduce las iteraciones de Newton de ocho a cuatro, pero multiplica "
        "por **2.6** el trabajo total. Eisenstat-Walker 2, que usa el doble de "
        "iteraciones externas, es el más barato de los cinco sin que se le indique "
        "nada: empieza con __η = 0.9__ y solo aprieta cuando el residuo demuestra que el "
        "modelo lineal ya sirve.")
    figura(doc, "fig2c_solucion_eta.png",
           "**Figura 3.** Izquierda: solución de Bratu para λ = 3. Derecha: el término "
           "de forzado elegido por Eisenstat-Walker sigue al residuo relativo sin "
           "intervención del usuario.", ancho_cm=15.0)
    doc.texto(
        "El efecto del precondicionamiento es de otro orden de magnitud. Tomando como "
        "__M__ el laplaciano discreto —tridiagonal, factorizado una sola vez con LU "
        "disperso— el trabajo se reduce en un factor de varios cientos, como muestra la "
        "Tabla 7.")
    doc.tabla([
        ["Término de forzado", "J·v sin M", "J·v con M", "Ganancia"],
        ["η = 10⁻¹ fijo", "9 854", "19", "519×"],
        ["η = 10⁻³ fijo", "21 332", "24", "889×"],
        ["η = 10⁻¹² fijo", "24 804", "24 603", "1.0×"],
        ["Eisenstat-Walker 2", "9 707", "22", "441×"],
    ], anchos=[3000, 2200, 2200, 1960], sz=17)
    doc.pie_figura("**Tabla 7.** Efecto del precondicionador laplaciano sobre el mismo "
                   "problema. Obsérvese la fila que no mejora.")
    doc.texto(
        "Dos observaciones sobre esta tabla. La primera es que la ganancia es tan "
        "grande porque el laplaciano **es** casi todo el Jacobiano en este problema: el "
        "término no lineal __−h²λ·exp(u)__ es una perturbación diagonal pequeña. En un "
        "problema dominado por convección el precondicionador tendría que incluir esa "
        "parte y la ganancia sería menor.")
    doc.texto(
        "La segunda es la fila de __η = 10⁻¹²__, que **no mejora en absoluto**. Se le "
        "está pidiendo a GMRES un residuo relativo por debajo de lo que la aritmética "
        "de doble precisión puede entregar sobre este operador; no puede alcanzarlo, "
        "agota su presupuesto de iteraciones en cada paso de Newton, y el "
        "precondicionador no lo salva. Pedir más precisión de la que existe no acelera "
        "nada: solo consume trabajo.")
    figura(doc, "fig2b_residual_vs_trabajo.png",
           "**Figura 4.** Convergencia frente al trabajo real (productos J·v) en ambos "
           "ejes logarítmicos. Con precondicionador todas las curvas colapsan a unas "
           "pocas decenas de productos, salvo la de η = 10⁻¹².", ancho_cm=15.5)
    doc.texto(
        "En cuanto al paso __ε__ de la diferencia finita, la medición confirmó el "
        "comportamiento en «V» previsto por la teoría, pero con un matiz: el mínimo "
        "medido no cayó exactamente en __√εₘ__, sino unos dos órdenes de magnitud "
        "por encima. La explicación está en la fórmula completa de la Sección 2.3: en "
        "este problema el operador está dominado por su parte lineal, __‖F″‖__ es "
        "pequeña frente a __‖F‖__ y el óptimo se desplaza. Usar __√εₘ__ de todos "
        "modos da un error unas cien veces mayor que el óptimo — y resulta "
        "**irrelevante**, porque Newton inexacto ya tolera un residuo lineal de "
        "__η‖F‖__ con __η ~ 10⁻²__: un Jacobiano con siete cifras correctas le sobra.")
    figura(doc, "fig2a_epsilon.png",
           "**Figura 5.** Error relativo de J·v frente al paso ε, medido contra el "
           "Jacobiano analítico. A la izquierda domina la cancelación por redondeo; a "
           "la derecha, el truncamiento.", ancho_cm=11.0)

    # ---------------------------------------------------------------
    doc.add(salto_pagina())
    doc.titulo_sub("5.3. Las cuatro estrategias sobre un problema convectivo")
    doc.texto(
        "El experimento central emplea la ecuación de Burgers estacionaria, "
        "__−ν u″ + u u′ = 0__ con __u(0) = 1__ y __u(1) = −1__ y __ν = 0.01__. Su "
        "solución presenta una capa límite interna: una transición casi vertical en "
        "__x = ½__ de anchura del orden de __2ν__. Es el prototipo de lo que aparece al "
        "resolver Navier-Stokes estacionario.")
    doc.texto(
        "Conviene subrayar que el número de Péclet de celda vale __h/ν = 0.5 < 2__, de "
        "modo que las diferencias centradas no producen oscilaciones espurias: **la "
        "dificultad del problema es genuinamente la no linealidad convectiva, no una "
        "discretización inestable**. Si lo fuera, ninguna estrategia de globalización "
        "lo arreglaría.")
    figura(doc, "fig3a_problema.png",
           "**Figura 6.** La capa límite se afila al disminuir ν. Las rectas grises son "
           "los siete puntos iniciales u₀ = a(1−2x) empleados en la comparación.",
           ancho_cm=12.5)
    doc.texto(
        "Se resolvió el problema desde siete puntos iniciales con cada una de las "
        "cuatro estrategias. Todas las corridas comparten el paso de Newton, GMRES, el "
        "precondicionador y el término de forzado Eisenstat-Walker 2; lo único que "
        "cambia es la globalización.")
    doc.tabla([
        ["Estrategia", "Converge desde", "Tasa", "Iters. (mediana)", "J·v (mediana)"],
        ["Newton inexacto, sin globalizar", "1 de 7", "14 %", "52", "10 631"],
        ["Búsqueda de línea (Armijo)", "5 de 7", "71 %", "88", "3 114"],
        ["Región de confianza (dogleg)", "3 de 7", "43 %", "54", "3 174"],
        ["Continuación pseudo-transitoria", "7 de 7", "100 %", "90", "863"],
    ], anchos=[3400, 1700, 900, 1700, 1660], sz=17)
    doc.pie_figura("**Tabla 8.** Burgers 1D, ν = 0.01, N = 200. El resultado central "
                   "del trabajo. Las medianas se calculan solo sobre los casos que "
                   "efectivamente convergen.")
    doc.texto(
        "El contraste es categórico. Newton sin globalizar resuelve uno de siete casos; "
        "la continuación pseudo-transitoria, los siete. Y —esto es lo que suele "
        "sorprender— **Ψtc no solo es la más robusta sino también la más barata**, con "
        "un factor de doce respecto de no globalizar. La robustez no se paga: en este "
        "problema se cobra. La razón es que las corridas sin globalización consumen "
        "enormes cantidades de trabajo divergiendo antes de agotar su presupuesto.")
    figura(doc, "fig3b_exito.png",
           "**Figura 7.** Robustez y costo son dos preguntas distintas. Izquierda: "
           "porcentaje de puntos iniciales desde los que cada estrategia converge. "
           "Derecha: trabajo mediano sobre los casos resueltos.", ancho_cm=15.5)
    doc.texto(
        "La Figura 8 muestra la garantía teórica de manera directa. Para __u₀ ≡ 0__, un "
        "punto inicial perfectamente razonable, el residuo de Newton sin globalizar "
        "**sube** de 1.4×10⁻² a 3.9×10¹: casi tres órdenes de magnitud peor que el "
        "punto de partida. Las otras tres curvas no suben nunca, y eso no es fortuna: "
        "es exactamente lo que imponen la condición de Armijo y su equivalente "
        "__ρ > 0__ en la región de confianza.")
    figura(doc, "fig3c_residuales.png",
           "**Figura 8.** Historias de residuo desde u₀ ≡ 0, por iteración y por "
           "trabajo. La monotonía de las curvas globalizadas es la garantía del "
           "teorema, escrita en el código.", ancho_cm=15.5)
    doc.texto(
        "Por último, se cuantificó la sensibilidad de Ψtc a su parámetro __δ₀__. La "
        "Tabla 9 muestra que el método falla tanto por defecto como por exceso, lo que "
        "confirma que se trata de un parámetro real y no de un detalle de "
        "implementación.")
    doc.tabla([
        ["δ₀", "Resultado", "Iters. de Newton", "Productos J·v"],
        ["10⁻²", "no converge", "250", "2 499"],
        ["1", "no converge", "250", "2 772"],
        ["70.7 = 1/‖F₀‖ (por defecto)", "converge", "141", "863"],
        ["10⁴", "converge", "8", "200"],
        ["10⁶", "converge", "9", "315"],
        ["10⁹", "no converge", "134", "15 586"],
    ], anchos=[3400, 2000, 2000, 1960], sz=17)
    doc.pie_figura("**Tabla 9.** Sensibilidad de Ψtc a δ₀ sobre Burgers 1D desde "
                   "u₀ ≡ 0. Demasiado pequeño sobre-amortigua; demasiado grande "
                   "equivale a Newton puro.")

    # ---------------------------------------------------------------
    doc.add(salto_pagina())
    doc.titulo_sub("5.4. Los límites: lo que la globalización no garantiza")
    doc.texto(
        "El enunciado exacto del teorema de convergencia global es más modesto de lo "
        "que su nombre sugiere. Si __F__ es continuamente diferenciable y los iterados "
        "permanecen en un conjunto acotado, la globalización garantiza que la sucesión "
        "converge a un **punto estacionario de f = ½‖F‖²**. Y un punto estacionario de "
        "__f__ no tiene por qué ser una raíz de __F__: basta con que "
        "__∇f = JᵀF = 0__ con __F ≠ 0__, lo que ocurre siempre que __J__ sea singular "
        "allí.")
    doc.texto("Los experimentos exhibieron los dos modos de fallo:")
    doc.vinieta(
        "**Mínimo local de la función de mérito.** La función de Freudenstein y Roth "
        "tiene uno en __x₂ ≈ −0.8968__, con __f ≈ 49__. Una sucesión que solo sabe bajar "
        "cae en ese valle y se detiene sin haber resuelto nada.")
    doc.vinieta(
        "**Jacobiano singular.** En el sistema 2×2 de la Sección 5.1, __J__ es singular "
        "sobre toda la recta __x₂ = 0__, porque su segunda columna se anula. Si el "
        "descenso arrastra la iteración hacia esa recta, el paso de Newton se vuelve "
        "enorme y mal condicionado, __λ__ tiende a cero y la búsqueda de línea se rinde.")
    doc.texto(
        "La consecuencia práctica es incómoda pero hay que decirla: **a veces Newton "
        "puro gana**, precisamente porque sus pasos salvajes saltan por encima del valle "
        "donde el método prudente queda atrapado. Sobre Freudenstein-Roth, Newton puro "
        "alcanza la raíz desde el 83 % de los puntos iniciales y Newton con Armijo solo "
        "desde el 37 %.")
    figura(doc, "clase1_cuencas.png",
           "**Figura 9.** Cuencas de convergencia. Arriba, un problema con exponencial "
           "empinada: la búsqueda de línea gana (100 % frente a 98 %). Abajo, "
           "Freudenstein-Roth: la búsqueda de línea pierde (37 % frente a 83 %), "
           "atrapada en el mínimo local de f. En rojo oscuro, los puntos desde los que "
           "no se alcanza ninguna raíz.", ancho_cm=13.5)
    doc.texto(
        "Ninguna estrategia de globalización puede evitar estos casos, porque son una "
        "limitación de la **formulación** —traducir un problema de raíces en uno de "
        "minimización— y no de los algoritmos. Saber exactamente qué se promete es lo "
        "que permite diagnosticar el fallo cuando ocurre en lugar de culpar al solver.")
    doc.texto(
        "Un último resultado, colateral pero instructivo, apareció al comparar la "
        "solución numérica de Burgers con la solución analítica. Pese a converger hasta "
        "__‖F‖ ≈ 6×10⁻¹⁰__, la diferencia máxima con la solución exacta muestreada es "
        "de 0.3, y el residuo de esa solución exacta vale 2.6×10⁻⁵ — es decir, **no es "
        "solución del sistema discreto**. La posición de la capa límite está "
        "exponencialmente mal determinada: el Jacobiano tiene un autovalor del orden de "
        "__e^(−1/ν)__ y desplazar la capa casi no altera el residuo. Es un recordatorio "
        "de que un residuo pequeño no equivale a una solución precisa: lo que acota el "
        "error es __‖J⁻¹‖·‖F‖__.")


def conclusiones(doc):
    doc.add(salto_pagina())
    doc.titulo_seccion("6. CONCLUSIONES")
    doc.texto(
        "**1. El método de Newton es local, y el teorema no dice cuán local.** Su "
        "convergencia cuadrática es real, pero solo dentro de una bola cuyo radio el "
        "teorema no permite estimar. Se comprobó con arctan(x) = 0, donde el umbral es "
        "|x₀| ≈ 1.3917.")
    doc.texto(
        "**2. El término de forzado decide el costo, no la velocidad.** Resolver el "
        "sistema lineal con más precisión de la necesaria multiplicó por 2.6 el trabajo "
        "total sin ganar precisión en la respuesta. Eisenstat-Walker resuelve el "
        "problema automáticamente y fue la opción más barata en todos los experimentos.")
    doc.texto(
        "**3. El enfoque matriz-libre es lo que hace posibles los problemas grandes**, "
        "y su precio es un Jacobiano con unas ocho cifras significativas. Ese ruido "
        "tiene consecuencias concretas: inutiliza el CG de Steihaug sobre las "
        "ecuaciones normales y fija un piso al residuo alcanzable.")
    doc.texto(
        "**4. El precondicionamiento vale más que la elección del algoritmo.** Un "
        "precondicionador laplaciano tridiagonal, factorizado una sola vez, redujo el "
        "trabajo en un factor de varios cientos. Ninguna otra decisión del método tuvo "
        "un efecto comparable.")
    doc.texto(
        "**5. Globalizar convierte un método frágil en uno utilizable.** Sobre la "
        "ecuación de Burgers, Newton sin globalizar resolvió 1 de 7 casos y la "
        "continuación pseudo-transitoria los 7, consumiendo además doce veces menos "
        "trabajo. La robustez no se pagó: se cobró.")
    doc.texto(
        "**6. Pero la garantía es más modesta de lo que el nombre sugiere.** La "
        "globalización asegura que ‖F‖ no aumente y que se converja a un punto "
        "estacionario de la función de mérito, que puede no ser raíz. Se documentaron "
        "los dos modos de fallo, incluido el caso en que Newton puro supera al método "
        "globalizado.")
    doc.texto(
        "**7. La métrica importa.** Contar iteraciones de Newton en lugar de "
        "evaluaciones de F llevó, en dos de los tres experimentos, a conclusiones "
        "exactamente invertidas.")


def referencias(doc):
    doc.titulo_seccion("7. REFERENCIAS")
    refs = [
        "Dembo, R. S., Eisenstat, S. C. & Steihaug, T. (1982). «Inexact Newton "
        "Methods». __SIAM Journal on Numerical Analysis__, 19(2), 400–408.",
        "Dennis, J. E. & Schnabel, R. B. (1996). __Numerical Methods for "
        "Unconstrained Optimization and Nonlinear Equations__. SIAM.",
        "Eisenstat, S. C. & Walker, H. F. (1996). «Choosing the forcing terms in an "
        "inexact Newton method». __SIAM Journal on Scientific Computing__, 17(1), 16–32.",
        "Kelley, C. T. (1995). __Iterative Methods for Linear and Nonlinear "
        "Equations__. SIAM, Frontiers in Applied Mathematics 16.",
        "Kelley, C. T. (2003). __Solving Nonlinear Equations with Newton's Method__. "
        "SIAM, Fundamentals of Algorithms 1.",
        "Kelley, C. T. & Keyes, D. E. (1998). «Convergence analysis of pseudo-transient "
        "continuation». __SIAM Journal on Numerical Analysis__, 35(2), 508–523.",
        "Knoll, D. A. & Keyes, D. E. (2004). «Jacobian-free Newton-Krylov methods: a "
        "survey of approaches and applications». __Journal of Computational Physics__, "
        "193(2), 357–397.",
        "Nocedal, J. & Wright, S. J. (2006). __Numerical Optimization__, 2ª ed. Springer.",
        "Pawlowski, R. P., Shadid, J. N., Simonis, J. P. & Walker, H. F. (2006). "
        "«Globalization techniques for Newton-Krylov methods and applications to the "
        "fully coupled solution of the Navier-Stokes equations». __SIAM Review__, "
        "48(4), 700–721.",
        "Saad, Y. (2003). __Iterative Methods for Sparse Linear Systems__, 2ª ed. SIAM.",
        "Virtanen, P. et al. (2020). «SciPy 1.0: fundamental algorithms for scientific "
        "computing in Python». __Nature Methods__, 17, 261–272.",
    ]
    for r in refs:
        doc.add(parrafo(runs_con_marcas("[" + str(refs.index(r) + 1) + "] " + r, sz=20),
                        jc="both", despues=120, sangria=360))


def anexo(doc):
    doc.add(salto_pagina())
    doc.titulo_seccion("ANEXO A. CÓMO EJECUTAR EL CÓDIGO")
    doc.texto("Requisitos: Python 3.9 o superior con NumPy, SciPy y Matplotlib.")
    doc.texto("`pip install numpy scipy matplotlib`", jc="left")
    doc.titulo_sub("A.1. Ejercicios de la exposición")
    doc.texto("`cd ejercicios_exposicion`", jc="left", despues=40)
    doc.texto("`python3 ej1_newton_vs_globalizado.py`      # ~20 s", jc="left", despues=40)
    doc.texto("`python3 ej2_bratu1d_newton_krylov.py`      # ~20 s", jc="left", despues=40)
    doc.texto("`python3 ej3_comparativa_globalizacion.py`  # ~40 s", jc="left")
    doc.texto(
        "Cada programa imprime tablas por consola y deja sus figuras en "
        "`ejercicios_exposicion/figuras/`. Todas las cifras de la Sección 5 salen de "
        "esas corridas.")
    doc.titulo_sub("A.2. Ejercicios para la clase")
    doc.texto("`cd ejercicios_clase`", jc="left", despues=40)
    doc.texto("`python3 clase1_armijo.py`           # ~5 s", jc="left", despues=40)
    doc.texto("`python3 clase2_matrix_free.py`      # ~25 s", jc="left", despues=40)
    doc.texto("`python3 clase3_bratu_forcing.py`    # ~35 s", jc="left")
    doc.texto(
        "Son autónomos: cada uno se puede copiar suelto y correr sin el resto del "
        "proyecto. Terminan con un bloque de preguntas de análisis.")
    doc.titulo_sub("A.3. Regenerar la presentación y este informe")
    doc.texto("`python3 presentacion/build_presentacion.py`", jc="left", despues=40)
    doc.texto("`python3 informe/generar_informe.py`", jc="left")


# =============================================================================
def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    doc = Documento(PLANTILLA)
    portada(doc)
    introduccion(doc)
    marco_teorico(doc)
    globalizacion(doc)
    implementacion(doc)
    resultados(doc)
    conclusiones(doc)
    referencias(doc)
    anexo(doc)
    ruta = doc.guardar(SALIDA)
    print(f"→ {ruta}  ({ruta.stat().st_size/1024:,.0f} KB, "
          f"{len(doc.medios)} imágenes empotradas)")


if __name__ == "__main__":
    main()
