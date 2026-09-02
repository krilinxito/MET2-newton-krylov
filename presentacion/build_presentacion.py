# -*- coding: utf-8 -*-
"""
Construye /home/max1/MET2/presentacion.html a partir de esta plantilla y de las
figuras generadas por los ejercicios.

Las imágenes se empotran como data URI en base64 para que el archivo final sea
UNO SOLO y funcione sin conexión a internet: se abre con doble clic en el aula
aunque no haya red y aunque no esté la carpeta de figuras.

    python3 presentacion/build_presentacion.py
"""

import base64
import pathlib
import re
import sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent
FIGURAS = [RAIZ / "ejercicios_exposicion" / "figuras"]
SALIDA = RAIZ / "presentacion.html"
PLANTILLA = pathlib.Path(__file__).with_name("plantilla.html")


def data_uri(nombre):
    for carpeta in FIGURAS:
        ruta = carpeta / nombre
        if ruta.exists():
            b64 = base64.b64encode(ruta.read_bytes()).decode("ascii")
            return f"data:image/png;base64,{b64}"
    raise FileNotFoundError(
        f"No encuentro la figura {nombre!r}. Corré primero los notebooks de\n"
        f"    ejercicios_exposicion/  (ej1, ej2 y ej3), que las generan.")


def main():
    html = PLANTILLA.read_text(encoding="utf-8")
    faltantes = []
    for nombre in sorted(set(re.findall(r"%%FIG:([^%]+)%%", html))):
        try:
            html = html.replace(f"%%FIG:{nombre}%%", data_uri(nombre))
        except FileNotFoundError as e:
            faltantes.append(str(e))
    if faltantes:
        print("\n".join(faltantes))
        return 1
    SALIDA.write_text(html, encoding="utf-8")
    kb = SALIDA.stat().st_size / 1024
    n = len(re.findall(r'class="slide', html))
    print(f"→ {SALIDA}  ({kb:,.0f} KB, {n} diapositivas)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
