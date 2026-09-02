# -*- coding: utf-8 -*-
"""
docx_min.py — Constructor mínimo de documentos .docx sobre una plantilla.

Por qué existe: `python-docx` no está instalado en esta máquina y no queremos
depender de que lo esté. Un .docx es un ZIP con XML adentro, así que se puede
construir con la biblioteca estándar. Lo único que hace este módulo es:

  1. abrir la plantilla Gnombres_Dat252.docx como ZIP,
  2. reemplazar word/document.xml por el que armamos aquí,
  3. agregar las imágenes en word/media/ y registrarlas en las relaciones,
  4. dejar intactos styles.xml, theme1.xml, settings.xml y el sectPr,
     que son los que definen el formato de la plantilla.

Así el documento generado hereda EXACTAMENTE la tipografía, los márgenes y el
tamaño de página de la plantilla institucional.
"""

from __future__ import annotations

import io
import pathlib
import re
import shutil
import zipfile

EMU_POR_PULGADA = 914400
EMU_POR_CM = 360000

NS_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
NS_PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
TIPO_IMAGEN = ("http://schemas.openxmlformats.org/officeDocument/2006/"
               "relationships/image")


def esc(t):
    """Escapa texto para XML."""
    return (str(t).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


# =============================================================================
# Runs y párrafos
# =============================================================================
def run(texto, negrita=False, cursiva=False, sz=None, color="000000",
        mono=False, sub=False, sup=False):
    """Un fragmento de texto con formato. sz va en medios puntos (24 = 12 pt)."""
    rpr = []
    if negrita:
        rpr.append("<w:b/>")
    if cursiva:
        rpr.append("<w:i/>")
    if mono:
        rpr.append('<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" '
                   'w:cs="Consolas"/>')
    if sz:
        rpr.append(f'<w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/>')
    rpr.append(f'<w:color w:val="{color}"/>')
    if sub:
        rpr.append('<w:vertAlign w:val="subscript"/>')
    if sup:
        rpr.append('<w:vertAlign w:val="superscript"/>')
    rpr.append('<w:lang w:val="es-BO"/>')
    return (f'<w:r><w:rPr>{"".join(rpr)}</w:rPr>'
            f'<w:t xml:space="preserve">{esc(texto)}</w:t></w:r>')


_MARCAS = re.compile(r"(\*\*.+?\*\*|__.+?__|`[^`]+?`|~[^~]+?~|\^[^\^]+?\^)")


def runs_con_marcas(texto, sz=None, color="000000"):
    """Convierte un mini-marcado en runs.

        **negrita**   __cursiva__   `monoespaciado`   ~subíndice~   ^superíndice^
    """
    salida = []
    for trozo in _MARCAS.split(texto):
        if not trozo:
            continue
        if "\n" in trozo and not trozo.startswith(("**", "__", "`", "~", "^")):
            # Un salto de línea dentro de un párrafo o de una celda de tabla
            # se traduce a <w:br/>, no a un carácter literal.
            piezas = trozo.split("\n")
            for i, pieza in enumerate(piezas):
                if i:
                    salida.append("<w:r><w:br/></w:r>")
                if pieza:
                    salida.append(run(pieza, sz=sz, color=color))
            continue
        if trozo.startswith("**") and trozo.endswith("**"):
            salida.append(run(trozo[2:-2], negrita=True, sz=sz, color=color))
        elif trozo.startswith("__") and trozo.endswith("__"):
            salida.append(run(trozo[2:-2], cursiva=True, sz=sz, color=color))
        elif trozo.startswith("`") and trozo.endswith("`"):
            salida.append(run(trozo[1:-1], mono=True,
                              sz=(sz - 2) if sz else 20, color="1F4E79"))
        elif trozo.startswith("~") and trozo.endswith("~"):
            salida.append(run(trozo[1:-1], sub=True, sz=sz, color=color))
        elif trozo.startswith("^") and trozo.endswith("^"):
            salida.append(run(trozo[1:-1], sup=True, sz=sz, color=color))
        else:
            salida.append(run(trozo, sz=sz, color=color))
    return "".join(salida)


def parrafo(contenido_runs, jc=None, antes=0, despues=200, linea=276,
            borde_inf=None, sangria=0, mantener=False):
    """Un párrafo. `contenido_runs` ya viene como XML de runs."""
    ppr = []
    if mantener:
        ppr.append("<w:keepNext/><w:keepLines/>")
    if borde_inf:
        ppr.append(f'<w:pBdr><w:bottom w:val="single" w:sz="{borde_inf[1]}" '
                   f'w:space="4" w:color="{borde_inf[0]}"/></w:pBdr>')
    ppr.append(f'<w:spacing w:before="{antes}" w:after="{despues}"'
               + (f' w:line="{linea}" w:lineRule="auto"' if linea else "") + "/>")
    if sangria:
        ppr.append(f'<w:ind w:left="{sangria}"/>')
    if jc:
        ppr.append(f'<w:jc w:val="{jc}"/>')
    return f"<w:p><w:pPr>{''.join(ppr)}</w:pPr>{contenido_runs}</w:p>"


def salto_pagina():
    return ('<w:p><w:pPr><w:spacing w:after="0"/></w:pPr>'
            '<w:r><w:br w:type="page"/></w:r></w:p>')


# =============================================================================
# Documento
# =============================================================================
class Documento:
    """Acumula cuerpo XML y las imágenes, y al final reescribe el ZIP."""

    def __init__(self, plantilla):
        self.plantilla = pathlib.Path(plantilla)
        self.cuerpo = []
        self.medios = []            # (nombre_en_zip, bytes)
        self._sig_rel = 100
        self._sig_doc = 1000

    # ---------------- bloques de alto nivel ----------------
    def add(self, xml):
        self.cuerpo.append(xml)

    def titulo_seccion(self, texto):
        """Formato de la plantilla: negrita, 12 pt (sz 24)."""
        self.add(parrafo(run(texto, negrita=True, sz=24),
                         antes=240, despues=120, mantener=True))

    def titulo_sub(self, texto):
        """Subsección: negrita, 11 pt (sz 22)."""
        self.add(parrafo(run(texto, negrita=True, sz=22),
                         antes=160, despues=100, mantener=True))

    def texto(self, contenido, jc="both", sz=22, despues=160, sangria=0):
        self.add(parrafo(runs_con_marcas(contenido, sz=sz), jc=jc,
                         despues=despues, sangria=sangria))

    def vinieta(self, contenido, sz=22):
        self.add(parrafo(runs_con_marcas("•\t" + contenido, sz=sz), jc="both",
                         despues=80, sangria=360))

    def espacio(self, alto=120):
        self.add(parrafo("", despues=alto))

    def regla(self, color="1F4E79", grosor=12, antes=0, despues=120):
        self.add(parrafo("", antes=antes, despues=despues,
                         borde_inf=(color, grosor)))

    # ---------------- imágenes ----------------
    def imagen(self, ruta, ancho_cm=None, ancho_emu=None, jc="center",
               despues=120):
        """Inserta una imagen PNG en línea, centrada."""
        ruta = pathlib.Path(ruta)
        datos = ruta.read_bytes()
        px_w, px_h = _tamano_png(datos)

        if ancho_emu is None:
            ancho_emu = int((ancho_cm or 15.0) * EMU_POR_CM)
        alto_emu = int(ancho_emu * px_h / px_w)

        rid = f"rIdImg{self._sig_rel}"
        self._sig_rel += 1
        nombre = f"media/img{self._sig_rel}{ruta.suffix.lower()}"
        self.medios.append((nombre, datos, rid))

        did = self._sig_doc
        self._sig_doc += 1
        dibujo = (
            f'<w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0">'
            f'<wp:extent cx="{ancho_emu}" cy="{alto_emu}"/>'
            f'<wp:effectExtent l="0" t="0" r="0" b="0"/>'
            f'<wp:docPr id="{did}" name="Imagen {did}"/>'
            f'<wp:cNvGraphicFramePr><a:graphicFrameLocks xmlns:a="{NS_A}" '
            f'noChangeAspect="1"/></wp:cNvGraphicFramePr>'
            f'<a:graphic xmlns:a="{NS_A}"><a:graphicData '
            f'uri="{NS_PIC}"><pic:pic xmlns:pic="{NS_PIC}">'
            f'<pic:nvPicPr><pic:cNvPr id="{did}" name="{esc(ruta.name)}"/>'
            f'<pic:cNvPicPr/></pic:nvPicPr>'
            f'<pic:blipFill><a:blip r:embed="{rid}"/>'
            f'<a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
            f'<pic:spPr><a:xfrm><a:off x="0" y="0"/>'
            f'<a:ext cx="{ancho_emu}" cy="{alto_emu}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
            f'</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r>')
        self.add(parrafo(dibujo, jc=jc, despues=despues))

    def pie_figura(self, texto):
        self.add(parrafo(runs_con_marcas(texto, sz=18, color="404040"),
                         jc="center", despues=220))

    # ---------------- tablas ----------------
    def tabla(self, filas, anchos, encabezado=True, sz=18,
              relleno_enc="1F4E79", relleno_alt="F2F6FA"):
        """Tabla con bordes, encabezado azul y filas alternadas."""
        grid = "".join(f'<w:gridCol w:w="{a}"/>' for a in anchos)
        xml = [
            '<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/>'
            '<w:tblBorders>'
            '<w:top w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '<w:left w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '<w:bottom w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '<w:right w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="BFC9D4"/>'
            '</w:tblBorders><w:tblLayout w:type="fixed"/>'
            '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" '
            'w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>'
            f'</w:tblPr><w:tblGrid>{grid}</w:tblGrid>']

        for i, fila in enumerate(filas):
            es_enc = encabezado and i == 0
            relleno = relleno_enc if es_enc else (
                relleno_alt if (i % 2 == 0) else "FFFFFF")
            color = "FFFFFF" if es_enc else "000000"
            celdas = []
            for j, celda in enumerate(fila):
                jc = "left" if j == 0 else "center"
                celdas.append(
                    f'<w:tc><w:tcPr><w:tcW w:w="{anchos[j]}" w:type="dxa"/>'
                    f'<w:shd w:val="clear" w:color="auto" w:fill="{relleno}"/>'
                    f'<w:tcMar><w:top w:w="60" w:type="dxa"/>'
                    f'<w:left w:w="100" w:type="dxa"/>'
                    f'<w:bottom w:w="60" w:type="dxa"/>'
                    f'<w:right w:w="100" w:type="dxa"/></w:tcMar>'
                    f'<w:vAlign w:val="center"/></w:tcPr>'
                    + parrafo(runs_con_marcas(str(celda), sz=sz, color=color)
                              if not es_enc else
                              run(str(celda), negrita=True, sz=sz, color=color),
                              jc=jc, despues=0, linea=240)
                    + '</w:tc>')
            trpr = '<w:trPr><w:tblHeader/><w:cantSplit/></w:trPr>' if es_enc else ''
            xml.append(f"<w:tr>{trpr}{''.join(celdas)}</w:tr>")
        xml.append("</w:tbl>")
        self.add("".join(xml))
        self.add(parrafo("", despues=160))

    # ---------------- escritura ----------------
    def guardar(self, destino, sect_pr=None):
        destino = pathlib.Path(destino)
        original = zipfile.ZipFile(self.plantilla)
        partes = {n: original.read(n) for n in original.namelist()}
        original.close()

        if sect_pr is None:
            doc_viejo = partes["word/document.xml"].decode("utf-8")
            sect_pr = re.search(r"<w:sectPr.*?</w:sectPr>", doc_viejo,
                                re.S).group(0)

        encabezado_ns = re.search(
            r"<w:document[^>]*>", partes["word/document.xml"].decode("utf-8")).group(0)

        documento = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                     + encabezado_ns + "<w:body>"
                     + "".join(self.cuerpo) + sect_pr + "</w:body></w:document>")
        partes["word/document.xml"] = documento.encode("utf-8")

        # --- relaciones ---
        rels = partes["word/_rels/document.xml.rels"].decode("utf-8")
        nuevas = "".join(
            f'<Relationship Id="{rid}" Type="{TIPO_IMAGEN}" Target="{nombre}"/>'
            for nombre, _datos, rid in self.medios)
        partes["word/_rels/document.xml.rels"] = rels.replace(
            "</Relationships>", nuevas + "</Relationships>").encode("utf-8")

        # --- tipos de contenido: registrar png ---
        ct = partes["[Content_Types].xml"].decode("utf-8")
        if 'Extension="png"' not in ct:
            ct = ct.replace(
                '<Default Extension="xml"',
                '<Default Extension="png" ContentType="image/png"/>'
                '<Default Extension="xml"')
        partes["[Content_Types].xml"] = ct.encode("utf-8")

        # --- imágenes ---
        for nombre, datos, _rid in self.medios:
            partes["word/" + nombre] = datos

        destino.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(destino, "w", zipfile.ZIP_DEFLATED) as z:
            # [Content_Types].xml debe ir primero en un OPC bien formado
            z.writestr("[Content_Types].xml", partes.pop("[Content_Types].xml"))
            for nombre, datos in partes.items():
                z.writestr(nombre, datos)
        return destino


def _tamano_png(datos):
    """Ancho y alto en píxeles leyendo la cabecera IHDR del PNG."""
    if datos[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("Solo se admiten archivos PNG")
    ancho = int.from_bytes(datos[16:20], "big")
    alto = int.from_bytes(datos[20:24], "big")
    return ancho, alto
