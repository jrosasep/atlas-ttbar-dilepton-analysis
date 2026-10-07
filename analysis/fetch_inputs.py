"""Descargar archivos públicos y comprobar su integridad.

Los ROOT proceden de CERN; este programa no genera eventos.
analysis.py importa verify. La descarga sólo comienza con la acción download.
"""
from __future__ import annotations
import argparse                    # Elegir la acción al ejecutar el programa.
import concurrent.futures          # Coordinar transferencias simultáneas.
import hashlib                     # Huella SHA256 del archivo completo.
from html.parser import HTMLParser # Extraer documentación HTML.
import json                        # Leer y guardar registros.
from pathlib import Path           # Ubicar archivos.
import time                        # Pausas entre intentos fallidos.
import urllib.request              # Solicitudes HTTPS con verificación TLS.
import zlib                        # Checksum Adler-32 de CERN.

ruta_del_codigo = Path(__file__)
ruta_completa_del_codigo = ruta_del_codigo.resolve()
BASE = ruta_completa_del_codigo.parent
OFFICIAL_REPOSITORY = 'atlas-outreach-data-tools/atlas-outreach-cpp-framework-13tev'
OFFICIAL_REVISION = 'ff71d6ba82f2afd5a45d2e9b7f80bf915ceb1c84'


class TextParser(HTMLParser):
    """Reunir texto y enlaces mientras HTMLParser recorre una página."""
    def __init__(self):
        super().__init__()
        self.parts = []
        self.links = []
    def handle_starttag(self, tag, attrs):
        if tag in ('p', 'tr', 'h1', 'h2', 'h3', 'li', 'br'):
            self.parts.append('\n')
        if tag == 'a':
            for nombre, valor in attrs:
                if nombre == 'href':
                    self.links.append(valor)
    def handle_data(self, data):
        self.parts.append(data)


def fetch(url):
    """Leer los bytes de una URL pública, con hasta tres intentos."""
    for intento in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as respuesta:
                return respuesta.read()
        except Exception:
            if intento == 2:
                raise
            time.sleep(2)


def catalog():
    """Guardar los catálogos y la documentación oficiales completos."""
    carpeta = BASE / 'provenance'
    carpeta.mkdir(parents=True, exist_ok=True)
    for numero in (93913, 93934):
        contenido = fetch(f'https://opendata.cern.ch/api/records/{numero}')
        ruta = carpeta / f'catalog_{numero}.json'
        ruta.write_bytes(contenido)
        catalogo = json.loads(contenido)
        metadata = catalogo['metadata']
        print(numero, metadata['title'], metadata['distribution'], flush=True)
    for nombre in ('13TeV25_details', '13TeV25_metadata'):
        url = f'https://opendata.atlas.cern/docs/data/for_education/{nombre}'
        contenido = fetch(url)
        (carpeta / f'{nombre}.html').write_bytes(contenido)
        parser = TextParser()
        parser.feed(contenido.decode())
        texto = ''.join(parser.parts)
        (carpeta / f'{nombre}.txt').write_text(texto, encoding='utf-8')
        enlaces = json.dumps(parser.links, indent=2)
        (carpeta / f'{nombre}_links.json').write_text(enlaces)


def references():
    """Archivar metadatos y la referencia educativa citada en la Pauta 3.

    Un commit identifica una revisión; su árbol identifica los archivos.
    La descarga de contenido usa el commit, no la huella del árbol.
    """
    carpeta = BASE / 'provenance'
    carpeta.mkdir(parents=True, exist_ok=True)
    contenido_csv = fetch('https://opendata.atlas.cern/files/metadata.csv')
    (carpeta / 'metadata.csv').write_bytes(contenido_csv)
    url_commit = f'https://api.github.com/repos/{OFFICIAL_REPOSITORY}/commits/{OFFICIAL_REVISION}'
    commit = json.loads(fetch(url_commit))
    revision = commit['sha']
    if revision != OFFICIAL_REVISION:
        raise ValueError('La referencia recibida no coincide con el commit fijado')
    sha_arbol = commit['commit']['tree']['sha']
    url_arbol = f'https://api.github.com/repos/{OFFICIAL_REPOSITORY}/git/trees/{sha_arbol}?recursive=1'
    arbol = json.loads(fetch(url_arbol))
    (carpeta / 'official_commit.json').write_text(json.dumps(commit, indent=2))
    (carpeta / 'official_tree.json').write_text(json.dumps(arbol, indent=2))
    for entrada in arbol['tree']:
        ruta_relativa = entrada['path']
        es_archivo = entrada['type'] == 'blob'
        es_analisis = 'TTbarDilepAnalysis' in ruta_relativa
        es_informacion = 'infofile' in ruta_relativa.lower()
        es_licencia = ruta_relativa == 'LICENSE'
        if es_archivo and (es_analisis or es_informacion or es_licencia):
            destino = carpeta / 'official' / ruta_relativa
            destino.parent.mkdir(parents=True, exist_ok=True)
            url = f'https://raw.githubusercontent.com/{OFFICIAL_REPOSITORY}/{revision}/{ruta_relativa}'
            destino.write_bytes(fetch(url))
            print(ruta_relativa, flush=True)


def verify(path, info):
    """Comprobar tamaño y Adler-32 y devolver la huella local SHA256.

    Devuelve None si el archivo falta o está incompleto. Un checksum diferente
    detiene el procesamiento. Lee bytes sin interpretar los eventos físicos.
    """
    if not path.exists():
        return None
    tamaño_actual = path.stat().st_size
    if tamaño_actual != info['size']:
        return None
    adler = 1
    sha = hashlib.sha256()
    with path.open('rb') as archivo:
        while True:
            bloque = archivo.read(8 * 1024 * 1024)
            if not bloque:
                break
            adler = zlib.adler32(bloque, adler)
            sha.update(bloque)
    adler_32bits = adler & 0xffffffff
    checksum = f'adler32:{adler_32bits:08x}'
    if checksum != info['checksum']:
        raise ValueError(f'Checksum incorrecto: {path.name}: {checksum}')
    return sha.hexdigest()


def ranged_download(url, partial, size):
    """Transferir rangos de bytes de un archivo grande, sin seleccionar eventos."""
    tamaño_bloque = 16 * 1024**2
    rangos = []
    for inicio in range(0, size, tamaño_bloque):
        final = min(inicio + tamaño_bloque, size) - 1
        rangos.append((inicio, final))

    def read_range(pair):
        """Comprobar que cada respuesta contenga exactamente el rango pedido."""
        inicio, final = pair
        for intento in range(3):
            try:
                solicitud = urllib.request.Request(
                    url, headers={'Range': f'bytes={inicio}-{final}'}
                )
                with urllib.request.urlopen(solicitud, timeout=180) as respuesta:
                    esperado = f'bytes {inicio}-{final}/{size}'
                    estado_correcto = respuesta.status == 206
                    rango_correcto = respuesta.headers.get('Content-Range') == esperado
                    if not estado_correcto or not rango_correcto:
                        raise ValueError('El servidor no respetó el rango solicitado')
                    datos = respuesta.read()
                    if len(datos) != final - inicio + 1:
                        raise ValueError('Rango truncado')
                return inicio, datos
            except Exception:
                if intento == 2:
                    raise
                time.sleep(2)

    with partial.open('wb') as archivo:
        archivo.truncate(size)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            # Cada Future representa una transferencia pendiente.
            # Sólo se conservan hasta ocho rangos simultáneos en memoria.
            siguientes = iter(rangos)
            pendientes = set()
            for _ in range(min(8, len(rangos))):
                rango = next(siguientes)
                pendientes.add(pool.submit(read_range, rango))
            while pendientes:
                terminadas, pendientes = concurrent.futures.wait(
                    pendientes, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for tarea in terminadas:
                    posicion, datos = tarea.result()
                    archivo.seek(posicion)
                    archivo.write(datos)
                    siguiente = next(siguientes, None)
                    if siguiente is not None:
                        pendientes.add(pool.submit(read_range, siguiente))


def download_one(info):
    """Reutilizar un ROOT íntegro o descargarlo primero a un archivo .part."""
    destino = BASE / 'data' / info['key']
    destino.parent.mkdir(exist_ok=True)
    huella = verify(destino, info)
    url = 'https://opendata.cern.ch/eos/opendata/atlas/rucio/opendata/' + info['key']
    if huella is None:
        parcial = destino.with_suffix('.root.part')
        for intento in range(3):
            try:
                if info['size'] > 64 * 1024**2:
                    ranged_download(url, parcial, info['size'])
                else:
                    with urllib.request.urlopen(url, timeout=180) as respuesta:
                        with parcial.open('wb') as archivo:
                            while True:
                                bloque = respuesta.read(4 * 1024 * 1024)
                                if not bloque:
                                    break
                                archivo.write(bloque)
                huella = verify(parcial, info)
                if huella is None:
                    raise ValueError(f'Descarga incompleta: {info["key"]}')
                parcial.replace(destino)
                break
            except Exception:
                if intento == 2:
                    raise
                time.sleep(3)
    tamaño_mb = info['size'] / 1e6
    print(f'VERIFICADO {info["key"]} ({tamaño_mb:.1f} MB)', flush=True)
    registro = dict(info)
    registro['https_access_url'] = url
    registro['sha256'] = huella
    registro['local_file'] = str(destino.relative_to(BASE))
    return registro


def download():
    """Descargar sólo la lista del manifiesto y guardar sus huellas."""
    texto = (BASE / 'manifest.json').read_text()
    manifiesto = json.loads(texto)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        tareas = pool.map(download_one, manifiesto['files'])
        archivos_completados = list(tareas)
    verificado = dict(manifiesto)
    verificado['files'] = archivos_completados
    carpeta = BASE / 'provenance'
    carpeta.mkdir(parents=True, exist_ok=True)
    ruta_verificada = carpeta / 'verified_inputs.json'
    ruta_verificada.write_text(json.dumps(verificado, indent=2), encoding='utf-8')


def main():
    """Elegir una acción explícita; importar este módulo no descarga archivos."""
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('catalog', 'download', 'references'))
    args = parser.parse_args()
    if args.action == 'catalog':
        catalog()
    elif args.action == 'download':
        download()
    elif args.action == 'references':
        references()


if __name__ == '__main__':
    main()
