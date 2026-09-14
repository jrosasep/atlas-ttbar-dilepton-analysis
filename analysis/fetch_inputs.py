"""Descarga insumos públicos y comprueba su integridad antes de analizarlos.

No se generan colisiones artificiales: los archivos vienen del catálogo CERN.
El manifiesto conserva URL, tamaño y checksum oficial de cada muestra.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import time
import urllib.request
import zlib

BASE = Path(__file__).resolve().parent


class TextParser(HTMLParser):
    """Extrae texto legible y enlaces de la documentación técnica oficial."""
    def __init__(self):
        super().__init__()
        self.parts, self.links = [], []
    def handle_starttag(self, tag, attrs):
        if tag in ('p', 'tr', 'h1', 'h2', 'h3', 'li', 'br'):
            self.parts.append('\n')
        if tag == 'a':
            self.links.extend(v for k, v in attrs if k == 'href')
    def handle_data(self, data):
        self.parts.append(data)


def fetch(url):
    """Lee una URL pública, con tiempo límite y un máximo de tres intentos."""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as response:
                return response.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)


def catalog():
    """Guarda el catálogo íntegro; no deduce nombres ni tamaños de archivos."""
    dest = BASE / 'provenance'
    dest.mkdir(parents=True, exist_ok=True)
    for rec in (93913, 93934):
        content = fetch(f'https://opendata.cern.ch/api/records/{rec}')
        (dest / f'catalog_{rec}.json').write_bytes(content)
        meta = json.loads(content)['metadata']
        print(rec, meta['title'], meta['distribution'], flush=True)
    for name in ('13TeV25_details', '13TeV25_metadata'):
        url = f'https://opendata.atlas.cern/docs/data/for_education/{name}'
        content = fetch(url)
        (dest / f'{name}.html').write_bytes(content)
        parser = TextParser()
        parser.feed(content.decode())
        (dest / f'{name}.txt').write_text(''.join(parser.parts), encoding='utf-8')
        (dest / f'{name}_links.json').write_text(json.dumps(parser.links, indent=2))


def references():
    """Archiva metadatos y código oficial con el commit de procedencia."""
    dest = BASE / 'provenance'
    (dest / 'metadata.csv').write_bytes(fetch('https://opendata.atlas.cern/files/metadata.csv'))
    repo = 'atlas-outreach-data-tools/atlas-outreach-cpp-framework-13tev'
    tree = json.loads(fetch(f'https://api.github.com/repos/{repo}/git/trees/master?recursive=1'))
    (dest / 'official_tree.json').write_text(json.dumps(tree, indent=2))
    for entry in tree['tree']:
        path = entry['path']
        if entry['type'] == 'blob' and ('TTbarDilepAnalysis' in path or 'infofile' in path.lower() or path == 'LICENSE'):
            print(path, flush=True)
            target = dest / 'official' / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(fetch(f'https://raw.githubusercontent.com/{repo}/{tree["sha"]}/{path}'))


def verify(path, info):
    """Adler-32 verifica identidad con CERN; SHA-256 deja una huella local."""
    if not path.exists() or path.stat().st_size != info['size']:
        return None
    adler, sha = 1, hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            adler = zlib.adler32(block, adler)
            sha.update(block)
    actual = f'adler32:{adler & 0xffffffff:08x}'
    if actual != info['checksum']:
        raise ValueError(f'Checksum incorrecto: {path.name}: {actual}')
    return sha.hexdigest()


def download_one(info):
    """Descarga a .part y publica el archivo local solo tras verificarlo."""
    dest = BASE / 'data' / info['key']
    dest.parent.mkdir(exist_ok=True)
    digest = verify(dest, info)
    transport_url = 'https://opendata.cern.ch/eos/opendata/atlas/rucio/opendata/' + info['key']
    if digest is None:
        partial = dest.with_suffix('.root.part')
        for attempt in range(3):
            try:
                if info['size'] > 64 * 1024**2:
                    ranged_download(transport_url, partial, info['size'])
                else:
                    with urllib.request.urlopen(transport_url, timeout=180) as response, partial.open('wb') as out:
                        while block := response.read(4 * 1024 * 1024):
                            out.write(block)
                digest = verify(partial, info)
                if digest is None:
                    raise ValueError(f'Descarga incompleta: {info["key"]}')
                partial.replace(dest)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(3)
    print(f'VERIFICADO {info["key"]} ({info["size"]/1e6:.1f} MB)', flush=True)
    return dict(info, https_access_url=transport_url, sha256=digest, local_file=str(dest.relative_to(BASE)))


def ranged_download(url, partial, size):
    """Acceso HTTPS directo oficial: descarga rangos y comprueba sus offsets.

    Es el mapeo de protocolo publicado por atlasopenmagic. Se verifica TLS y no
    se abre el ROOT hasta comprobar tamaño y checksum oficial del archivo entero.
    Los rangos no cambian los eventos ni seleccionan una fracción de la muestra.
    """
    block_size = 16 * 1024**2
    ranges = [(start, min(start+block_size, size)-1) for start in range(0, size, block_size)]
    def read_range(pair):
        start, end = pair
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={'Range': f'bytes={start}-{end}'})
                with urllib.request.urlopen(req, timeout=180) as response:
                    expected = f'bytes {start}-{end}/{size}'
                    if response.status != 206 or response.headers.get('Content-Range') != expected:
                        raise ValueError('El servidor no respetó el rango solicitado')
                    data = response.read()
                    if len(data) != end-start+1:
                        raise ValueError('Rango truncado')
                return start, data
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2)
    with partial.open('wb') as stream:
        stream.truncate(size)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            # Como máximo ocho bloques en vuelo: no acumular todo el archivo
            # de varios GB dentro de los resultados de los Future.
            todo = iter(ranges)
            pending = {pool.submit(read_range, pair) for pair in [next(todo) for _ in range(min(8, len(ranges)))]}
            while pending:
                done, pending = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    offset, data = future.result()
                    stream.seek(offset)
                    stream.write(data)
                    pair = next(todo, None)
                    if pair is not None:
                        pending.add(pool.submit(read_range, pair))


def download():
    """Ejecuta únicamente la lista explícita y revisada del manifiesto."""
    manifest = json.loads((BASE / 'manifest.json').read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        completed = list(pool.map(download_one, manifest['files']))
    (BASE / 'provenance' / 'verified_inputs.json').write_text(
        json.dumps(dict(manifest, files=completed), indent=2), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('catalog', 'download', 'references'))
    args = parser.parse_args()
    globals()[args.action]()
