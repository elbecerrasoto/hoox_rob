#!/usr/bin/env python3
import os
import re
import sys
from collections import namedtuple
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Union

import numpy as np
from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERIES_DIR = Path(sys.argv[1])
GENOMES_FILE = sys.argv[2]
OUT_FILE = Path(sys.argv[3])

GENOME_REGEX = re.compile(r"(GC[FA]_\d+\.\d)\.faa$")
FIELDS = (
    "genome",
    "pid",
    "query",
    "score",
    "evalue",
    "start",
    "end",
    "pid_txt",
    "query_txt",
)
Results = namedtuple("Results", FIELDS)


class HMMFiles(Iterable[HMM]):
    def __init__(self, *files: Union[str, bytes, os.PathLike]):
        self.files = files

    def __iter__(self):
        for file in self.files:
            with HMMFile(file) as hmm_file:
                yield from hmm_file


def get_hmms(queries_path):
    queries_path = Path(queries_path)
    hmms_files = HMMFiles(*queries_path.iterdir())
    return hmms_files


def parse_genome(genome_path):
    genome_path = str(genome_path)
    genome = re.search(GENOME_REGEX, genome_path).group(1)
    return genome


def parse_hit(hit, genome_id):

    out = [genome_id] + [
        (pid := hit.name.decode("utf-8")),
        (query := hit.hits.query_accession.decode("utf-8")),
        (score := hit.score),
        (evalue := hit.evalue),
        (start := hit.best_domain.env_from),
        (end := hit.best_domain.env_to),
        (pid_txt := hit.description.decode("utf-8")),
        (query_txt := hit.hits.query_name.decode("utf-8")),
    ]

    out = [str(i) for i in out]

    return Results(*out)


def run_genomes(genome_paths, hmms_files):

    def run_genome(genome_path):
        genome_id = parse_genome(genome_path)

        with SequenceFile(genome_path, digital=True) as genome_file:
            genome = genome_file.read_block()
            # Esta funcion como tal abre hilos en el caso manejado, o eso me parece, significa que por si misma podria estar
            # abriendo nuevos subprocesos https://pyhmmer.readthedocs.io/en/stable/api/hmmer/profile.html#pyhmmer.hmmer.hmmsearch
            results = hmmsearch(hmms_files, genome, bit_cutoffs="trusted")

        hittup = []
        # results parece un generador, dado que la documentaci[on nos trae un "Yields"
        for top_hits in results:
            for hit in top_hits:
                # y dice que propiedad esta es un iterador de Hit ?
                # este comparador un es Truthy entonces? cuidado, tal vez verificar su longitud
                #if hit.included:
                #    parsed = parse_hit(hit, genome_id)
                #    hittup.append(parsed)
                if len(hit.included) < 1:
                    continue
                parsed = parse_hit(hit, genome_id)
                hittup.append(parsed)

        out = {}
        out[genome_id] = hittup

        return out

    merged = {}
    # Entonces map es secuencial dentro de esta llamada
    for resultD in map(run_genome, genome_paths):
        # no me conocia esta, como funciona en un diccionario? agrega las llaves nuevas? Como un setdefault?
        merged |= resultD

    return merged


if __name__ == "__main__":

    with open(GENOMES_FILE, "r") as genomes_file:
        genomes_paths = [Path(line.rstrip()) for line in genomes_file]
        # No se confirma que los paths existen?
        genomes_paths = [path for path in genomes_paths if path.exists() and path.is_file()]
        n_batches = len(genomes_paths)
        # no esta esto redundante? dado que haces un array shape (n,) y lo divides en n partes... mejor crea una lista de items individuales de arrays
        # batches = np.array_split(np.array(genomes_paths), n_batches)
        batches = [np.array(item) for item in genomes_paths]

    hmms_files = get_hmms(QUERIES_DIR)
    worker = partial(run_genomes, hmms_files=hmms_files)
    
    # added
    files_kwds = dict(hmms_files=hmms_files)
    with Pool() as pool:
        # sin determinar un chunksize no ganas mucho de imap contra map
        # results = pool.imap_unordered(worker, batches)
        # alternativa tambien puedes utilizar algo como apply para pasar args y kwargs sin crear el objeto partial
        promises = [pool.apply(run_genomes, args=(genome_path), kwds=files_kwds) for genome_path in batches]
        
        # resuelvelos
        for promise in promises: promise.ready()
        # Colectalos
        results = [promise.get() for promise in promises]
        
        # por que hacer esto cuando Pool se abre como un recurso, el chiste del scope es que se cierre y esas weas solo al terminar el scope
        # pool.close()
        # pool.join()

    merged = {}
    for result in results:
        merged |= result

    with open(OUT_FILE, "w") as tsv:
        w = lambda itsv: "\t".join(itsv) + "\n"

        tsv.write(w(FIELDS))

        for genome_id in merged:
            top_hits = merged[genome_id]

            for hittup in top_hits:
                tsv.write(w(hittup))
