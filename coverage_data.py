"""timdr_dna/coverage_data.py

Ladowanie/generowanie danych GLEBOKOSCI POKRYCIA (depth-of-coverage) -
standardowy format wejsciowy dla wykrywania CNV (copy-number variants:
delecji/duplikacji) w bioinformatyce, np. wyjscie `samtools depth`
(chrom, pozycja, glebokosc).

To NIE jest analiza surowej sekwencji ACGT - to analiza SZEREGU
CZASOWEGO glebokosci pokrycia wzdluz pozycji genomowej, dokladnie tak
jak realne narzedzia CNV (CNVnator, ReadDepth, CNVkit) - patrz README.
"""
from __future__ import annotations

import csv
from typing import Tuple

import numpy as np


def generate_synthetic_coverage(
    length: int = 6000,
    seed: int = 42,
    mean_depth: float = 30.0,
    deletion_region: Tuple[int, int] | None = (1500, 1700),
    duplication_region: Tuple[int, int] | None = (4200, 4450),
) -> Tuple[np.ndarray, np.ndarray]:
    """Generuje SYNTETYCZNA serie glebokosci pokrycia - NIE prawdziwe dane
    pacjenta. Symuluje: szum Poissona (typowy dla NGS), lagodny GC-bias
    (sinusoida) i opcjonalne wstrzykniete regiony delecji (glebokosc ~0)
    i duplikacji (glebokosc ~2-3x wyzsza) - do testowania i demonstracji.

    `deletion_region`/`duplication_region` = (start, end) w jednostkach
    pozycji (indeksy tablicy); None wylacza dany region.
    """
    rng = np.random.default_rng(seed)
    positions = np.arange(length)

    gc_bias = 1.0 + 0.15 * np.sin(2 * np.pi * positions / 800.0)
    base_depth = mean_depth * gc_bias
    depth = rng.poisson(lam=np.clip(base_depth, 1, None)).astype(float)

    if deletion_region is not None:
        s, e = deletion_region
        depth[s:e] = rng.poisson(lam=1.5, size=e - s).astype(float)

    if duplication_region is not None:
        s, e = duplication_region
        depth[s:e] = rng.poisson(lam=np.clip(base_depth[s:e] * 2.5, 1, None)).astype(float)

    return positions, depth


class CoverageFileError(Exception):
    """Czytelny blad parsowania pliku glebokosci - nigdy cichy pusty wynik
    (ta sama konwencja co DataLoaderError w analizator-gieldowy-v3)."""


def load_depth_tsv(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Wczytuje plik TSV w formacie `samtools depth` (chrom, pos, depth),
    bez naglowka. WYMAGA dokladnie jednego chromosomu w pliku - jesli
    plik zawiera wiecej niz jeden, rzuca czytelny blad (analiza
    zaklada CIAGLA os pozycji wzdluz jednego chromosomu; mieszanie
    chromosomow w jednej serii dalyby fikcyjne "skoki" na granicach).
    """
    positions, depths, chroms = [], [], set()
    try:
        with open(path, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or len(row) < 3:
                    continue
                chrom, pos, depth = row[0], int(row[1]), float(row[2])
                chroms.add(chrom)
                positions.append(pos)
                depths.append(depth)
    except (OSError, ValueError) as e:
        raise CoverageFileError(f"Nie udalo sie wczytac '{path}': {e}") from e

    if not positions:
        raise CoverageFileError(f"Plik '{path}' jest pusty lub w zlym formacie (oczekiwano TSV: chrom\\tpos\\tdepth).")
    if len(chroms) > 1:
        raise CoverageFileError(
            f"Plik '{path}' zawiera {len(chroms)} chromosomow ({sorted(chroms)}) - "
            "ta funkcja wymaga DOKLADNIE jednego (analizuj kazdy chromosom osobno)."
        )

    return np.asarray(positions, dtype=float), np.asarray(depths, dtype=float)
