"""tests/test_timdr_dna.py - testy TIMDR-DNA (wykrywanie kandydatow CNV
z glebokosci pokrycia), na syntetycznych danych z wstrzykniete znanymi
regionami delecji/duplikacji."""
from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

from timdr_dna.coverage_data import (
    generate_synthetic_coverage,
    load_depth_tsv,
    CoverageFileError,
)
from timdr_dna.cnv_analyzer import bin_coverage, analyze_coverage


# ---------------------------------------------------------------
# generate_synthetic_coverage
# ---------------------------------------------------------------
def test_generate_synthetic_coverage_ksztalt():
    positions, depth = generate_synthetic_coverage(length=2000, seed=1)
    assert len(positions) == len(depth) == 2000
    assert np.all(depth >= 0)


def test_generate_synthetic_coverage_delecja_ma_nizsza_glebokosc():
    positions, depth = generate_synthetic_coverage(
        length=2000, deletion_region=(500, 700), duplication_region=None,
    )
    poza_delecja = np.median(np.concatenate([depth[:500], depth[700:]]))
    w_delecji = np.median(depth[500:700])
    assert w_delecji < poza_delecja * 0.3  # wyraznie nizsza


def test_generate_synthetic_coverage_duplikacja_ma_wyzsza_glebokosc():
    positions, depth = generate_synthetic_coverage(
        length=2000, deletion_region=None, duplication_region=(1000, 1200),
    )
    poza = np.median(np.concatenate([depth[:1000], depth[1200:]]))
    w_dup = np.median(depth[1000:1200])
    assert w_dup > poza * 1.5


# ---------------------------------------------------------------
# bin_coverage
# ---------------------------------------------------------------
def test_bin_coverage_dlugosc_i_srednia():
    positions = np.arange(1000)
    depth = np.full(1000, 10.0)
    pos_w, depth_w = bin_coverage(positions, depth, window=100)
    assert len(pos_w) == len(depth_w) == 10
    assert np.allclose(depth_w, 10.0)


def test_bin_coverage_pusta_seria():
    pos_w, depth_w = bin_coverage(np.array([]), np.array([]), window=50)
    assert len(pos_w) == 0 and len(depth_w) == 0


# ---------------------------------------------------------------
# analyze_coverage - GLOWNY test integracyjny
# ---------------------------------------------------------------
def test_analyze_coverage_wykrywa_wstrzykniete_delecje_i_duplikacje():
    positions, depth = generate_synthetic_coverage(
        length=6000, seed=42,
        deletion_region=(1500, 1700),
        duplication_region=(4200, 4450),
    )
    result = analyze_coverage(positions, depth, window=50, rezonans_min=2)

    assert result["candidates"], "powinien byc co najmniej jeden kandydat CNV"

    delecje = [c for c in result["candidates"] if c["kind"] == "mozliwa delecja"
               and 1500 <= c["position"] <= 1700]
    duplikacje = [c for c in result["candidates"] if c["kind"] == "mozliwa duplikacja"
                  and 4200 <= c["position"] <= 4450]

    assert delecje, f"nie wykryto delecji w regionie 1500-1700, kandydaci: {result['candidates']}"
    assert duplikacje, f"nie wykryto duplikacji w regionie 4200-4450, kandydaci: {result['candidates']}"


def test_analyze_coverage_bez_wariantow_ma_malo_lub_zero_kandydatow():
    positions, depth = generate_synthetic_coverage(
        length=3000, seed=7, deletion_region=None, duplication_region=None,
    )
    result = analyze_coverage(positions, depth, window=50, rezonans_min=2)
    # gladka, "czysta" seria - nie powinna dawac lawiny kandydatow
    assert len(result["candidates"]) <= 3


def test_analyze_coverage_za_malo_danych_rzuca_czytelny_blad():
    positions, depth = generate_synthetic_coverage(
        length=50, deletion_region=None, duplication_region=None,
    )
    with pytest.raises(ValueError, match="Za malo danych"):
        analyze_coverage(positions, depth, window=50)


def test_analyze_coverage_zerowa_mediana_rzuca_czytelny_blad():
    positions = np.arange(500)
    depth = np.zeros(500)
    with pytest.raises(ValueError, match="(?i)mediana"):
        analyze_coverage(positions, depth, window=50)


# ---------------------------------------------------------------
# load_depth_tsv
# ---------------------------------------------------------------
def test_load_depth_tsv_poprawny_plik():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        for pos in range(10):
            f.write(f"chr1\t{pos}\t{20 + pos}\n")
        path = f.name
    try:
        positions, depth = load_depth_tsv(path)
        assert len(positions) == 10
        assert depth[0] == 20.0
    finally:
        os.unlink(path)


def test_load_depth_tsv_wiele_chromosomow_rzuca_blad():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("chr1\t1\t20\n")
        f.write("chr2\t1\t25\n")
        path = f.name
    try:
        with pytest.raises(CoverageFileError, match="chromosomow"):
            load_depth_tsv(path)
    finally:
        os.unlink(path)


def test_load_depth_tsv_pusty_plik_rzuca_blad():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        path = f.name
    try:
        with pytest.raises(CoverageFileError):
            load_depth_tsv(path)
    finally:
        os.unlink(path)


def test_load_depth_tsv_brakujacy_plik_rzuca_blad():
    with pytest.raises(CoverageFileError):
        load_depth_tsv("/nieistniejacy/plik.tsv")
