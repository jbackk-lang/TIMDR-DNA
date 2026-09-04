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
    # deletion_region/duplication_region jawnie None - ten test sprawdza
    # tylko ksztalt/nieujemnosc, nie warianty, wiec nie powinien zalezec
    # od domyslnych regionow (1500,1700)/(4200,4450), ktore i tak nie
    # miesza sie w length=2000. PRZED DODANIEM WALIDACJI ZAKRESU REGIONU
    # (ta sesja) ten test przechodzil przez przypadek: deletion mieści
    # się w 2000 pozycjach, duplication NIE mieści się, ale poza-zakresowe
    # przypisanie dla duplication milczaco nie robilo nic (rng.poisson()
    # z lam=tablica dopasowuje ksztalt do pustej tablicy, w przeciwienstwie
    # do lam=liczba+size, ktore twardo rzuca blad przy niezgodnosci) -
    # asymetryczne, przypadkowe zachowanie, nie zamierzony kontrakt.
    positions, depth = generate_synthetic_coverage(
        length=2000, seed=1, deletion_region=None, duplication_region=None,
    )
    assert len(positions) == len(depth) == 2000
    assert np.all(depth >= 0)


def test_generate_synthetic_coverage_delecja_ma_nizsza_glebokosc():
    positions, depth = generate_synthetic_coverage(
        length=2000, deletion_region=(500, 700), duplication_region=None,
    )
    poza_delecja = np.median(np.concatenate([depth[:500], depth[700:]]))
    w_delecji = np.median(depth[500:700])
    assert w_delecji < poza_delecja * 0.3  # wyraznie nizsza


def test_generate_synthetic_coverage_region_poza_zakresem_rzuca_czytelny_blad():
    # znalezione przy pisaniu testow referencji w tej sesji: domyslny
    # deletion_region=(1500,1700) z length=1000 kiedys dawal nieczytelny
    # blad numpy o niezgodnosci ksztaltow - teraz jasny komunikat
    with pytest.raises(ValueError, match="deletion_region"):
        generate_synthetic_coverage(length=1000, deletion_region=(1500, 1700), duplication_region=None)


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


# ---------------------------------------------------------------
# analyze_coverage(reference_depth=...) - probka referencyjna
# (matched-normal), naprawia self-baseline blind spot dla wariantow
# obejmujacych wiekszosc analizowanego okna
# ---------------------------------------------------------------
def test_analyze_coverage_z_referencja_wykrywa_wiecej_niz_self_baseline_przy_ogromnej_delecji():
    """GLOWNY test tej funkcji: delecja obejmujaca ~90% calej serii jest
    znanym slepym punktem normalizacji wzgledem wlasnej mediany (mediana
    przesuwa sie razem z delecja, delecja wyglada jak "normalny" poziom -
    patrz docstring analyze_coverage() i timdr_core/baseline.py w
    universal-state-analyzer). Probka referencyjna BEZ tej delecji
    (ten sam seed, wiec identyczna poza wstrzykniętym regionem - patrz
    generate_synthetic_coverage()) powinna wykryc wyraznie wiecej
    delecji WEWNATRZ tego regionu niz self-baseline."""
    length = 6000
    big_deletion = (300, 5700)  # ~90% z 6000 pozycji

    positions, case_depth = generate_synthetic_coverage(
        length=length, seed=42, deletion_region=big_deletion, duplication_region=None,
    )
    _, reference_depth = generate_synthetic_coverage(
        length=length, seed=42, deletion_region=None, duplication_region=None,
    )

    self_result = analyze_coverage(positions, case_depth, window=50, rezonans_min=2)
    ref_result = analyze_coverage(
        positions, case_depth, window=50, rezonans_min=2, reference_depth=reference_depth,
    )

    def delecje_w_regionie(result):
        return [
            c for c in result["candidates"]
            if c["kind"] == "mozliwa delecja" and big_deletion[0] <= c["position"] <= big_deletion[1]
        ]

    self_delecje = delecje_w_regionie(self_result)
    ref_delecje = delecje_w_regionie(ref_result)

    assert self_result["reference_used"] is False
    assert ref_result["reference_used"] is True
    assert ref_result["median_reference_depth"] is not None and ref_result["median_reference_depth"] > 0

    assert len(ref_delecje) > len(self_delecje), (
        "probka referencyjna powinna wykryc WIECEJ delecji wewnatrz ogromnego "
        f"regionu niz self-baseline (self={len(self_delecje)}, ref={len(ref_delecje)}) - "
        "jesli to sie nie powiodlo, sprawdz czy baseline= faktycznie jest przekazywane "
        "do OBU kanalow (depth i log2_ratio) w analyze_coverage(), nie tylko do jednego"
    )
    assert len(ref_delecje) >= 20, (
        f"probka referencyjna powinna wykryc wiekszosc z ~108 okien w regionie delecji, "
        f"wykryto tylko {len(ref_delecje)}"
    )


def test_analyze_coverage_referencja_zlej_dlugosci_rzuca_czytelny_blad():
    # deletion_region/duplication_region jawnie None - domyslne (1500,1700)/
    # (4200,4450) z generate_synthetic_coverage() nie mieszcza sie w
    # length=1000 (blad tego typu zlapany przy pierwszym uruchomieniu tych
    # testow - patrz historia gita/sesji)
    positions, depth = generate_synthetic_coverage(
        length=1000, deletion_region=None, duplication_region=None,
    )
    reference_zla_dlugosc = np.zeros(500)
    with pytest.raises(ValueError, match="reference_depth"):
        analyze_coverage(positions, depth, window=50, reference_depth=reference_zla_dlugosc)


def test_analyze_coverage_referencja_zerowa_mediana_rzuca_czytelny_blad():
    positions, depth = generate_synthetic_coverage(
        length=1000, deletion_region=None, duplication_region=None,
    )
    reference_pusta = np.zeros(1000)
    with pytest.raises(ValueError, match="(?i)referencyjnej"):
        analyze_coverage(positions, depth, window=50, reference_depth=reference_pusta)


def test_analyze_coverage_bez_referencji_ma_reference_used_false():
    positions, depth = generate_synthetic_coverage(
        length=1000, deletion_region=None, duplication_region=None,
    )
    result = analyze_coverage(positions, depth, window=50)
    assert result["reference_used"] is False
    assert result["median_reference_depth"] is None


def test_analyze_coverage_kandydaci_maja_pole_magnitude():
    positions, depth = generate_synthetic_coverage(length=6000)
    result = analyze_coverage(positions, depth, window=50, rezonans_min=2)
    for c in result["candidates"]:
        assert c["magnitude"] == pytest.approx(abs(c["log2_ratio"]))
