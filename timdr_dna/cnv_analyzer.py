"""timdr_dna/cnv_analyzer.py

Wykrywanie KANDYDATOW na warianty liczby kopii (CNV: delecje/duplikacje)
z serii glebokosci pokrycia, tym samym duchem co reszta projektow TIMDR
w tym srodowisku (pogoda, gielda, akcelerator) - ZERO wiedzy
biologicznej zaszytej w silniku, tylko adaptacyjne progi liczone na
zywo z okna danych.

NAPRAWA/ODDZIELENIE (ta sesja): ten modul wczesniej WYMAGAL folderu-
siostry `universal-state-analyzer` (ten sam poziom katalogow) i rzucal
ImportError przy jego braku - patrz `timdr_dna/_engine.py` za pelne
uzasadnienie i koszt tej zmiany. TIMDR-DNA jest teraz samodzielny:
`pip install -r requirements.txt && pytest` dziala od razu po
sklonowaniu SAMEGO tego repo, bez niczego obok.

NOWA FUNKCJA (ta sesja): opcjonalna normalizacja wzgledem PROBKI
REFERENCYJNEJ (`reference_depth=`, np. dopasowana probka "normalna" w
terminologii onkologicznej - matched normal / panel of normals, ten
sam pomysl co w CNVkit) zamiast wzgledem wlasnej mediany. Patrz
docstring `analyze_coverage()` nizej za pelne uzasadnienie - to
naprawia realny, udokumentowany (ale wczesniej NIE zaadresowany w tym
repo) slepy punkt: wariant obejmujacy CALE albo prawie cale analizowane
okno jest niewidoczny dla normalizacji wzgledem wlasnej mediany, bo
mediana przesuwa sie razem z nim (ten sam "self-baseline blind spot"
opisany w skillu timdr-signal-framework i w
universal-state-analyzer/timdr_core/baseline.py).
"""
from __future__ import annotations

from typing import Tuple

import numpy as np

from timdr_dna._engine import anomalies, rezonans


def bin_coverage(positions: np.ndarray, depth: np.ndarray, window: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """Agreguje glebokosc w oknach po `window` pozycji (srednia) - standardowa
    praktyka w wykrywaniu CNV (redukuje szum pojedynczej zasady/odczytu,
    patrz README). Zwraca (window_center_position, window_mean_depth).
    """
    positions = np.asarray(positions, float)
    depth = np.asarray(depth, float)
    n = len(depth)
    if n == 0:
        return np.array([]), np.array([])

    n_windows = n // window
    if n_windows == 0:
        return np.array([positions.mean()]), np.array([depth.mean()])

    trimmed = n_windows * window
    depth_windows = depth[:trimmed].reshape(n_windows, window).mean(axis=1)
    pos_windows = positions[:trimmed].reshape(n_windows, window).mean(axis=1)
    return pos_windows, depth_windows


def analyze_coverage(
    positions: np.ndarray,
    depth: np.ndarray,
    window: int = 50,
    rezonans_min: int = 2,
    reference_depth: np.ndarray | None = None,
    anomaly_factor: float = 3.0,
) -> dict:
    """Pelny pipeline: bin -> log2-ratio -> anomalie(depth, log2_ratio) -> rezonans -> kandydaci CNV.

    Kanaly analizowane niezaleznie, potem laczone przez rezonans:
    - "depth"      : srednia glebokosc per okno (surowa skala)
    - "log2_ratio" : log2(okno / odniesienie) - standardowa normalizacja
                     w narzedziach CNV (CNVnator/CNVkit)

    Region flagowany jako kandydat, gdy OBA kanaly zgadzaja sie
    (rezonans >= rezonans_min) - pojedynczy kanal moze byc szumem,
    zgodnosc dwoch niezaleznych spojrzen na te sama dana jest silniejszym
    sygnalem (ten sam duch, co rezonans w Synoptyku/gieldzie).

    `reference_depth` (opcjonalne, NOWE): glebokosc pokrycia probki
    REFERENCYJNEJ na TYCH SAMYCH pozycjach co `depth` (ta sama dlugosc,
    ta sama os pozycji - np. dopasowana probka "normalna" dla tego
    samego regionu/panelu). Gdy podane, "odniesienie" powyzej to
    glebokosc referencyjna per okno zamiast wlasnej mediany calej serii.

    PO CO: bez referencji, log2_ratio liczy sie wzgledem WLASNEJ mediany
    analizowanej serii - to ma znany slepy punkt (patrz
    universal-state-analyzer/timdr_core/baseline.py i skill
    timdr-signal-framework, "self-baseline blind spot"): wariant
    obejmujacy CALE albo prawie cale okno przesuwa mediane razem ze soba
    i wychodzi jako "normalny", bo nie ma w danych nic, z czym go
    porownac. Z probka referencyjna, kazde okno jest porownywane do
    NIEZALEZNEGO pomiaru tej samej pozycji (nie do reszty tej samej
    serii) - dokladnie tak dziala matched-normal / panel-of-normals w
    realnych narzedziach CNV (CNVkit i inne).

    WAZNY SZCZEGOL: samo przeliczenie log2_ratio wzgledem referencji NIE
    WYSTARCZA - krok wykrywania anomalii (progowanie obu kanalow) MUSI
    TEZ liczyc prog wzgledem referencji, nie wzgledem wlasnej mediany
    log2_ratio/depth, inaczej ten sam slepy punkt wraca na poziomie
    progu (patrz komentarz w kodzie nizej i test dowodzacy tego w
    tests/test_timdr_dna.py). To zostalo tu zaimplementowane poprawnie -
    ta uwaga jest dla kogos, kto bedzie to modyfikowac w przyszlosci.

    OGRANICZENIE (uczciwie): to NIE rozwiazuje problemu, jesli sama
    probka referencyjna tez ma wariant w tym samym miejscu (np. oboje
    "case" i "reference" pochodza od tej samej osoby/linii komorkowej z
    tym samym wariantem dziedzicznym) - to fundamentalne ograniczenie
    kazdej metody porownawczej, nie cos, co da sie naprawic w kodzie.
    """
    if len(depth) < window * 4:
        raise ValueError(
            f"Za malo danych ({len(depth)} pozycji) na sensowna analize z window={window} "
            f"- potrzeba co najmniej {window * 4} pozycji (>=4 okna)."
        )

    pos_w, depth_w = bin_coverage(positions, depth, window=window)
    median_depth = np.median(depth_w)
    if median_depth <= 0:
        raise ValueError("Mediana glebokosci w oknach wynosi 0 - dane wygladaja na puste/zle sekwencjonowanie.")

    reference_used = reference_depth is not None
    median_reference_depth = None

    if reference_used:
        reference_depth = np.asarray(reference_depth, float)
        if len(reference_depth) != len(depth):
            raise ValueError(
                f"reference_depth ma {len(reference_depth)} pozycji, depth ma {len(depth)} - "
                "obie serie musza pokrywac ta sama, wyrownana os pozycji (probka referencyjna "
                "dla tego samego regionu/panelu co probka analizowana)."
            )
        _, reference_w = bin_coverage(positions, reference_depth, window=window)
        median_reference_depth = float(np.median(reference_w))
        if median_reference_depth <= 0:
            raise ValueError(
                "Mediana glebokosci w oknach probki referencyjnej wynosi 0 - referencja "
                "wyglada na pusta/zle sekwencjonowanie, nie nadaje sie jako odniesienie."
            )
        # podloga na referencji (nie na wyniku) - okno z prawie-zerowym
        # pokryciem referencyjnym dalioby niestabilny/eksplodujacy
        # stosunek, nie realny sygnal biologiczny
        reference_floor = max(median_reference_depth * 0.05, 1e-9)
        reference_safe = np.clip(reference_w, reference_floor, None)
        log2_ratio = np.log2(np.clip(depth_w, 1e-9, None) / reference_safe)
    else:
        log2_ratio = np.log2(np.clip(depth_w, 1e-9, None) / median_depth)

    if reference_used:
        # WAZNE (odkryte przy pisaniu testu na duzej delecji, patrz
        # tests/test_timdr_dna.py): samo przeliczenie log2_ratio wzgledem
        # referencji NIE WYSTARCZA, jesli krok wykrywania anomalii dalej
        # liczy mediane/MAD z SAMEGO log2_ratio (albo z samego depth_w) -
        # to odtwarza DOKLADNIE TEN SAM slepy punkt na poziomie progu:
        # region obejmujacy wiekszosc okna dominowalby wtedy WLASNA
        # mediane log2_ratio (blisko 0 = "normalny"), a maly normalny
        # fragment na brzegach wygladalby jak anomalia. Zeby probka
        # referencyjna faktycznie cos dala, OBA kanaly musza liczyc prog
        # z odniesienia do referencji, nie do samego siebie:
        # - depth: baseline = (mediana, MAD) glebokosci referencyjnej per okno
        # - log2_ratio: baseline = (0.0, oszacowany rozrzut) - "0" to
        #   dokladna wartosc oczekiwana przy braku zmiany (log2(x/x)=0),
        #   rozrzut oszacowany przez propagacje wzglednego szumu
        #   referencji (mad/mediana w skali liniowej) do skali log2
        #   (dzielenie przez ln(2)) - nie liczony z samego log2_ratio
        ref_med = float(np.median(reference_w))
        ref_mad = float(np.median(np.abs(reference_w - ref_med)) * 1.4826)
        if ref_mad == 0 or not np.isfinite(ref_mad):
            ref_std = float(np.std(reference_w))
            ref_mad = ref_std if ref_std > 0 else max(abs(ref_med) * 0.05, 1e-9)
        an_idx_depth, _, _ = anomalies(depth_w, factor=anomaly_factor, baseline=(ref_med, ref_mad))

        log2_mad = max(ref_mad / (ref_med * np.log(2)), 1e-6) if ref_med > 0 else 1e-6
        an_idx_log2, _, _ = anomalies(log2_ratio, factor=anomaly_factor, baseline=(0.0, log2_mad))
    else:
        an_idx_depth, _, _ = anomalies(depth_w, factor=anomaly_factor)
        an_idx_log2, _, _ = anomalies(log2_ratio, factor=anomaly_factor)

    rez_idx, rez_counts = rezonans([an_idx_depth, an_idx_log2], n=len(depth_w), min_count=rezonans_min)

    candidates = []
    for idx in rez_idx:
        kind = "mozliwa duplikacja" if log2_ratio[idx] > 0 else "mozliwa delecja"
        candidates.append({
            "window_index": int(idx),
            "position": float(pos_w[idx]),
            "depth": float(depth_w[idx]),
            "log2_ratio": float(log2_ratio[idx]),
            "magnitude": float(abs(log2_ratio[idx])),
            "kind": kind,
            "rezonans_count": int(rez_counts[idx]),
        })

    return {
        "window_positions": pos_w,
        "window_depth": depth_w,
        "log2_ratio": log2_ratio,
        "median_depth": float(median_depth),
        "reference_used": reference_used,
        "median_reference_depth": median_reference_depth,
        "candidates": candidates,
        "raw_signal_result": {
            "anomaly_idx": {"depth": an_idx_depth, "log2_ratio": an_idx_log2},
            "rezonans_idx": rez_idx,
            "rezonans_counts": rez_counts,
        },
    }
