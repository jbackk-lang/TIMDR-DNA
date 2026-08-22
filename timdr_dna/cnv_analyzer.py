"""timdr_dna/cnv_analyzer.py

Wykrywanie KANDYDATOW na warianty liczby kopii (CNV: delecje/duplikacje)
z serii glebokosci pokrycia, tym samym generycznym silnikiem TIMDR co
reszta projektow w tym srodowisku (pogoda, gielda, akcelerator) -
ZERO wiedzy biologicznej zaszytej w silniku, tylko adaptacyjne progi
liczone na zywo z okna danych (patrz timdr_core w universal-state-analyzer).

Wymaga folderu-siostry `universal-state-analyzer` (ten sam poziom
katalogow) - patrz _ensure_timdr_core_on_path().
"""
from __future__ import annotations

import os
import sys
from typing import List, Tuple

import numpy as np


def _ensure_timdr_core_on_path() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    sibling = os.path.abspath(os.path.join(here, "..", "..", "universal-state-analyzer"))
    if not os.path.isdir(sibling):
        raise ImportError(
            "timdr_dna wymaga folderu-siostry 'universal-state-analyzer' "
            f"(szukano w: {sibling}) - stamtad pochodzi generyczny silnik TIMDRCore."
        )
    if sibling not in sys.path:
        sys.path.insert(0, sibling)


_ensure_timdr_core_on_path()
from timdr_core import TIMDRCore  # noqa: E402


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
) -> dict:
    """Pelny pipeline: bin -> log2-ratio -> TIMDRCore.analyze_multi -> kandydaci CNV.

    Kanaly przekazywane do analyze_multi:
    - "depth"      : srednia glebokosc per okno (surowa skala)
    - "log2_ratio" : log2(okno / mediana calej serii) - standardowa
                     normalizacja w narzedziach CNV (CNVnator/CNVkit),
                     zeby wartosc byla porownywalna niezaleznie od
                     bezwzglednej glebokosci sekwencjonowania probki.

    Region flagowany jako kandydat, gdy OBA kanaly zgadzaja sie
    (rezonans >= rezonans_min) - pojedynczy kanal moze byc szumem,
    zgodnosc dwoch niezaleznych spojrzen na te sama dana jest silniejszym
    sygnalem (ten sam duch, co rezonans w Synoptyku/gieldzie).
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

    log2_ratio = np.log2(np.clip(depth_w, 1e-9, None) / median_depth)

    core = TIMDRCore()
    result = core.analyze_multi(
        pos_w,
        {"depth": depth_w, "log2_ratio": log2_ratio},
        rezonans_min=rezonans_min,
    )

    candidates = []
    for idx in result["rezonans_idx"]:
        kind = "mozliwa duplikacja" if log2_ratio[idx] > 0 else "mozliwa delecja"
        candidates.append({
            "window_index": int(idx),
            "position": float(pos_w[idx]),
            "depth": float(depth_w[idx]),
            "log2_ratio": float(log2_ratio[idx]),
            "kind": kind,
            "rezonans_count": int(result["rezonans_counts"][idx]),
        })

    return {
        "window_positions": pos_w,
        "window_depth": depth_w,
        "log2_ratio": log2_ratio,
        "median_depth": float(median_depth),
        "candidates": candidates,
        "raw_signal_result": result,
    }
