"""timdr_dna/_engine.py

Silnik anomalii/rezonansu UZYWANY przez cnv_analyzer.py - CELOWO
zdublowany, zminimalizowany fork dwoch funkcji z
`universal-state-analyzer/timdr_core/core.py` (TIMDRCore.anomalies()
i TIMDRCore.rezonans()), zeby TIMDR-DNA dzialalo SAMODZIELNIE, bez
wymogu klonowania folderu-siostry `universal-state-analyzer` obok
niego.

DLACZEGO TAK (i jaki jest tego koszt) - patrz tez notatka o
"duplication-drift" w skillu timdr-signal-framework tego srodowiska:
zamierzone zdublowanie kodu miedzy repo (w odroznieniu od PRZYPADKOWEGO
zdublowania W JEDNYM repo, ktore jest bledem) jest tu swiadomym
kompromisem. TIMDR-DNA zyskuje: dziala samodzielnie po sklonowaniu
(`pip install -r requirements.txt && pytest`), bez ImportError przy
braku sasiedniego repo, latwiej dystrybuowac/publikowac osobno. Traci:
NIE dostanie automatycznie przyszlych poprawek/ulepszen silnika z
`universal-state-analyzer/timdr_core/core.py` (np. gdyby ktos naprawil
tam blad w MAD-scale albo w podlodze floor_frac) - trzeba by je
recznie przeniesc tutaj. To swiadomy wybor, nie przeoczenie.

CO ZOSTALO POMINIETE wzgledem oryginalnego TIMDRCore.analyze_multi()
(defekt/twist/flow/trm): analyze_coverage() w cnv_analyzer.py NIGDY
nie czytal wynikow defekt_idx/twist_idx z wyniku analyze_multi() -
tylko rezonans_idx/rezonans_counts (patrz stara wersja cnv_analyzer.py
w historii gita commitu tego repo). Zamiast kopiowac nieuzywany kod,
przeniesiono tylko to, co faktycznie zasila kandydatow CNV: anomalies()
per kanal + rezonans() miedzy kanalami. Jesli w przyszlosci CNV-analyzer
zacznie potrzebowac defekt()/twist() (np. do wykrywania ostrych granic
segmentow), doloz je tu w ten sam sposob (skopiuj z core.py, zacytuj
zrodlo w docstringu funkcji).

Zrodlo tego forka: universal-state-analyzer/timdr_core/core.py, metody
`anomalies()` (WLACZNIE z `baseline=` - patrz cnv_analyzer.py za nowy
sposob jego uzycia przez probke referencyjna/dopasowana-normalna, co
omija slepy punkt `baseline=None`) i `rezonans()`, skopiowane 1:1 co do
logiki (zweryfikowane recznie linia po linii przy tworzeniu tego
pliku).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np


def anomalies(
    s: np.ndarray,
    factor: float = 3.0,
    floor_frac: float = 0.05,
    mad_scale: float = 1.4826,
    baseline: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """MAD-owy z-score z podloga.

    `baseline=None` (domyslnie): median/MAD liczone z TEGO SAMEGO `s`,
    ktore jest oceniane ("self"). ZNANY slepy punkt (identyczny jak w
    TIMDRCore.anomalies() z baseline=None - patrz jego docstring w
    universal-state-analyzer): wariant obejmujacy CALE albo prawie cale
    ocenianie okno wychodzi jako "normalny", bo mediana/MAD samego okna
    przesuwaja sie razem z nim.

    `baseline=(median, mad)`: jawnie podane odniesienie, policzone
    GDZIE INDZIEJ niz `s` (w cnv_analyzer.py: z probki referencyjnej,
    nie z ocenianej probki) - to jest sposob na obejscie powyzszego
    slepego punktu, tak samo jak `baseline=` w oryginalnym
    TIMDRCore.anomalies().

    Zwraca (idx_powyzej_progu, z_scores, prog_bezwzgledny).
    """
    s = np.asarray(s, float)
    n = len(s)
    if n == 0:
        return np.array([], dtype=int), np.array([]), 0.0
    if baseline is not None:
        med, mad = baseline
        if mad == 0 or not np.isfinite(mad):
            mad = max(abs(med) * floor_frac, 1e-9)
    else:
        med = np.median(s)
        mad = np.median(np.abs(s - med)) * mad_scale
        if mad == 0 or not np.isfinite(mad):
            std = np.std(s)
            mad = std if std > 0 else max(abs(med) * floor_frac, 1e-9)
    z = (s - med) / mad
    idx = np.where(np.abs(z) > factor)[0]
    return idx, z, mad * factor


def rezonans(
    anomaly_index_lists: Sequence[np.ndarray],
    n: int,
    min_count: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """>=min_count list z anomaly_index_lists flaguje ten sam indeks -
    zgodnosc kilku niezaleznych kanalow/spojrzen na te sama dana jest
    silniejszym sygnalem niz pojedynczy kanal (ktory moze byc szumem).

    Zwraca (idx_rezonansowe, counts) gdzie counts[i] = ile kanalow
    flagowalo pozycje i.
    """
    counts = np.zeros(n, dtype=int)
    for idxs in anomaly_index_lists:
        counts[np.asarray(idxs, dtype=int)] += 1
    idx = np.where(counts >= min_count)[0]
    return idx, counts
