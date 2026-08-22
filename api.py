"""api.py — TIMDR-DNA: lokalne REST API + dashboard (Flask).

Ten sam wzorzec co analizator-gieldowy-v3/api.py: pojedynczy proces
Flask serwujacy dashboard (Canvas 2D, bez CDN) + endpointy JSON.

WAZNE: patrz README.md - to NIE jest narzedzie diagnostyczne/kliniczne.
Domyslne zrodlo danych to SYNTETYCZNE dane demo (generate_synthetic_coverage) -
zadne prawdziwe dane pacjenta nie sa tu przetwarzane, chyba ze
uzytkownik jawnie poda sciezke do wlasnego pliku (source=file&path=...).
"""
from __future__ import annotations

import re

import numpy as np
from flask import Flask, jsonify, request, send_from_directory

from timdr_dna.coverage_data import (
    generate_synthetic_coverage,
    load_depth_tsv,
    CoverageFileError,
)
from timdr_dna.cnv_analyzer import analyze_coverage

app = Flask(__name__, static_folder="static", static_url_path="")

DISCLAIMER = (
    "Narzedzie badawczo-edukacyjne. NIE jest narzedziem diagnostycznym ani "
    "klinicznym i nie zastepuje zwalidowanych narzedzi CNV (CNVnator/CNVkit/GATK) "
    "ani konsultacji z genetykiem/lekarzem. Domyslnie pokazuje SYNTETYCZNE dane "
    "demo, nie dane pacjenta."
)

_REGION_RE = re.compile(r"^\s*(\d+)\s*-\s*(\d+)\s*$")


def _parse_region(raw: str | None):
    """'1500-1700' -> (1500, 1700); pusty/brak -> None. Rzuca ValueError
    z czytelnym komunikatem przy zlym formacie (nigdy cichy fallback)."""
    if not raw:
        return None
    m = _REGION_RE.match(raw)
    if not m:
        raise ValueError(f"Zly format regionu '{raw}' - oczekiwano 'start-koniec', np. '1500-1700'.")
    start, end = int(m.group(1)), int(m.group(2))
    if end <= start:
        raise ValueError(f"Region '{raw}': koniec musi byc wiekszy niz start.")
    return start, end


def _clean(obj):
    if isinstance(obj, float):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return _clean(float(obj))
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    return obj


@app.route("/")
def dashboard():
    return send_from_directory(app.static_folder, "dashboard.html")


@app.route("/api/health")
def health():
    return jsonify({"status": "ok", "disclaimer": DISCLAIMER})


@app.route("/api/analyze")
def analyze():
    source = request.args.get("source", "demo")
    window = int(request.args.get("window", 50))
    rezonans_min = int(request.args.get("rezonans_min", 2))

    try:
        if source == "file":
            path = request.args.get("path")
            if not path:
                return jsonify({"error": "parametr 'path' wymagany dla source=file"}), 400
            positions, depth = load_depth_tsv(path)
        else:
            length = int(request.args.get("length", 6000))
            seed = int(request.args.get("seed", 42))
            deletion_region = _parse_region(request.args.get("deletion", "1500-1700"))
            duplication_region = _parse_region(request.args.get("duplication", "4200-4450"))
            positions, depth = generate_synthetic_coverage(
                length=length, seed=seed,
                deletion_region=deletion_region,
                duplication_region=duplication_region,
            )

        result = analyze_coverage(positions, depth, window=window, rezonans_min=rezonans_min)
    except (CoverageFileError, ValueError) as e:
        return jsonify({"error": str(e)}), 400

    result.pop("raw_signal_result", None)  # zawiera np.ndarray zagniezdzone w slownikach - niepotrzebne dashboardowi
    result["source"] = source
    result["disclaimer"] = DISCLAIMER
    return jsonify(_clean(result))


def _open_browser_when_ready(url: str, delay: float = 1.2) -> None:
    """Otwiera przegladarke z krotkim opoznieniem, w osobnym watku.

    NAPRAWIONE wzgledem wzorca z analizator-gieldowy-v3/run.bat (`start ""
    http://...` PRZED odpaleniem serwera) - to jest dokladnie ten sam wyscig,
    ktory juz raz byl znaleziony i naprawiony w Synoptyk-v2.0/run.bat
    (przegladarka otwierala sie, zanim serwer zdazyl zaczac nasluchiwac,
    dajac "nie mozna polaczyc"). Tutaj otwieramy przegladarke z poziomu
    Pythona, PO starcie `app.run()` (w osobnym watku z opoznieniem), nie z
    poziomu .bat przed nim - ten sam duch naprawy, co `inbrowser=True` w
    gui_app.py Synoptyka.
    """
    import threading
    import webbrowser

    def _do_open():
        webbrowser.open(url)

    threading.Timer(delay, _do_open).start()


if __name__ == "__main__":
    # Port 8070 (nie 8060, ktorego uzywa analizator-gieldowy-v3, ani 8000/5060 -
    # patrz komentarze o kolizjach portow w innych repo tego zestawu) - zeby
    # dashboard TIMDR-DNA mogl dzialac obok innych bez konfliktu.
    _open_browser_when_ready("http://127.0.0.1:8070")
    app.run(host="127.0.0.1", port=8070, debug=False, threaded=True)
