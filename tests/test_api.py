"""tests/test_api.py - testy dla api.py (Flask), na danych demo (bez sieci,
bez plikow zewnetrznych)."""
import os
import tempfile

import pytest

from api import app


@pytest.fixture()
def client():
    app.config["TESTING"] = True
    with app.test_client() as c:
        yield c


def test_health(client):
    res = client.get("/api/health")
    assert res.status_code == 200
    assert res.get_json()["status"] == "ok"


def test_dashboard_root_serwuje_html(client):
    res = client.get("/")
    assert res.status_code == 200
    assert "text/html" in res.content_type


def test_analyze_demo_wykrywa_wstrzykniete_warianty(client):
    res = client.get(
        "/api/analyze?source=demo&length=6000&deletion=1500-1700"
        "&duplication=4200-4450&window=50&rezonans_min=2"
    )
    assert res.status_code == 200
    data = res.get_json()
    assert len(data["candidates"]) > 0
    assert any(c["kind"] == "mozliwa delecja" for c in data["candidates"])
    assert any(c["kind"] == "mozliwa duplikacja" for c in data["candidates"])
    assert "disclaimer" in data


def test_analyze_zly_format_regionu_zwraca_400(client):
    res = client.get("/api/analyze?source=demo&deletion=zle_dane")
    assert res.status_code == 400
    assert "error" in res.get_json()


def test_analyze_source_file_bez_path_zwraca_400(client):
    res = client.get("/api/analyze?source=file")
    assert res.status_code == 400
    assert "path" in res.get_json()["error"]


def test_analyze_source_file_niedostepny_plik_zwraca_400(client):
    res = client.get("/api/analyze?source=file&path=/nieistniejacy/plik.tsv")
    assert res.status_code == 400
    assert "error" in res.get_json()


def test_analyze_demo_z_referencja_zwraca_reference_used(client):
    res = client.get(
        "/api/analyze?source=demo&length=6000&deletion=300-5700&duplication="
        "&window=50&rezonans_min=2&reference=1"
    )
    assert res.status_code == 200
    data = res.get_json()
    assert data["reference_used"] is True
    assert data["median_reference_depth"] is not None and data["median_reference_depth"] > 0


def test_analyze_demo_bez_referencji_zwraca_reference_used_false(client):
    res = client.get("/api/analyze?source=demo&length=6000&window=50&rezonans_min=2")
    assert res.status_code == 200
    data = res.get_json()
    assert data["reference_used"] is False
    assert data["median_reference_depth"] is None


def test_analyze_source_file_niedostepny_plik_referencyjny_zwraca_400(client):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        for pos in range(300):
            f.write(f"chr1\t{pos}\t{20 + pos % 5}\n")
        path = f.name
    try:
        res = client.get(
            "/api/analyze?source=file&path=" + path + "&reference_path=/nieistniejacy/ref.tsv"
        )
        assert res.status_code == 400
        assert "error" in res.get_json()
    finally:
        os.unlink(path)
