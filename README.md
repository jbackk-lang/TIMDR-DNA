# TIMDR-DNA

> **TO NIE JEST NARZĘDZIE DIAGNOSTYCZNE ANI KLINICZNE.** Wykrywa
> statystyczne odchylenia w syntetycznych/testowych danych głębokości
> pokrycia sekwencjonowania - nie jest zwalidowane klinicznie, nie
> zastępuje zatwierdzonych narzędzi CNV (CNVnator, CNVkit, GATK) ani
> konsultacji z genetykiem klinicznym/lekarzem. Wyłącznie do celów
> badawczych/edukacyjnych nad frameworkiem sygnałowym TIMDR. Nie zawiera
> ani nie przetwarza żadnych prawdziwych danych pacjentów.

## Co to robi

**TIMDR-DNA** wykrywa kandydatów na warianty
liczby kopii (CNV - copy-number variants: delecje/duplikacje) z
prawdziwego formatu danych bioinformatycznych - głębokości pokrycia
(depth-of-coverage) wzdłuż pozycji genomowej, tym samym generycznym
silnikiem adaptacyjnych progów (TIMDRCore), który stoi za projektami
pogodowym (Synoptyk) i giełdowym (analizator-gieldowy-v3) w tym
środowisku.

## Grunt bioinformatyczny (uczciwie o pierwowzorach)

Wykrywanie CNV z głębokości pokrycia to REALNA, ugruntowana technika -
nie wymyślona tutaj. Istniejące, zwalidowane narzędzia oparte na tej
samej idei: **CNVnator**, **ReadDepth**, **CNVkit** - wszystkie liczą
statystyki głębokości pokrycia w oknach wzdłuż chromosomu i szukają
regionów odstających od normy. Ten projekt NIE konkuruje z nimi jako
narzędzie kliniczne - używa tej samej ogólnej idei (okno + statystyka
odstająca) jako demonstracji generycznego silnika TIMDR na realnym
kształcie danych bioinformatycznych, z uczciwym podejściem: adaptacyjny
próg liczony z samego analizowanego okna (MAD-z, rozstęp p90-p10),
zamiast jednego wytrenowanego modelu ML.

## Jak to działa

1. **`coverage_data.py`** - wczytuje głębokość pokrycia z pliku TSV w
   formacie `samtools depth` (chrom, pozycja, głębokość) ALBO generuje
   syntetyczną serię demo (szum Poissona + łagodny GC-bias + opcjonalne
   wstrzyknięte regiony delecji/duplikacji).
2. **`cnv_analyzer.bin_coverage()`** - agreguje głębokość w oknach
   (domyślnie 50 pozycji) - standardowa praktyka w CNV, redukuje szum
   pojedynczej zasady.
3. **`cnv_analyzer.analyze_coverage()`** - liczy `log2(okno / mediana
   całej serii)` (standardowa normalizacja z CNVnator/CNVkit), puszcza
   DWA kanały (`depth`, `log2_ratio`) przez `TIMDRCore.analyze_multi()`
   z `universal-state-analyzer` (folder-siostra, ten sam poziom
   katalogów) i flaguje **kandydata CNV** tam, gdzie OBA kanały
   zgadzają się (rezonans) - pojedynczy kanał może być szumem, zgodność
   dwóch niezależnych spojrzeń na te same dane jest silniejszym
   sygnałem.
4. Kandydat jest oznaczany jako "możliwa delecja" (log2_ratio < 0) albo
   "możliwa duplikacja" (log2_ratio > 0).

## Dashboard + API (Flask)

```bash
run.bat                      # Windows: instaluje zaleznosci, odpala serwer
# lub recznie:
pip install -r requirements.txt
python api.py                 # http://127.0.0.1:8070
```

Endpointy:
- `GET /` - dashboard (Canvas 2D, ciemny motyw, bez CDN - ten sam styl co `analizator-gieldowy-v3`)
- `GET /api/health`
- `GET /api/analyze?source=demo&length=6000&deletion=1500-1700&duplication=4200-4450&window=50&rezonans_min=2`
- `GET /api/analyze?source=file&path=chr1_depth.tsv&window=50` - realny plik `samtools depth`

Dashboard domyslnie laduje dane demo z wstrzyknieta delecja i duplikacja,
rysuje `log2_ratio` wzdluz pozycji i oznacza wykrytych kandydatow
kolorowymi punktami (czerwony = delecja, pomaranczowy = duplikacja) +
tabela z lista kandydatow. Port 8070 (nie 8060 jak `analizator-gieldowy-v3`,
zeby oba dzialaly rownolegle bez kolizji).

## Użycie (programowo / bez dashboardu)

```bash
pip install -r requirements.txt
pytest tests/ -q   # 19 testów (13 rdzen + 6 API)
```

```python
from timdr_dna.coverage_data import generate_synthetic_coverage, load_depth_tsv
from timdr_dna.cnv_analyzer import analyze_coverage

# demo (syntetyczne dane, NIE prawdziwy pacjent)
positions, depth = generate_synthetic_coverage(length=6000)
result = analyze_coverage(positions, depth, window=50)
for c in result["candidates"]:
    print(c["kind"], "przy pozycji", c["position"], "log2_ratio=", round(c["log2_ratio"], 2))

# na realnym pliku (samtools depth chr1.bam > chr1_depth.tsv)
positions, depth = load_depth_tsv("chr1_depth.tsv")
result = analyze_coverage(positions, depth, window=50)
```

## Testy

19 testów łącznie: 13 w `tests/test_timdr_dna.py` (silnik/dane) + 6 w `tests/test_api.py` (endpointy Flask, w tym błędne parametry -> 400).

`tests/test_timdr_dna.py`: generowanie syntetycznych danych
(kształt, że delecja/duplikacja faktycznie mają inną głębokość),
`bin_coverage()` (agregacja, pusta seria), `analyze_coverage()`
integracyjnie - **wykrywa wstrzyknięte, znane z góry regiony delecji i
duplikacji** (test najważniejszy: `test_analyze_coverage_wykrywa_wstrzykniete_delecje_i_duplikacje`),
za mało danych / zerowa mediana (czytelne błędy, nie cichy crash),
`load_depth_tsv()` (poprawny plik, wiele chromosomów = błąd, pusty plik
= błąd, brakujący plik = błąd - nigdy `except: pass`).

## Czego to NIE jest

- Nie jest to narzędzie diagnostyczne/kliniczne (patrz zastrzeżenie na
  górze). Prawdziwe wykrywanie CNV do celów klinicznych wymaga
  zwalidowanych narzędzi (CNVnator/CNVkit/GATK), kontroli jakości
  sekwencjonowania, korekty na mapowalność/GC-bias znacznie
  dokładniejszej niż tu, i interpretacji przez klinicystę/genetyka.
- Progi (rezonans_min, okno=50) są arbitralne, nie skalibrowane na
  żadnym realnym zbiorze danych klinicznych - jak wszędzie w TIMDR,
  wymagają dostrojenia do konkretnego zastosowania.
- Nie zawiera i nie przetwarza żadnych prawdziwych danych pacjentów -
  wszystkie testy i przykłady używają syntetycznych danych
  wygenerowanych przez `generate_synthetic_coverage()`.

## Licencja

MIT — patrz [LICENSE](LICENSE).
