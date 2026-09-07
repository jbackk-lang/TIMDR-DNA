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
(depth-of-coverage) wzdłuż pozycji genomowej, tym samym duchem
adaptacyjnych progów co projekty pogodowy (Synoptyk) i giełdowy
(analizator-gieldowy-v3) w tym środowisku.

**Samodzielny.** Wcześniej ten silnik był importowany z folderu-siostry
`universal-state-analyzer` (wymagał sklonowania obu repo obok siebie).
Od tej sesji `timdr_dna/_engine.py` zawiera własny, zminimalizowany
fork dwóch potrzebnych funkcji (`anomalies()`, `rezonans()`) - patrz
docstring tego pliku za pełne uzasadnienie i świadomy koszt tej decyzji
(nie dostaje automatycznie przyszłych poprawek silnika z tamtego repo).
`pip install -r requirements.txt && pytest` działa od razu po
sklonowaniu samego TIMDR-DNA.

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
3. **`cnv_analyzer.analyze_coverage()`** - liczy `log2(okno / odniesienie)`
   (standardowa normalizacja z CNVnator/CNVkit), puszcza DWA kanały
   (`depth`, `log2_ratio`) przez `anomalies()`+`rezonans()` z
   `timdr_dna/_engine.py` (własny, samodzielny silnik - patrz wyżej) i
   flaguje **kandydata CNV** tam, gdzie OBA kanały zgadzają się
   (rezonans) - pojedynczy kanał może być szumem, zgodność dwóch
   niezależnych spojrzeń na te same dane jest silniejszym sygnałem.
4. Kandydat jest oznaczany jako "możliwa delecja" (log2_ratio < 0) albo
   "możliwa duplikacja" (log2_ratio > 0), z polem `magnitude` (=
   |log2_ratio|) do szybkiego sortowania/triage.
5. **Opcjonalnie: próbka referencyjna** (`reference_depth=`, NOWE) -
   "odniesienie" w kroku 3 to głębokość próbki referencyjnej (matched
   normal / panel of normals, jak w CNVkit) zamiast własnej mediany
   analizowanej serii. To naprawia realny, znany "self-baseline blind
   spot" (patrz `timdr-signal-framework`/`universal-state-analyzer/timdr_core/baseline.py`):
   wariant obejmujący całe albo prawie całe okno przesuwa własną
   medianę razem ze sobą i wychodzi jako "normalny" - z niezależną
   próbką referencyjną każde okno jest porównywane do zewnętrznego
   pomiaru tej samej pozycji, nie do reszty tej samej serii. Pełne
   uzasadnienie (i dlaczego samo przeliczenie ratio nie wystarcza, próg
   detekcji też musi liczyć się względem referencji) - w docstringu
   `analyze_coverage()` i w teście
   `test_analyze_coverage_z_referencja_wykrywa_wiecej_niz_self_baseline_przy_ogromnej_delecji`.

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
- `GET /api/analyze?source=demo&...&reference=1` - jak wyzej + syntetyczna probka referencyjna (bez wariantow) do normalizacji, patrz sekcja o probce referencyjnej
- `GET /api/analyze?source=file&path=chr1_depth.tsv&window=50` - realny plik `samtools depth`
- `GET /api/analyze?source=file&path=chr1_case.tsv&reference_path=chr1_normal.tsv&window=50` - jak wyzej + realny plik probki referencyjnej (musi miec te sama liczbe pozycji co `path`)

Dashboard domyslnie laduje dane demo z wstrzyknieta delecja i duplikacja,
rysuje `log2_ratio` wzdluz pozycji i oznacza wykrytych kandydatow
kolorowymi punktami (czerwony = delecja, pomaranczowy = duplikacja) +
tabela z lista kandydatow (z kolumna `magnitude` do szybkiego triage).
Pole "próba referencyjna" w formularzu (checkbox przy danych demo,
sciezka do pliku przy `source=file`) wlacza normalizacje wzgledem
referencji zamiast wlasnej mediany - karta "Odniesienie" i notka pod
formularzem pokazuja, ktory tryb byl uzyty. Port 8070 (nie 8060 jak
`analizator-gieldowy-v3`, zeby oba dzialaly rownolegle bez kolizji).

## Użycie (programowo / bez dashboardu)

```bash
pip install -r requirements.txt
pytest tests/ -q   # 28 testow (19 rdzen + 9 API)
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

28 testów łącznie: 19 w `tests/test_timdr_dna.py` (silnik/dane) + 9 w `tests/test_api.py` (endpointy Flask, w tym błędne parametry -> 400).

> UWAGA O WYKONANIU: **zweryfikowane** - użytkownik uruchomił `pytest` trzy
> razy podczas dodawania tych zmian (odseparowanie silnika + próbka referencyjna).
> 1. podejście: 25/28 - 3 nowe testy referencji użyły `length=1000` bez
> wyłączenia domyślnych `deletion_region=(1500,1700)`/`duplication_region=(4200,4450)`
> z `generate_synthetic_coverage()`, które nie mieszczą się w tak krótkiej serii -
> naprawione w testach + dodana jawna walidacja zakresu regionu w
> `generate_synthetic_coverage()` (`ValueError` czytelny zamiast błędu numpy).
> 2. podejście: 27/28 - ta nowa walidacja ujawniła, że ISTNIEJĄCY test
> `test_generate_synthetic_coverage_ksztalt` (sprzed tej sesji, length=2000) po
> cichu polegał na tym, że domyślny `duplication_region=(4200,4450)` POZA
> zakresem dla length=2000 wcześniej nie robił NIC (bo `rng.poisson(lam=tablica)`
> dopasowuje kształt do pustej tablicy, w przeciwieństwie do `deletion_region`,
> które twardo się wywala przy tej samej sytuacji) - asymetryczne, przypadkowe
> zachowanie, nie zamierzony kontrakt. Naprawione: ten test teraz jawnie przekazuje
> `deletion_region=None, duplication_region=None`. 3. podejście: **28/28 - wszystkie
> testy przechodzą.**

`tests/test_timdr_dna.py`: generowanie syntetycznych danych
(kształt, że delecja/duplikacja faktycznie mają inną głębokość),
`bin_coverage()` (agregacja, pusta seria), `analyze_coverage()`
integracyjnie - **wykrywa wstrzyknięte, znane z góry regiony delecji i
duplikacji** (test najważniejszy bez referencji: `test_analyze_coverage_wykrywa_wstrzykniete_delecje_i_duplikacje`),
za mało danych / zerowa mediana (czytelne błędy, nie cichy crash),
`load_depth_tsv()` (poprawny plik, wiele chromosomów = błąd, pusty plik
= błąd, brakujący plik = błąd - nigdy `except: pass`). Plus próbka
referencyjna: **`test_analyze_coverage_z_referencja_wykrywa_wiecej_niz_self_baseline_przy_ogromnej_delecji`**
(test najważniejszy dla tej funkcji - dowodzi na konkretnych liczbach,
że normalizacja względem referencji wykrywa istotnie więcej niż własna
mediana przy delecji obejmującej ~90% serii), plus błędna
długość/zerowa mediana referencji (czytelne błędy).

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
- Normalizacja względem próbki referencyjnej (patrz wyżej) NIE
  rozwiązuje problemu, jeśli sama referencja ma wariant w tym samym
  miejscu (np. wariant dziedziczny wspólny dla obu próbek) - to
  fundamentalne ograniczenie każdej metody porównawczej z referencją,
  nie coś specyficznego dla tej implementacji.
- Bez próbki referencyjnej (domyślny tryb) analiza nadal ma znany
  "self-baseline blind spot": wariant obejmujący całe albo prawie całe
  analizowane okno może zostać przeoczony, bo normalizuje się względem
  własnej (już przesuniętej) mediany. To jest teraz opcjonalnie
  adresowalne (`reference_depth=`), ale NIE jest to naprawione w trybie
  domyślnym - patrz sekcja o próbce referencyjnej wyżej.

## Licencja

MIT — patrz [LICENSE](LICENSE).
