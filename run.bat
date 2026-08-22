@echo off
title TIMDR-DNA -- API + Dashboard
color 0A
cls

:: Zawsze pracuj w katalogu, w ktorym faktycznie lezy ten plik .bat,
:: niezaleznie od tego, skad zostal uruchomiony (skrot, terminal, itp.)
cd /d "%~dp0"

echo ============================================================
echo   TIMDR-DNA: Uruchamianie API + Dashboard
echo   Katalog roboczy: %cd%
echo ============================================================
echo.

echo Sprawdzanie, czy port 8070 jest juz zajety...
netstat -ano | findstr ":8070" >nul
if %ERRORLEVEL% EQU 0 (
    echo.
    echo [UWAGA] Port 8070 jest JUZ zajety przez inny proces - to on
    echo         bedzie odpowiadal w przegladarce, nie ten skrypt:
    netstat -ano | findstr ":8070"
    echo         Ostatnia liczba w linii to PID procesu. Sprawdz go w
    echo         Menedzerze zadan i zamknij, albo zmien port w api.py.
    echo.
    pause
)

if exist "venv\Scripts\activate.bat" (
    echo [OK] Aktywacja venv...
    call venv\Scripts\activate.bat
) else if exist ".venv\Scripts\activate.bat" (
    echo [OK] Aktywacja .venv...
    call .venv\Scripts\activate.bat
) else (
    echo [INFO] Uzywanie systemowej instalacji Pythona.
)

echo.
echo [1/2] Weryfikacja i instalacja pakietow pip...
python -m pip install --upgrade pip --disable-pip-version-check
python -m pip install -r requirements.txt

if %ERRORLEVEL% NEQ 0 (
    echo [BLAD] Nie udalo sie zainstalowac wymaganych pakietow.
    pause
    exit /b 1
)

if not exist "api.py" (
    echo [BLAD] Nie znaleziono pliku "%cd%\api.py"
    pause
    exit /b 1
)

echo.
echo [2/2] Uruchamianie API + dashboard pod http://127.0.0.1:8070 ...
echo Przegladarka otworzy sie sama za chwile (moze zajac 1-2 sekundy).
echo Nacisnij CTRL+C, aby zatrzymac serwer.
echo ---------------------------------------------------
echo.

python api.py

echo.
echo ============================================================
echo Serwer zostal zamkniety.
echo ============================================================
pause
