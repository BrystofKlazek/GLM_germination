@echo off
title GLM Germination Analysis

echo.
echo  ============================================================
echo     GLM Germination Analysis - Launcher
echo  ============================================================
echo.

:: ── Find R ──
where Rscript.exe >nul 2>&1
if not errorlevel 1 goto :use_path

dir /b "C:\Program Files\R\R-*" >nul 2>&1
if errorlevel 1 goto :no_r

for /d %%d in ("C:\Program Files\R\R-*") do set "RSCRIPT=%%d\bin\Rscript.exe"
if exist "%RSCRIPT%" goto :found_r

for /d %%d in ("C:\Program Files\R\R-*") do set "RSCRIPT=%%d\bin\x64\Rscript.exe"
if exist "%RSCRIPT%" goto :found_r

:no_r
echo  R nebylo nalezeno. Nainstalujte z https://cran.r-project.org/
pause
exit /b 1

:use_path
for /f "delims=" %%i in ('where Rscript.exe') do set "RSCRIPT=%%i"

:found_r
echo  R: %RSCRIPT%

:: ── Check app.R ──
if not exist "%~dp0app.R" goto :no_app
echo  app.R nalezen
echo.
goto :run

:no_app
echo  app.R nenalezen! Dejte ho do stejne slozky jako tento soubor.
pause
exit /b 1

:run
:: Convert backslashes to forward slashes for R
set "APPDIR=%~dp0"
set "APPDIR=%APPDIR:\=/%"

echo  Spoustim (pri prvnim spusteni se instaluji balicky, muze to trvat)...
echo.
"%RSCRIPT%" --no-save --no-restore -e "setwd('%APPDIR%');if(!require('shiny',quietly=TRUE)){install.packages(c('shiny','bslib','readxl','readr','emmeans','ggplot2','dplyr','tidyr','broom','DT','multcomp','multcompView','car','scales'),repos='https://cloud.r-project.org')};shiny::runApp('app.R',launch.browser=TRUE)"

echo.
echo  Aplikace ukoncena.
pause
