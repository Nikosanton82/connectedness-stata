@echo off
rem Runs make_fixtures.R with R 4.5.2; output in make_fixtures.log
cd /d "%~dp0"
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" make_fixtures.R > make_fixtures.log 2>&1
