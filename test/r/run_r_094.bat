@echo off
rem Runs make_fixtures_094.R with R 4.5.2; output in make_fixtures_094.log
cd /d "%~dp0"
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" make_fixtures_094.R > make_fixtures_094.log 2>&1
