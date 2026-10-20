@echo off
rem Runs make_fixtures.R and make_fixtures_ext.R with R 4.5.2;
rem output in make_fixtures.log and make_fixtures_ext.log
cd /d "%~dp0"
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" make_fixtures.R > make_fixtures.log 2>&1
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" make_fixtures_ext.R > make_fixtures_ext.log 2>&1
