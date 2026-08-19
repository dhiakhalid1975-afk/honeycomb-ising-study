@echo off
setlocal
cd /d "%~dp0"
set "PY=python"
if exist PYTHON_EXE.local.txt set /p PY=<PYTHON_EXE.local.txt
"%PY%" -m pytest -q tests
if errorlevel 1 exit /b 1
echo PACKAGE TESTS PASS.
endlocal
