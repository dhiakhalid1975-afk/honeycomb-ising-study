@echo off
setlocal
cd /d "%~dp0"
set "PY=python"
if exist PYTHON_EXE.local.txt set /p PY=<PYTHON_EXE.local.txt
"%PY%" run_final_hardening.py all --config USER_CONFIG.json 
if errorlevel 1 (
  echo.
  echo FAIL-CLOSED. Do not continue.
  exit /b 1
)
echo.
echo PASS: all
endlocal
