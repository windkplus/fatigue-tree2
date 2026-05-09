@echo off
REM ============================================================
REM   Extract-only run (re-uses already-computed .odb files)
REM   Update the cd path below if your .odb files are elsewhere.
REM ============================================================

REM cd /d "C:\Users\lcapk"
REM ^^^ if you are running this batch from a different folder,
REM     uncomment the line above with the folder containing the .odb files.

REM Delete stale CSV so we start fresh
if exist validation_results.csv del validation_results.csv

call abaqus python extract_validation.py val_rd50_e000
call abaqus python extract_validation.py val_rd50_e050
call abaqus python extract_validation.py val_rd50_e095
call abaqus python extract_validation.py val_rd70_e000
call abaqus python extract_validation.py val_rd70_e050
call abaqus python extract_validation.py val_rd70_e095

echo.
echo ============================================================
echo   Done. validation_results.csv should now have 6 rows.
echo   Bring this CSV back to the main project.
echo ============================================================
