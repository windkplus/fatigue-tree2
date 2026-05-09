@echo off
echo ============================================================
echo   ABAQUS 2D Validation Run (paper #2 Appendix A)
echo ============================================================
cd /d "C:/ClaudeCode_Projects/fatigue-tree2/data/abaqus_validation"

echo Running val_rd50_e000
call abaqus job=val_rd50_e000 cpus=1 interactive ask_delete=off
echo Running val_rd50_e050
call abaqus job=val_rd50_e050 cpus=1 interactive ask_delete=off
echo Running val_rd50_e095
call abaqus job=val_rd50_e095 cpus=1 interactive ask_delete=off
echo Running val_rd70_e000
call abaqus job=val_rd70_e000 cpus=1 interactive ask_delete=off
echo Running val_rd70_e050
call abaqus job=val_rd70_e050 cpus=1 interactive ask_delete=off
echo Running val_rd70_e095
call abaqus job=val_rd70_e095 cpus=1 interactive ask_delete=off

echo.
echo ============================================================
echo   Extracting mass properties from each .odb
echo ============================================================

call abaqus python extract_validation.py val_rd50_e000
call abaqus python extract_validation.py val_rd50_e050
call abaqus python extract_validation.py val_rd50_e095
call abaqus python extract_validation.py val_rd70_e000
call abaqus python extract_validation.py val_rd70_e050
call abaqus python extract_validation.py val_rd70_e095

echo.
echo ============================================================
echo   Validation runs complete. See validation_results.csv
echo ============================================================
