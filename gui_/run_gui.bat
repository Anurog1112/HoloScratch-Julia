@echo off
echo ========================================
echo   HoloScratch GUI - Starting...
echo ========================================
set JULIA_DEPOT_PATH=D:\Julia-Depot
set PATH=%PATH%;C:\Users\Anurag Namgat\AppData\Local\Programs\Julia-1.12.4\bin
cd /d D:\HoloScratch-Julia
julia --project=. holoscratch_gui.jl
pause
