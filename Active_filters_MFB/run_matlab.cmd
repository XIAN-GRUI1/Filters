@echo off
rem =====================================================================
rem  run_matlab.cmd - start MATLAB under the DSH file sandbox
rem
rem  The sandbox (workspace-write) blocks MATLAB's startup because it
rem  writes outside the writable areas (preferences, AppData, temp ...).
rem  This wrapper redirects ALL of MATLAB's writable locations into the
rem  local folder ".matlab_sandbox" (which lies inside the writable
rem  session workspace), so MATLAB starts normally.
rem
rem  Usage (from this folder):
rem     run_matlab.cmd                     -> interactive MATLAB (desktop)
rem     run_matlab.cmd -batch "sktest"     -> run a MATLAB batch command
rem     run_matlab.cmd -batch "skDemo"
rem     run_matlab.cmd -r "skDesign('Type','cauer','Order',5,'Fc',1e3)"
rem
rem  Any argument is passed through to matlab.exe.
rem =====================================================================
setlocal
set "BASE=%~dp0.matlab_sandbox"
if not exist "%BASE%"                 mkdir "%BASE%"
set "PROF=%BASE%\user"
if not exist "%PROF%"                 mkdir "%PROF%"
if not exist "%PROF%\AppData\Roaming" mkdir "%PROF%\AppData\Roaming"
if not exist "%PROF%\AppData\Local"   mkdir "%PROF%\AppData\Local"
if not exist "%BASE%\tmp"             mkdir "%BASE%\tmp"

set "USERPROFILE=%PROF%"
set "HOMEDRIVE=%~d0"
set "HOMEPATH=%~p0.matlab_sandbox\user"
set "APPDATA=%PROF%\AppData\Roaming"
set "LOCALAPPDATA=%PROF%\AppData\Local"
set "MATLAB_PREFDIR=%PROF%\mlpref"
set "TEMP=%BASE%\tmp"
set "TMP=%BASE%\tmp"

cd /d "%~dp0"
matlab %*
endlocal
