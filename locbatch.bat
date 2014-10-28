@echo off

rem find the path for Anaconda
for /f "delims=" %%a in ('where Anaconda') do @set TMP=%%~dpa

rem run the actual command
echo Running: "%TMP%..\python.exe" "%~dp0\locbatch\locbatch.py" -c %*
"%TMP%..\python.exe" "%~dp0\locbatch_code\locbatch.py" -c %*

pause