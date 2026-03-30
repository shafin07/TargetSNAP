@echo off
cd /d "%~dp0"
echo Installing Flask...
pip install flask flask-cors
echo.
echo Starting TargetSNAP Web Server...
echo Open: http://localhost:5000
echo.
python app.py
pause
