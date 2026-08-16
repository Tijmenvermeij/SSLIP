@echo off
echo ========================================================
echo Pushing SSLIP Updates to GitHub
echo ========================================================

:: Change directory to where the batch file is located
cd /d "%~dp0"

:: Check git status
git status
echo.

:: Stage all changes
git add -A

:: Prompt user for commit message
set /p commit_msg="Enter commit message (or press Enter for default 'Update'): "
if "%commit_msg%"=="" set commit_msg=Update

:: Commit and Push
git commit -m "%commit_msg%"
git push origin main

echo.
echo ========================================================
echo Push complete!
pause
