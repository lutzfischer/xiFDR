SET mypath=%~dp0
java -jar %mypath%xiFDR-1.0.22-jar-with-dependencies.jar
if NOT ["%errorlevel%"]==["0"] (
    pause
    exit /b %errorlevel%
)
exit
