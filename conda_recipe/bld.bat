:: BPM solver 

if exist bpm\bin rd /s /q bpm\bin
mkdir bpm\bin
cd bpm\bin

cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ..
if errorlevel 1 exit 1

nmake install python-install
if errorlevel 1 exit 1

cd ..\..


:: RTM solver

if exist rtm\bin rd /s /q rtm\bin
mkdir rtm\bin
cd rtm\bin

cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ..
if errorlevel 1 exit 1

nmake install
if errorlevel 1 exit 1

cd ../..



:: DTMM solver 

cd dtmm

"%PYTHON%" setup.py install
if errorlevel 1 exit 1

cd ..



:: High-level interface

cd nemaktis

"%PYTHON%" setup.py install
if errorlevel 1 exit 1
