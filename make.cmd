cl /O2 /GL /std:c++latest /EHs /fp:fast /arch:AVX2 /I lib\lodepng src\buddhabrot.cpp lib\lodepng\lodepng.cpp /Fo:bin\ /link /out:bin\buddhabrot.exe
