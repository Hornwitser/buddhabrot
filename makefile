CXXFLAGS = -Wall -Wextra -Wpedantic -std=c++20 -march=native -Ilib/lodepng
HEADERS = lib/lodepng/lodepng.h src/parse_args.hpp src/srgb.hpp
UNITS = lib/lodepng/lodepng.cpp src/buddhabrot.cpp
SOURCES = ${HEADERS} ${UNITS}


run: bin/buddhabrot
	./bin/buddhabrot

debug: bin/buddhabrot-debug
	gdb --args ./buddhabrot-debug

bin:
	mkdir bin

bin/buddhabrot: ${SOURCES} | bin
	g++ ${CXXFLAGS} -O3 ${UNITS} -o bin/buddhabrot

bin/buddhabrot-debug: ${SOURCES} | bin
	g++ ${CXXFLAGS} -g -O0 ${UNITS} -o bin/buddhabrot-debug
