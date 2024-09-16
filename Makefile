CXX=g++
AR=ar
CXXFLAGS_BASE= -std=c++20 -O2 -march=native -flto=auto -Wno-maybe-uninitialized
CXXFLAGS= ${CXXFLAGS_BASE} -Wall -Wextra -pedantic -Wno-unused-parameter

all:test speed gen gamma

gamma:libkrcrand.a tests/gamma.cpp  headers/gammaDistribution.hpp
	$(CXX) ${CXXFLAGS} tests/gamma.cpp libkrcrand.a -o gamma

gen:libkrcrand.a tests/gen.cpp
	$(CXX) ${CXXFLAGS} tests/gen.cpp libkrcrand.a -o gen

speed:libkrcrand.a tests/speed.cpp headers/exponentialDistribution.hpp headers/gammaDistribution.hpp
	$(CXX) ${CXXFLAGS} tests/speed.cpp libkrcrand.a -o speed

test:libkrcrand.a tests.o
	$(CXX) ${CXXFLAGS} tests.o libkrcrand.a -lboost_unit_test_framework -o test

tests.o:tests/tests.cpp
	$(CXX) ${CXXFLAGS} -c tests/tests.cpp -o tests.o

libkrcrand.a:tools.o Xoshiro256mm.o math_generic.o math_sse2.o math_avx2.o math_avx512.o gammaDistribution.o cdflib.o functions.o 
	$(AR) rc libkrcrand.a tools.o Xoshiro256mm.o math_generic.o math_sse2.o math_avx2.o math_avx512.o gammaDistribution.o cdflib.o functions.o

tools.o:src/tools.cpp src/tools.hpp
	$(CXX) ${CXXFLAGS} -c src/tools.cpp -o tools.o

math_generic.o:src/math/generic.cpp headers/math.hpp
	$(CXX) ${CXXFLAGS} -c src/math/generic.cpp -o math_generic.o

math_sse2.o:src/math/sse2.cpp headers/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/sse2.cpp -o math_sse2.o

math_avx2.o:src/math/avx2.cpp headers/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/avx2.cpp -o math_avx2.o

math_avx512.o:src/math/avx512.cpp headers/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/avx512.cpp -o math_avx512.o

Xoshiro256mm.o:src/generators/Xoshiro256mm.cpp headers/generators.hpp
	$(CXX) ${CXXFLAGS} -c src/generators/Xoshiro256mm.cpp -o Xoshiro256mm.o

gammaDistribution.o:src/distributions/gammaDistribution.cpp headers/math.hpp headers/gammaDistribution.hpp headers/linearlyApproximatedDistribution.hpp
	$(CXX) ${CXXFLAGS} -c src/distributions/gammaDistribution.cpp -o gammaDistribution.o

cdflib.o:src/3rdparty/cdflib.cpp src/3rdparty/cdflib.hpp src/3rdparty/raw/cdflib.cpp src/3rdparty/raw/cdflib.hpp
	$(CXX) ${CXXFLAGS_BASE} -Wno-maybe-uninitialized -c src/3rdparty/cdflib.cpp -o cdflib.o -Wno-maybe-uninitialized

functions.o:src/math/functions.cpp headers/math.hpp
	$(CXX) ${CXXFLAGS} -c src/math/functions.cpp -o functions.o

check:
	cppcheck --enable=all ./src/* ./headers/* ./tests/*

clean:
	rm *.o *.a test speed gen
