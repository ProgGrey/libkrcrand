CXX=clang++
CXXFLAGS=-std=c++17 -O3 -march=native  -flto -Wall -Wextra -Wpedantic

all:test speed gen

gen:libkrcrand.a tests/gen.cpp
	$(CXX) ${CXXFLAGS} tests/gen.cpp libkrcrand.a -o gen

speed:libkrcrand.a tests/speed.cpp headers/exponentialDistribution.hpp
	$(CXX) ${CXXFLAGS} tests/speed.cpp libkrcrand.a -o speed

test:libkrcrand.a tests.o
	$(CXX) ${CXXFLAGS} tests.o libkrcrand.a -lboost_unit_test_framework -o test

tests.o:tests/tests.cpp headers/generators.hpp
	$(CXX) ${CXXFLAGS} -c tests/tests.cpp -o tests.o

libkrcrand.a:tools.o Xoshiro256mm.o math_generic.o math_sse2.o math_avx2.o math_avx512.o
	ar crf libkrcrand.a tools.o Xoshiro256mm.o math_generic.o math_sse2.o math_avx2.o math_avx512.o

tools.o:src/tools.cpp src/tools.hpp
	$(CXX) ${CXXFLAGS} -c src/tools.cpp -o tools.o

math_generic.o:src/math/generic.cpp src/math.hpp
	$(CXX) ${CXXFLAGS} -c src/math/generic.cpp -o math_generic.o

math_sse2.o:src/math/sse2.cpp src/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/sse2.cpp -o math_sse2.o

math_avx2.o:src/math/avx2.cpp src/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/avx2.cpp -o math_avx2.o

math_avx512.o:src/math/avx512.cpp src/math.hpp src/math/amd64.hpp
	$(CXX) ${CXXFLAGS} -c src/math/avx512.cpp -o math_avx512.o

Xoshiro256mm.o:src/generators/Xoshiro256mm.cpp headers/generators.hpp
	$(CXX) ${CXXFLAGS} -c src/generators/Xoshiro256mm.cpp -o Xoshiro256mm.o

check:
	cppcheck --enable=all ./src/* ./headers/* ./tests/*

clean:
	rm *.o *.a test speed gen
