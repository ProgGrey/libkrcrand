CXX=g++
CXXFLAGS=-std=c++17 -O2 -march=x86-64-v3  -flto -Wall -Wextra -Wpedantic

all:test speed

speed:libkrcrand.a tests/speed.cpp
	$(CXX) ${CXXFLAGS} tests/speed.cpp libkrcrand.a -o speed

test:libkrcrand.a tests.o
	$(CXX) ${CXXFLAGS} tests.o libkrcrand.a -lboost_unit_test_framework -o test

tests.o:tests/tests.cpp headers/generators.hpp
	$(CXX) ${CXXFLAGS} -c tests/tests.cpp -o tests.o

libkrcrand.a:tools.o Xoshiro256mm.o
	ar crf libkrcrand.a tools.o Xoshiro256mm.o

tools.o:src/tools.cpp
	$(CXX) ${CXXFLAGS} -c src/tools.cpp -o tools.o

Xoshiro256mm.o:src/generators/Xoshiro256mm.cpp headers/generators.hpp
	$(CXX) ${CXXFLAGS} -c src/generators/Xoshiro256mm.cpp -o Xoshiro256mm.o

check:
	cppcheck --enable=all ./src/* ./headers/* ./tests/*

clean:
	rm *.o *.a test speed
