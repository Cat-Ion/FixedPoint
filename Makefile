UNITTEST_CPP_INCLUDE ?= /usr/include/UnitTest++

CXXFLAGS += -Wall -pedantic -std=c++14 -g -I$(UNITTEST_CPP_INCLUDE)
LDFLAGS = -lUnitTest++


all: tests

FixedPointTests: FixedPointTests.o
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

FixedPointTests.o: FixedPointTests.cpp FixedPoint.hpp

tests: FixedPointTests
	./FixedPointTests

.PHONY: clean tests

clean:
	rm -rf FixedPointTests.o FixedPointTests
