CXXFLAGS += -Wall -pedantic -std=c++14 -g
LDFLAGS = -lUnitTest++

all: tests

FixedPointTests: FixedPointTests.o
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

FixedPointTests.o: FixedPointTests.cpp FixedPoint.hpp

tests: FixedPointTests
	ls /usr/include
	./FixedPointTests

.PHONY: clean tests

clean:
	rm -rf FixedPointTests.o FixedPointTests
