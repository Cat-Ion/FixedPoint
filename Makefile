CXXFLAGS += -Wall -pedantic -std=c++17
LDFLAGS = -lcpptest

all: tests

FixedPointTests: FixedPointTests.o
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

FixedPointTests.o: FixedPointTests.cpp FixedPoint.hpp

tests: FixedPointTests
	./FixedPointTests

.PHONY: clean

clean:
	rm -rf FixedPointTests.o FixedPointTests
