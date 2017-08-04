UNITTEST_CPP_INCLUDE ?= /usr/include/UnitTest++

ifneq (,$(findstring g++,$(CXX)))
	CXXFLAGS += --coverage
endif

CXXFLAGS += -Wall -pedantic -std=c++14 -g -I$(UNITTEST_CPP_INCLUDE)
LDFLAGS = -lUnitTest++ -lgcov

all: tests

FixedPointTests: FixedPointTests.o
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

FixedPointTests.o: FixedPointTests.cpp FixedPoint.hpp

tests: FixedPointTests
	./FixedPointTests

coverage: FixedPoint.hpp.gcov

FixedPoint.hpp.gcov: FixedPointTests
	rm -rf *.gcov *.gcda
	./FixedPointTests
	gcov FixedPointTests.gcda

.PHONY: clean tests coverage

clean:
	rm -rf FixedPointTests.o FixedPointTests *.gcov *.gcda *.gcno
