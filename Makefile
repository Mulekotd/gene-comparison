CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -pthread
CPPFLAGS ?= -Isrc

.PHONY: all test run run-panel extract clean

all: build/gene-comparison

build:
	mkdir -p build

build/gene-comparison: src/main.cpp src/cli.cpp src/cli.hpp src/fasta.cpp src/fasta.hpp src/report.cpp src/report.hpp src/utils.cpp src/utils.hpp src/alignment.cpp src/alignment.hpp src/multiple_alignment.cpp src/multiple_alignment.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/main.cpp src/cli.cpp src/fasta.cpp src/report.cpp src/utils.cpp src/alignment.cpp src/multiple_alignment.cpp -o $@

build/test-alignment: tests/test_alignment.cpp src/alignment.cpp src/alignment.hpp src/multiple_alignment.cpp src/multiple_alignment.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/test_alignment.cpp src/alignment.cpp src/multiple_alignment.cpp -o $@

test: build/test-alignment build/gene-comparison
	./build/test-alignment
	python3 tests/test_integration.py

run: build/gene-comparison
	./build/gene-comparison --benchmark

extract:
	python3 scripts/extract_proteins.py

run-panel: build/gene-comparison
	./build/gene-comparison --protein recA --benchmark
	./build/gene-comparison --protein trhO --benchmark
	./build/gene-comparison --protein rpoB --benchmark

clean:
	rm -rf build
