CXX := g++
CXXFLAGS := -std=c++17 -O3
NVFLAGS := $(CXXFLAGS)
OMPFLAGS := -std=c++17 -O3 -fopenmp
TARGET := final 
SEQUENTIAL := sequential 
TCGEN := testcase_generator
BFS := bfs
SWEEP := sweep
EXEDIR = ../bin

.PHONY: all
all: $(TARGET)

.PHONY: seq 
seq: sequential.cc
	$(CXX) $(CXXFLAGS) -o $(EXEDIR)/sequential sequential.cc

.PHONY: bfs 
bfs: bfs.cc
	$(CXX) $(CXXFLAGS) -o $(EXEDIR)/bfs bfs.cc

.PHONY: sweep 
sweep: sweep.cc
	$(CXX) $(CXXFLAGS) -o $(EXEDIR)/sweep sweep.cc

.PHONY: tcgen
tcgen: testcase_generator.cc
	$(CXX) $(CXXFLAGS) -o $(EXEDIR)/testcase_generator testcase_generator.cc

.PHONY: cps 
cps: sweepCPS.cu
	nvcc $(NVFLAGS) -o $(EXEDIR)/cps sweepCPS.cu

.PHONY: ps 
ps: sweepWithPS.cu
	nvcc $(NVFLAGS) -o $(EXEDIR)/ps sweepWithPS.cu

.PHONY: cpsomp
cpsomp: sweepCPS_omp.cc
	$(CXX) $(OMPFLAGS) -o $(EXEDIR)/cpsomp sweepCPS_omp.cc

.PHONY: psomp
ps: sweepWithPS_omp.cc
	$(CXX) $(OMPFLAGS) -o $(EXEDIR)/psomp sweepWithPS_omp.cc

.PHONY: clean
clean:
	rm -f $(EXEDIR)/$(TARGET) $(EXEDIR)/$(SEQUENTIAL) $(EXEDIR)/$(TCGEN) $(EXEDIR)/$(BFS) $(EXEDIR)/$(SWEEP)
