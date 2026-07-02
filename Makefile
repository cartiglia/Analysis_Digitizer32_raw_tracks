# Makefile for the Digitizer32 raw-tracks analysis
#
#   make        -> builds Analysis_Digitizer32_raw_tracks_flt
#   make clean  -> removes the executable

CXX      := clang++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra $(shell root-config --cflags)
LDFLAGS  := $(shell root-config --glibs)

SRC    := Analysis_Digitizer32_raw_tracks_flt.cpp
TARGET := Analysis_Digitizer32_raw_tracks_flt

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC) Analysis_Digitizer16.h
	$(CXX) $(CXXFLAGS) $(SRC) -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o *~
