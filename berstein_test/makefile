


SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXX = g++
CXXFLAGS := $(shell pkg-config --cflags /home/bumquist/Tools/ibex/OUT/share/pkgconfig/ibex.pc) 
LIBS	 := $(shell pkg-config --libs  /home/bumquist/Tools/ibex/OUT/share/pkgconfig/ibex.pc)
#LIBDIR	 := $(shell pkg-config --libdir  ibex) 
LFLAGS = -lcln -lginac -std=c++11 -I/user/local/lib


ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall -frounding-math -std=c++11
else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-deprecated -frounding-math -std=c++11
endif

all: $(BINS)
      
% :	%.cpp
	echo $(BINS) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS) $(LFLAGS)

clean:
	rm -f $(BINS)

	
