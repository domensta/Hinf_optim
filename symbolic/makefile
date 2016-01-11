SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXX = g++
LFLAGS = -lcln -lginac -std=c++11 -I/user/local/lib

all: $(BINS)

% :	%.cpp
	echo $(BINS) $(OBJS)
	$(CXX) $< -o $@ $(LFLAGS)

clean:
	rm -f $(BINS)

	
