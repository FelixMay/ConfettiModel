CXX=g++

OBJS = \
	Main.o \
	Forest.o \
	Square.o \
	stoc1.o \
	mersenne.o \
	userintf.o

default: install

%.o: %.cpp
	$(CXX) -c -std=c++11 -g -O3 -Wall -Wno-sign-compare -o $@ $<

ConfettiSquare: $(OBJS)
	$(CXX) -static -o ConfettiSquare $(OBJS)
	strip -s ConfettiSquare

clean:
	rm -f *.o
	rm -f ConfettiSquare

install: ConfettiSquare
	mkdir -p $(HOME)/bin
	cp ConfettiSquare $(HOME)/bin

.PHONY: clean default install
