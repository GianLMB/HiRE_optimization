CC=g++
cfiles = $(wildcard *.c)
cppfiles = $(wildcard *.cpp)

MDDRIVER=../mddriver/protocol
LIBS= -I$(MDDRIVER)/include -L$(MDDRIVER)/lib

all: Interactor.o InteractorMDDriver.o OPEP_wrapper.o

%.o: %.c
	$(CC) -c  -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CC) -c  -o $@ $^ $(LIBS)

clean:
	@rm -rf *.o
