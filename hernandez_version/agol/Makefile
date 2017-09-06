CC = cc

CCFLAGS = -Wall -O3

all: nbody

universal.o: universal.h universal.c
	$(CC) $(CCFLAGS) -c -o universal.o universal.c

test-universal: test-universal.c universal.o
	$(CC) $(CCFLAGS) -o test-universal test-universal.c universal.o -lm

nbody: nbody.c universal.o
	$(CC) $(CCFLAGS) -o nbody nbody.c universal.o -lm	

clean:
	rm nbody universal.o 


