CC = c++
CFLAGS = -g -Wall -Wno-deprecated 
OBJECTS =  main.o distR.o matrixOperations.o fileio.o distanceList.o
OBJECTS2 = wrapper.o phylogeny.o global.o bit_set.o simple_nexus.o
INCDIR = -I../include 
PROGNAME = distR

all: main calcdistR matrixops fileio distancelist wrapper phylogeny global bitset simplenexus
	$(CC) $(CFLAGS) $(INCDIR) -o $(PROGNAME) $(OBJECTS) $(OBJECTS2) -lm

main: main.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c main.cpp

calcdistR: distR.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c distR.cpp

matrixops: matrixOperations.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c matrixOperations.cpp

fileio: fileio.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c fileio.cpp

distancelist: distanceList.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c distanceList.cpp

wrapper:   wrapper.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c wrapper.cpp

phylogeny:  phylogeny.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c phylogeny.cpp

global:  global.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c global.cpp

bitset: bit_set.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c bit_set.cpp

simplenexus: simple_nexus.cpp
	$(CC) $(CFLAGS) $(INCDIR) -c simple_nexus.cpp

clean:
	rm -f *.o *~ distR
