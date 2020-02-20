CC = g++
CFLAGS = -v -Wall -Wextra  -Werror -g

lattice : main.o
	$(CC) $(CFLAGS) main.o -o lattice

polymer : mainPolymer.o
	$(CC) $(CFLAGS) mainPolymer.o -o polymer

test : polymerTest.o
	$(CC) $(CFLAGS) polymerTest.o -o test

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o lattice polymer test
