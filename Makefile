OBJECTS=
CFLAGS= -std=c99 -Wall -g -Wall -O3
LDLIBS= -lm -lpthread
CC= gcc

all: plutoindex

plutoindex: $(OBJECTS) indexgen.c
	$(CC) $(CFLAGS) $(OBJECTS) indexgen.c $(LDLIBS) -o plutoindex

clean:
	rm -f $(OBJECTS) plutoindex
