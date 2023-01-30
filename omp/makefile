CC= g++
CFLAGS=-fopenmp -lm -O2
SOURCE= src/*.C
OBJECTS=$(SOURCE:.c=.o)
EXECUTABLE=solver
all: $(EXECUTABLE)
$(EXECUTABLE):$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(CFLAGS)
.c.o:
	$(CC) $(CFLAGS) $< 
.PHONY: all clean
clean:
	rm -rf *o solver
