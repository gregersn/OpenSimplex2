CC=g++
CFLAGS=-g
OBJ = main.o OpenSimplex2s.o

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

demo: $(OBJ)
	$(CC) $(CFLAGS) -o demo main.o OpenSimplex2s.o

.PHONY: clean
clean:
	rm *.o demo
