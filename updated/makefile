CC=g++
CFLAGS=-g
DEPS = Dispersion.h Coord_Num.h
OBJ = Main.o Dispersion.o Coord_Num.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MP2D: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY: clean

clean: 
	rm -f *.o MP2D
