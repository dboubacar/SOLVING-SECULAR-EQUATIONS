CC=gcc
CFLAGS= -g -O3 -Wall
LDFLAGS=
EXEC=solv_secular
LIB= -lm -fopenmp


all: $(EXEC)

solv_secular: param.o gragg.o main.o hybrid.o
	$(CC) -o $@ $^ $(LDFLAGS)  $(LIB)

main.o: main.c gragg.h param.h hybrid.h
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)

gragg.o: gragg.c param.h
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)

hybrid.o: hybrid.c param.h
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)

param.o: param.c
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)


clean:
	rm -rf *.o solv_secular

mrproper: clean
	rm -rf $(EXEC)
