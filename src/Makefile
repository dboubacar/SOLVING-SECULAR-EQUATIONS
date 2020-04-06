CC=gcc
CFLAGS= -g -O3 -Wall
LDFLAGS=
EXEC=gragg
LIB= -lm


all: $(EXEC)

gragg: param.o gragg.o main.o
	$(CC) -o $@ $^ $(LDFLAGS)  $(LIB)

main.o: main.c gragg.h param.h
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)

gragg.o: gragg.c param.h
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)


param.o: param.c
	$(CC) -o $@ -c $< $(CFLAGS) $(LIB)


clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
