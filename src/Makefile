CC= gcc
CFLAGS= -O3 -Wall -Wextra -g
EXEC= solv_secular
OBJ= param.o gragg.o main.o hybrid.o
LIB= -lm -fopenmp

all: $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -o $(EXEC) $(LIB)

gragg.o: gragg.c
	$(CC) -c $(CFLAGS) $< $(LIB)

hybrid.o: hybrid.c
		$(CC) -c $(CFLAGS) $< $(LIB)

param.o:param.c param.h
	$(CC) -c $(CFLAGS) param.c  $(LIB)


.PHONY: clean

clean :
	rm -rf $(EXEC) *~ *.o
