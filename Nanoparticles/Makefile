CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -O3 -ggdb
LIBS = -lm

#SOURCES = ../AnalysisTools.c

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = bin/C60

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) bin

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

bin/C60: $(OBJ) $(ODIR)/C60.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
