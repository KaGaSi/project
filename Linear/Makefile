CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -O3 -ggdb
LIBS = -lm

SOURCE = ../AnalysisTools.c
OBJ = ../obj/AnalysisTools.o

ODIR = obj
#_OBJ = $(SOURCES:.c=.o)
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = bin/EndToEnd bin/PersistenceLength

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) ../obj bin

$(OBJ): $(OBJ) $(SOURCE)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

bin/EndToEnd: $(OBJ) $(ODIR)/EndToEnd.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/PersistenceLength: $(OBJ) $(ODIR)/PersistenceLength.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
