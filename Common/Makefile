CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -O3 -ggdb
LIBS = -lm

SOURCES = ../AnalysisTools.c

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = bin/Aggregates bin/SelectedVcf bin/TransformVsf bin/Config

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) bin

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

bin/TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/Aggregates: $(OBJ) $(ODIR)/Aggregates.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/Config: $(OBJ) $(ODIR)/Config.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
