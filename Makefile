CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3 -ggdb
LIBS = -lm

SOURCES = AnalysisTools.c

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = Aggregates SelectedVcf TransformVsf

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR)

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Aggregates: $(OBJ) $(ODIR)/Aggregates.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
