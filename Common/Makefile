CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -O3 -ggdb
LIBS = -lm

SOURCE = ../AnalysisTools.c
OBJ = ../obj/AnalysisTools.o

ODIR = obj
#_OBJ = $(SOURCES:.c=.o)
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = bin/Aggregates bin/SelectedVcf bin/TransformVsf bin/Config bin/traject

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) ../obj bin

$(OBJ): $(OBJ) $(SOURCE)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

bin/TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/Aggregates: $(OBJ) $(ODIR)/Aggregates.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/Config: $(OBJ) $(ODIR)/Config.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/traject.o:
	f95 -c -o $(ODIR)/traject.o traject.f90 -O3

bin/traject: $(ODIR)/traject.o
	f95 -o $@ traject.f90 -O3

clean:
	rm -rf $(ODIR) *.o *~ core
