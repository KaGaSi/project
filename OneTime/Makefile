CC = gcc
#CFLAGS = -std=c11 -Wall -Werror -O3
CFLAGS = -std=c11 -Wall -Werror -ggdb
LIBS = -lm

SOURCE = ../AnalysisTools.c
OBJ = ../obj/AnalysisTools.o

ODIR = obj
#_OBJ = $(SOURCES:.c=.o)
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = ../bin/ChargeRatio \
			../bin/DistrAgg2D \
			../bin/AggDensity-CICheck \
			../bin/DistrAgg-NoNanoAggs

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) ../obj ../bin

$(OBJ): $(OBJ) $(SOURCE)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

../bin/ChargeRatio: $(OBJ) $(ODIR)/ChargeRatio.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DistrAgg2D: $(OBJ) $(ODIR)/DistrAgg2D.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/AggDensity-CICheck: $(OBJ) $(ODIR)/AggDensity-CICheck.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DistrAgg-NoNanoAggs: $(OBJ) $(ODIR)/DistrAgg-NoNanoAggs.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
