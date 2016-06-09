CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -ggdb
LIBS = -lm

SOURCE = ../AnalysisTools.c
OBJ = ../obj/AnalysisTools.o

ODIR = obj
#_OBJ = $(SOURCES:.c=.o)
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = ../bin/SelectedVcf \
			../bin/Aggregates \
			../bin/Aggregates-NotSameBeads \
			../bin/AggDensity \
			../bin/MolDensity \
			../bin/TransformVsf \
			../bin/Config \
			../bin/BondLength \
			../bin/DistrAgg \
			../bin/JoinAggregates \
			../bin/RgAgg \
			../bin/Average \
			../bin/traject-v2_5 \
			../bin/traject-v2_6

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR) ../obj ../bin

$(OBJ): $(OBJ) $(SOURCE)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

$(ODIR)/traject-v2_5.o:
	f95 -c -o $(ODIR)/traject-v2_5.o traject-v2_5.f90 -O3

$(ODIR)/traject-v2_6.o:
	f95 -c -o $(ODIR)/traject-v2_6.o traject-v2_6.f90 -O3

../bin/TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/Aggregates: $(OBJ) $(ODIR)/Aggregates.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/Aggregates-NotSameBeads: $(OBJ) $(ODIR)/Aggregates-NotSameBeads.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/AggDensity: $(OBJ) $(ODIR)/AggDensity.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/MolDensity: $(OBJ) $(ODIR)/MolDensity.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/Config: $(OBJ) $(ODIR)/Config.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/BondLength: $(OBJ) $(ODIR)/BondLength.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/JoinAggregates: $(OBJ) $(ODIR)/JoinAggregates.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/RgAgg: $(OBJ) $(ODIR)/RgAgg.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DistrAgg: $(OBJ) $(ODIR)/DistrAgg.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/Average: $(OBJ) $(ODIR)/Average.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/traject-v2_5: $(ODIR)/traject-v2_5.o
	f95 -o $@ traject-v2_5.f90 -O3

../bin/traject-v2_6: $(ODIR)/traject-v2_6.o
	f95 -o $@ traject-v2_6.f90 -O3

clean:
	rm -rf $(ODIR) *.o *~ core
