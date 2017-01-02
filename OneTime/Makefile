CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
#CFLAGS = -std=c11 -Wall -Werror -ggdb
LIBS = -lm

SOURCE = ../AnalysisTools.c
OBJ = ../obj/AnalysisTools.o

ODIR = obj
#_OBJ = $(SOURCES:.c=.o)
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = ../bin/ChargeRatio \
			../bin/DistrAgg2D \
			../bin/AggDensity-CICheck \
			../bin/DistrAgg-NumberOfAbsorbedHomopolymerChains \
			../bin/GyrationAggregates-NotHomopol \
			../bin/GyrationAggregates-NotHomopolNotB \
			../bin/GyrationAggregates-NotBBlock \
			../bin/DensityAggregates-NotHomop \
			../bin/DensityMolecules-NearNanoSurf \
			../bin/DistrAgg-NotHomopol-NumberOfAbsorbedHomopolymerChains \
			../bin/DistrAgg-OnlyPols \
			../bin/NanoparticleDistances

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

../bin/DistrAgg-NumberOfAbsorbedHomopolymerChains: $(OBJ) $(ODIR)/DistrAgg-NumberOfAbsorbedHomopolymerChains.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/GyrationAggregates-NotHomopol: $(OBJ) $(ODIR)/GyrationAggregates-NotHomopol.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/GyrationAggregates-NotHomopolNotB: $(OBJ) $(ODIR)/GyrationAggregates-NotHomopolNotB.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/GyrationAggregates-NotBBlock: $(OBJ) $(ODIR)/GyrationAggregates-NotBBlock.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DensityAggregates-NotHomop: $(OBJ) $(ODIR)/DensityAggregates-NotHomop.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DensityMolecules-NearNanoSurf: $(OBJ) $(ODIR)/DensityMolecules-NearNanoSurf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DistrAgg-NotHomopol-NumberOfAbsorbedHomopolymerChains: $(OBJ) $(ODIR)/DistrAgg-NotHomopol-NumberOfAbsorbedHomopolymerChains.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/DistrAgg-OnlyPols: $(OBJ) $(ODIR)/DistrAgg-OnlyPols.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

../bin/NanoparticleDistances: $(OBJ) $(ODIR)/NanoparticleDistances.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf $(ODIR) *.o *~ core
