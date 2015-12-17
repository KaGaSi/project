CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
LIBS = -lm

SOURCES = Aux.c Structure.c

DEPS = Aux.h Structure.h CStructs.h

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = SelectedVcf TransformVsf

all: dir $(BIN) $(OBJ)

dir:
	mkdir -p $(ODIR)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#copy:
#	cp $(BIN) bin

clean:
	rm -rf $(ODIR) *.o *~ core
