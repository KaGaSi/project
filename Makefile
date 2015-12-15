CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3 -g
LIBS = -lm

SOURCES = Aux.c Structure.c

DEPS = Aux.h Structure.h CStructs.h

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = SelectedVcf TransformVsf

all: $(BIN) $(OBJ)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TransformVsf: $(OBJ) $(ODIR)/TransformVsf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

SelectedVcf: $(OBJ) $(ODIR)/SelectedVcf.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#copy:
#	cp $(BIN) bin

clean:
	rm -f $(ODIR)/*.o *.o *~ core
