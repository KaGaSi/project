CC = gcc
CFLAGS = -std=c11 -Wall -Werror -O3
LIBS = -lm

SOURCES = TransformVsf.c Aux.c Structure.c

DEPS = Aux.h Structure.h CStructs.h

ODIR = obj
_OBJ = $(SOURCES:.c=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

BIN = TransformVsf #Trans

all: $(BIN) $(OBJ)

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

TransformVsf: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#Trans: $(OBJ)
#	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#copy:
#	cp $(BIN) bin

clean:
	rm -f $(ODIR)/*.o *.o *~ core
