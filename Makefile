CC=gcc
OUTPUT_OPTION=-Wall -O2 -ferror-limit=50 -MMD -MP -o $@
LDLIBS=-lm
 
SRC=$(wildcard src/*.c)
OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)
 
.PHONY: clean
 
sumatra: $(OBJ)
	$(LINK.o) -O2 $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f $(OBJ) $(DEP) sumatra

_temp:
	gcc -c src/smtr.c
	ar rcs libsmtr.a smtr.o
	mv libsmtr.a ~/Dropbox/tmp/brewer/macOS
	cp src/smtr.h ~/Dropbox/tmp/brewer/macOS
	cp src/vec3.h ~/Dropbox/tmp/brewer/macOS

-include $(DEP)
