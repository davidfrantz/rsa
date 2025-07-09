# Libraries

GDAL=-I/usr/include/gdal -L/usr/lib -Wl,-rpath=/usr/lib

# Linked libs

LDGDAL=-lgdal

### COMPILER

GCC=gcc

CFLAGS=-O3 -Wall -fopenmp
#CFLAGS=-g -Wall -fopenmp


### TARGETS

all: max-ndvi rtm-inversion install clean
utils: alloc dir string stats table
.PHONY: all install clean



### UTILS COMPILE UNITS

alloc: utils/alloc.c
	$(GCC) $(CFLAGS) -c utils/alloc.c -o alloc.o

dir: utils/dir.c
	$(GCC) $(CFLAGS) -c utils/dir.c -o dir.o

string: utils/string.c
	$(GCC) $(CFLAGS) -c utils/string.c -o string.o

stats: utils/stats.c
	$(GCC) $(CFLAGS) -c utils/stats.c -o stats.o -lm

table: utils/table.c
	$(GCC) $(CFLAGS) -c utils/table.c -o table.o


### EXECUTABLES

max-ndvi: utils max-ndvi.c
	$(GCC) $(CFLAGS) $(GDAL) -o max-ndvi max-ndvi.c *.o $(LDGDAL) -lm

rtm-inversion: utils rtm-inversion.c
	$(GCC) $(CFLAGS) $(GDAL) -o rtm-inversion rtm-inversion.c *.o $(LDGDAL) -lm

  
### MISC

install:
	chmod 0755 max-ndvi
	chmod 0755 rtm-inversion
	cp max-ndvi /usr/local/bin/
	cp rtm-inversion /usr/local/bin/

clean:
	rm *.o

