CFLAGS=		-g -Wall -O2 -ftree-vectorize
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE

profiler: profiler.o distance.o simpleRepeat.o
	gcc $(CFLAGS) $(DFLAGS) profiler.o distance.o simpleRepeat.o \
	libhts.a -lz -lstdc++ -pthread -o profiler

profiler.o: profiler.cpp
	gcc $(CFLAGS) $(DFLAGS) -c profiler.cpp

distance.o: distance.cpp
	gcc $(CFLAGS) $(DFLAGS) -c distance.cpp

simpleRepeat.o: simpleRepeat.cpp
	gcc $(CFLAGS) $(DFLAGS) -c simpleRepeat.cpp

clean:
	rm -f profiler
	rm -f profiler.o
	rm -f distance.o
