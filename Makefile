all: plschain

plschain: main.c index.h kseq.h khash.h kmer.h query.h
	gcc -Wall -g -O3 -fsanitize=address main.c -o plschain -lz

clean:
	rm -f *.o plschain
