all:
	gcc src/optimpoly.c -lnlopt -lm -o test

clean:
	rm -rf *.o test
