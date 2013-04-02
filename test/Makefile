# You must use PrgEnv-cray.

CC = CC -O3
UPCC = cc -h upc -O -q
# If you change this, also change the mppwidth parameter in "job-knapsack" accordingly
NTHREADS = 4

TARGETS=serial knapsack

all: $(TARGETS)

serial: serial.cpp
	$(CC) -o $@ $<

knapsack: knapsack.upc
	$(UPCC) -o $@ $<

clean:
	rm -f *.o $(TARGETS)
