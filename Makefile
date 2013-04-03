# You must use PrgEnv-cray.

CC = CC -O3
UPCC = cc -h upc -O
#UPCC = upcc -O
# If you change this, also change the mppwidth parameter in "job-knapsack" accordingly
NTHREADS = 4

TARGETS=serial knapsack knapsack_orig

all: $(TARGETS)

serial: serial.cpp
	$(CC) -o $@ $<

knapsack: knapsack.upc
	$(UPCC) -o $@ $<

knapsack_orig: knapsack_orig.upc
	$(UPCC) -o $@ $<

clean:
	rm -f *.o $(TARGETS)