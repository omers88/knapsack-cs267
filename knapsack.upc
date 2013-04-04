#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>
#include <upc_collective.h>

#define CAPACITY (1536 - 1)
#define NITEMS 6144
#define CYCLES 1
#define BLOCK_WIDTH ((CAPACITY + 1) / THREADS)
#define BLOCK_HEIGHT (NITEMS / THREADS / CYCLES)
#define PROCESSOR_BLOCK_SIZE (BLOCK_HEIGHT * (CAPACITY + 1))
#define NUM_LOCKS (THREADS * THREADS * CYCLES)

//
// auxiliary functions
//
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
upc_lock_t* shared *locks_shared;
double read_timer( ) {
    static int initialized = 0;
    static struct timeval start;
    struct timeval end;
    if( !initialized ) {
        gettimeofday( &start, NULL );
        initialized = 1;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  solvers
//
int build_table( shared [PROCESSOR_BLOCK_SIZE] int *T, shared [BLOCK_HEIGHT] int *w, shared [BLOCK_HEIGHT] int *v, upc_lock_t** locks ) {
    int w_item, v_item, col_index, row_start, prev_row_start, lock_index, T_start, T_start_prev;
    int *T_local, *w_local, *v_local;

    // Calculate top row
    if (MYTHREAD == 0) {
        for (int i = 0; i < w[0]; i++) T[i] = 0;
        for (int i = w[0]; i <= CAPACITY; i++) T[i] = v[0];
    }
    
    for (int cycle = 0; cycle < CYCLES; cycle++) {
        T_start = (CAPACITY+1) * (BLOCK_HEIGHT * (MYTHREAD + cycle*THREADS));
        T_start_prev = T_start - (CAPACITY+1);
        T_local = (int *) (T + T_start);
        w_local = (int *) (w + BLOCK_HEIGHT * (MYTHREAD + cycle*THREADS));
        v_local = (int *) (v + BLOCK_HEIGHT * (MYTHREAD + cycle*THREADS));

        for (int block = 0; block < THREADS; block++) {
            col_index = BLOCK_WIDTH * block;

            // Acquire lock
            lock_index = THREADS * (MYTHREAD + cycle*THREADS) + block;
            upc_lock(locks[lock_index]);

            // Calculate top row of block
            if (!(cycle == 0 && MYTHREAD == 0)) {
                w_item = w_local[0];
                v_item = v_local[0];
                for (int col = col_index; col < col_index + BLOCK_WIDTH; col++) {
                    if (col < w_item) {
                        T_local[col] = T[T_start_prev + col];
                    } else {
                        T_local[col] = max( T[T_start_prev + col], T[T_start_prev + col - w_item] + v_item );
                    }
                }
            }
            // Do rest of block with private pointer
            for (int row = 1; row < BLOCK_HEIGHT; row++) {
                w_item = w_local[row];
                v_item = v_local[row];
                row_start = row * (CAPACITY + 1);
                prev_row_start = row_start - (CAPACITY+1);
                for (int col = col_index; col < col_index + BLOCK_WIDTH; col++) {
                    if (col < w_item) {
                        T_local[row_start + col] = T_local[prev_row_start + col];
                    } else {
                        T_local[row_start + col] = max( T_local[prev_row_start + col], T_local[prev_row_start + col - w_item] + v_item );
                    }
                }
            }

            // Free next lock
            upc_unlock(locks[lock_index]);
            if (!(MYTHREAD == THREADS - 1 && cycle == CYCLES - 1)) {
                upc_unlock(locks[lock_index + THREADS]);
            }
        }
    }

    upc_barrier;
    return T[(CAPACITY+1) * NITEMS - 1];
}

void backtrack( shared [PROCESSOR_BLOCK_SIZE] int *T, shared [BLOCK_HEIGHT] int *w, shared int *u ) {
    int i, j;

    if( MYTHREAD != 0 ) return;

    i = NITEMS*(CAPACITY+1) - 1;
    for( j = NITEMS-1; j > 0; j-- ) {
        u[j] = T[i] != T[i-CAPACITY-1];
        i -= CAPACITY+1 + (u[j] ? w[j] : 0 );
    }
    u[0] = T[i] != 0;
}

//
//  serial solver to check correctness
//
int solve_serial( shared [BLOCK_HEIGHT] int *w, shared [BLOCK_HEIGHT] int *v ) {
    int i, j, best, *allocated, *T, wj, vj;

    // alloc local resources
    T = allocated = malloc( NITEMS*(CAPACITY+1)*sizeof(int) );
    if( !allocated ) {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    // build_table locally
    wj = w[0];
    vj = v[0];
    for( i = 0;  i <  wj;  i++ ) T[i] = 0;
    for( i = wj; i <= CAPACITY; i++ ) T[i] = vj;
    for( j = 1; j < NITEMS; j++ ) {
        wj = w[j];
        vj = v[j];
        for( i = 0;  i <  wj;  i++ ) T[i+CAPACITY+1] = T[i];
        for( i = wj; i <= CAPACITY; i++ ) T[i+CAPACITY+1] = max( T[i], T[i-wj]+vj );
        T += CAPACITY+1;
    }
    best = T[CAPACITY];

    // free resources
    free( allocated );

    return best;
}

upc_lock_t **allocate_lock_array(unsigned int count) {
    const unsigned int blksize = ((count + THREADS - 1) / THREADS); // Roundup
    upc_lock_t* shared *tmp = upc_all_alloc(blksize*THREADS, sizeof(upc_lock_t*));
    locks_shared = upc_all_alloc(blksize*THREADS*THREADS, sizeof(upc_lock_t*));

    // Allocate lock pointers into a temporary array.
    // This code overlays an array of blocksize [*] over the cyclic one.
    upc_lock_t** ptmp = (upc_lock_t**)(&tmp[MYTHREAD]); // Local array "slice"
    const int my_count = upc_affinitysize(count,blksize,MYTHREAD);

    for (int i=0; i<my_count; ++i) {
        ptmp[i] = upc_global_lock_alloc();
    }

    // Replicate the temporary array THREADS times into the locks_shared array
    // IN_MYSYNC:   Since each thread generates its local portion of input.
    // OUT_ALLSYNC: Ensures upc_free() occurs only after tmp is unneeded.
    upc_all_gather_all(locks_shared, tmp, blksize * sizeof(upc_lock_t*), UPC_IN_MYSYNC|UPC_OUT_ALLSYNC);

    if (MYTHREAD == 0) {
        upc_free(tmp);  // Free the temporary array
    }
    // Return a pointer-to-private for local piece of replicated locks_shared
    return (upc_lock_t**)(&locks_shared[MYTHREAD]);
}

//
//  benchmarking program
//
int main( int argc, char** argv ) {
    int i, best_value, best_value_serial, total_weight, nused, total_value;
    double seconds;
    shared [BLOCK_HEIGHT] int *weight;
    shared [BLOCK_HEIGHT] int *value;
    shared int *used;
    shared [PROCESSOR_BLOCK_SIZE] int *total;
    upc_lock_t **locks;

    // these constants have little effect on runtime
    int max_value  = 1000;
    int max_weight = 1000;

    srand48( (unsigned int)time(NULL) + MYTHREAD );

    // allocate distributed arrays, use cyclic distribution
    weight = (shared [BLOCK_HEIGHT] int *) upc_all_alloc( THREADS * CYCLES, BLOCK_HEIGHT * sizeof(int) );
    value  = (shared [BLOCK_HEIGHT] int *) upc_all_alloc( THREADS * CYCLES, BLOCK_HEIGHT * sizeof(int) );
    used   = (shared int *) upc_all_alloc( NITEMS, sizeof(int) );
    total  = (shared [PROCESSOR_BLOCK_SIZE] int *) upc_all_alloc( THREADS * CYCLES, PROCESSOR_BLOCK_SIZE * sizeof(int) );

    if( !weight || !value || !total || !used ) {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    // init
    max_weight = min( max_weight, CAPACITY ); // don't generate items that don't fit into bag
    upc_forall( i = 0; i < NITEMS; i++; i ) {
        weight[i] = 1 + (lrand48()%max_weight);
        value[i]  = 1 + (lrand48()%max_value);
    }
    upc_barrier;

    // locks = (shared upc_lock_t**) upc_all_alloc( NUM_LOCKS, sizeof(upc_lock_t*) );
    locks = allocate_lock_array(NUM_LOCKS);
    for (int cycle = 0; cycle < CYCLES; cycle++) {
        upc_forall( int row = 0; row < THREADS; row++; row-1+THREADS ) {
            if (!(cycle == 0 && row == 0)) {
                for ( int col = 0; col < THREADS; col++ ) {
                    upc_lock(locks[THREADS * (row + cycle*THREADS) + col]);
                }
            }
        }
    }
    upc_barrier;

    // time the solution
    seconds = read_timer( );

    best_value = build_table( total, weight, value, locks );
    backtrack( total, weight, used );

    seconds = read_timer( ) - seconds;

    // check the result
    if( MYTHREAD == 0 ) {
        printf( "%d items, capacity: %d, time: %g\n", NITEMS, CAPACITY, seconds );

        best_value_serial = solve_serial( weight, value );

        total_weight = nused = total_value = 0;
        for( i = 0; i < NITEMS; i++ ) {
            if( used[i] ) {
                nused++;
                total_weight += weight[i];
                total_value += value[i];
            }
        }

        printf( "%d items used, value %d, weight %d\n", nused, total_value, total_weight );

        if( best_value != best_value_serial || best_value != total_value || total_weight > CAPACITY ) {
            printf( "WRONG SOLUTION\n" );
        }
    }

    // release resources
    if( MYTHREAD == 0 ) {
        for( int i = 0; i < NUM_LOCKS; i++ ) {
            upc_lock_free(locks[i]);
        }
        upc_free( locks_shared );
        upc_free( weight );
        upc_free( value );
        upc_free( total );
        upc_free( used );
    }
    return 0;
}
