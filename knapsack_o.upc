#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>
#include <bupc_collectivev.h>

//
// auxiliary functions
//
inline int max( int a, int b ) { return a > b ? a : b; }
inline int min( int a, int b ) { return a < b ? a : b; }
double read_timer( )
{
    static int initialized = 0;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = 1;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  solvers
//
int build_table( int nitems, int cap, shared int *T, int *w, int *v )
{
    int wj, vj;

    wj = w[0];
    vj = v[0];
    upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) T[i] = 0;
    upc_forall( int i = wj; i <= cap; i++; &T[i] ) T[i] = vj;
    upc_barrier;

    for( int j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];
        upc_forall( int i = 0;  i <  wj;  i++; &T[i] ) T[i+cap+1] = T[i];
        upc_forall( int i = wj; i <= cap; i++; &T[i] ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        upc_barrier;

        T += cap+1;
    }

    return T[cap];
}

void backtrack( int nitems, int cap, shared int *T, int *w, shared [0] int *u )
{
    int i, j;

    if( MYTHREAD != 0 )
        return;

    i = nitems*(cap+1) - 1;
    for( j = nitems-1; j > 0; j-- )
    {
        u[j] = T[i] != T[i-cap-1];
        i -= cap+1 + (u[j] ? w[j] : 0 );
    }
    u[0] = T[i] != 0;
}

//
//  serial solver to check correctness
//
int solve_serial( int nitems, int cap, int *w, int *v )
{
    int i, j, best, *allocated, *T, wj, vj;

    //alloc local resources
    T = allocated = malloc( nitems*(cap+1)*sizeof(int) );
    if( !allocated )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    //build_table locally
    wj = w[0];
    vj = v[0];
    for( i = 0;  i <  wj;  i++ ) T[i] = 0;
    for( i = wj; i <= cap; i++ ) T[i] = vj;
    for( j = 1; j < nitems; j++ )
    {
        wj = w[j];
        vj = v[j];
        for( i = 0;  i <  wj;  i++ ) T[i+cap+1] = T[i];
        for( i = wj; i <= cap; i++ ) T[i+cap+1] = max( T[i], T[i-wj]+vj );
        T += cap+1;
    }
    best = T[cap];

    //free resources
    free( allocated );

    return best;
}

//
//  benchmarking program
//
int main( int argc, char** argv )
{
    int i, best_value, best_value_serial, total_weight, nused, total_value;
    double seconds;
    // shared int *weight;
    // shared int *value;
    shared [0] int *used;
    int *weight;
    int *value;
    // int *used;
    shared int *total;

    //these constants have little effect on runtime
    int max_value  = 1000;
    int max_weight = 1000;

    //these set the problem size
    int capacity   = 999;
    int nitems     = 5000;

    srand48( (unsigned int)time(NULL) + MYTHREAD );

    //allocate distributed arrays, use cyclic distribution
    // weight = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    // value  = (shared int *) upc_all_alloc( nitems, sizeof(int) );
    used   = (shared [0] int *) upc_all_alloc( nitems, sizeof(int) );
    weight = (int *) malloc( nitems * sizeof(int) );
    value = (int *) malloc( nitems * sizeof(int) );
    // used = (int *) malloc( nitems * sizeof(int) );
    total  = (shared int *) upc_all_alloc( nitems * (capacity+1), sizeof(int) );
    if( !weight || !value || !total || !used )
    {
        fprintf( stderr, "Failed to allocate memory" );
        upc_global_exit( -1 );
    }

    // init
    max_weight = min( max_weight, capacity );//don't generate items that don't fit into bag
    // upc_forall( i = 0; i < nitems; i++; i )
    // {
    //     weight[i] = 1 + (lrand48()%max_weight);
    //     value[i]  = 1 + (lrand48()%max_value);
    // }

    // upc_barrier;

    if( MYTHREAD == 0 ) {
        for( int i = 0; i < nitems; i++ ) {
            weight[i] = 1 + (lrand48()%max_weight);
            value[i]  = 1 + (lrand48()%max_value);
        }
    }

    for( int i = 0; i < nitems; i++ ) {
        weight[i] = bupc_allv_broadcast(int, weight[i], 0);
        value[i] = bupc_allv_broadcast(int, value[i], 0);
    }

    // time the solution
    seconds = read_timer( );

    best_value = build_table( nitems, capacity, total, weight, value );
    backtrack( nitems, capacity, total, weight, used );

    seconds = read_timer( ) - seconds;

    // check the result
    if( MYTHREAD == 0 )
    {
        printf( "%d items, capacity: %d, time: %g\n", nitems, capacity, seconds );

        best_value_serial = solve_serial( nitems, capacity, weight, value );

        total_weight = nused = total_value = 0;
        for( i = 0; i < nitems; i++ )
            if( used[i] )
            {
                nused++;
                total_weight += weight[i];
                total_value += value[i];
            }

        printf( "%d items used, value %d, weight %d\n", nused, total_value, total_weight );

        if( best_value != best_value_serial || best_value != total_value || total_weight > capacity )
            printf( "WRONG SOLUTION\n" );

        //release resources
        // upc_free( weight );
        // upc_free( value );
        upc_free( total );
        upc_free( used );
    }

    free( weight );
    free( value );
    // free( used );

    return 0;
}
