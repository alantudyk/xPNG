#include "../until_fork.h"

static inline int dcmp(const void *a, const void *b) {
    return (*(const s31_t *)a > *(const s31_t *)b) -
           (*(const s31_t *)a < *(const s31_t *)b);
}

#define N (size_t)1e8
static s31_t R[N], A[N], Q[N];

#define ERR_MSG "Error, line %d\n\n"

MAIN_VOID {
    
    const size_t T = sysconf(_SC_NPROCESSORS_ONLN);
    
    TIME_PAIR; size_t ms;
    srand(time(NULL)); fin(N) R[i] = rand();
    
    pf("%15d thread  ", 1); memcpy(A, R, N * 4);
    TIME_DIFF_EXEC(
        if (rsort(1, "s31_t", "<", A, N, NULL) ==  NULL) { pf(ERR_MSG, 22); return 1; },
        MS, ms);
    pf("%13lu ms\n", ms);
    
    memcpy(Q, R, N * 4); qsort(Q, N, 4, dcmp);
    if (memcmp(A, Q, N * 4)) { pf(ERR_MSG, 27); return 1; }
    
    if (T < 2) goto r;
    
    pf("%15zu threads ", T); memcpy(A, R, N * 4);
    TIME_DIFF_EXEC(
        if (rsort(T, "s31_t", "<", A, N, NULL) ==  NULL) { pf(ERR_MSG, 33); return 1; },
        MS, ms);
    pf("%13lu ms\n", ms);
    
    memcpy(Q, R, N * 4); qsort(Q, N, 4, dcmp);
    if (memcmp(A, Q, N * 4)) { pf(ERR_MSG, 38); return 1; }
    
r:  pf("\n");

    ret 0;
}
