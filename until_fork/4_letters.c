#include "until_fork.h"

/* I plan to bring this functionality into the kernel space
   (write a module and add a system call) to (at least)
   avoid unnecessary context switches. */
   
/* we -> v -> vector / batch processing of threads */

_Bool spawn_and_wait(size_t T, void *d_arr, size_t d_elem_sz, void* (*f)(void *)) {
    
    _Bool R = 0; pthread_t th[T];
    
    fin(T) { PTHC(th + i, f, d_arr) { R = 1; T = i; break; } d_arr += d_elem_sz; }
    fin(T) PTHJ(th[i]) R = 1;
    
    ret R;
}

_Bool we_spawn(void) {
    
    ret 0;
}

_Bool we_wait(void) {
    
    ret 0;
}
