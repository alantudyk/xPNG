#ifndef ___UNTIL_FORK_H
#define ___UNTIL_FORK_H

#define static_assert(CND) _Static_assert(CND, "");

static_assert(__SIZEOF_POINTER__ == 8);

#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>

#define MAIN_VOID int main(void)
#define MAIN_ARGS int main(const int argc, const char *const *const argv)

#define TIME_PAIR      struct timespec ___t1, ___t2
#define TIME_GET_START clock_gettime(CLOCK_REALTIME, &___t1)
#define TIME_GET_STOP  clock_gettime(CLOCK_REALTIME, &___t2)
#define TIME_DIFF_NS   ((___t2.tv_nsec + ___t2.tv_sec * (int64_t)1e9) - \
                        (___t1.tv_nsec + ___t1.tv_sec * (int64_t)1e9))
#define TIME_DIFF_MS   (TIME_DIFF_NS / (int64_t)1e6)
#define TIME_DIFF_MCS  (TIME_DIFF_NS / (int64_t)1e3)
#define TIME_DIFF_EXEC(EXEC, DIFF_SUFFIX, DIFF_LVAL) \
    TIME_GET_START; EXEC; TIME_GET_STOP; DIFF_LVAL = TIME_DIFF_##DIFF_SUFFIX

typedef  int8_t  s7_t; typedef  uint8_t  u8_t;
typedef int16_t s15_t; typedef uint16_t u16_t;
typedef int32_t s31_t; typedef uint32_t u32_t;
typedef int64_t s63_t; typedef uint64_t u64_t;

typedef struct reg_t { u8_t *p; u64_t s; } reg_t;

#define CHECK    __attribute__((warn_unused_result))
#define NOINLINE __attribute__((noinline))

_Bool  f_size(const char *fn,       u64_t *s);
_Bool  f_read(const char *fn,       reg_t *f);
_Bool f_write(const char *fn, const reg_t *f);

#define  F_SIZE(filename,   s_ptr) if ( f_size(filename,   s_ptr))
#define  F_READ(filename, reg_ptr) if ( f_read(filename, reg_ptr))
#define F_WRITE(filename, reg_ptr) if (f_write(filename, reg_ptr))

#define MALLOC(pointer, reg_size) if ((pointer = malloc(reg_size))    == NULL)
#define CALLOC(pointer, reg_size) if ((pointer = calloc(reg_size, 1)) == NULL)

#define pf printf
#define ret return
#define nil NULL
#define not(CND) !(CND)
#define unless(CND) if (not(CND))
#define BETWEEN(l, v, r) ((l) <= (v) && (v) <= (r))
#define BITMASK_SHL(n, shl) (((1LU << (n)) - 1) << (shl))
#define BITMASK(n) BITMASK_SHL(n, 0)

#define fin(stop) for (ssize_t i = -1, _stop = stop; ++i < _stop;)
#define fiN(iter_name, stop) for (ssize_t iter_name = -1, _stop = stop; ++iter_name < _stop;)
#define fix(start, stop, increment) \
    for (ssize_t i = start, _stop = stop; i < _stop; i += increment)
#define fiX(iter_name, start, stop, increment) \
    for (ssize_t iter_name = start, _stop = stop; iter_name < _stop; iter_name += increment)
#define r_fin(stop) for (ssize_t i = stop; --i >= 0;)
#define r_fiN(iter_name, stop) for (ssize_t iter_name = stop; --iter_name >= 0;)
#define r_fix(start, stop, decrement) \
    for (ssize_t i = stop, _start = start; (i -= decrement) >= _start;)
#define r_fiX(iter_name, start, stop, decrement) \
    for (ssize_t iter_name = stop, _start = start; (iter_name -= decrement) >= _start;)

#define PTHSPI(pointer_to_spinlock) pthread_spin_init(pointer_to_spinlock, 0);
#define PTHSPL(pointer_to_spinlock) pthread_spin_lock(pointer_to_spinlock);
#define PTHSPU(pointer_to_spinlock) pthread_spin_unlock(pointer_to_spinlock);
#define PTHSPD(pointer_to_spinlock) pthread_spin_destroy(pointer_to_spinlock);
#define PTHSPT pthread_spinlock_t
#define PTHTF(func_name) void* func_name(void *data)
#define PTHC(id_ptr, thread_func, data_ptr) \
    if (pthread_create(id_ptr, NULL, thread_func, data_ptr) != 0)
#define PTHJ(id) if (pthread_join(id, NULL) != 0)

_Bool spawn_and_wait(size_t T, void *d_arr, size_t d_elem_sz, void* (*f)(void *));

void* rsort(size_t T, const char *Type, const char *Ord, void *A, size_t L, void *Ext_Buff);

#endif
