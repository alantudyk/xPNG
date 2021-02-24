#ifndef _UNTIL_FORK_H
#define _UNTIL_FORK_H

#include <pthread.h>
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

#define ret return

#define fin(stop) for (s63_t i = 0, S_T_O_P_ = (stop); i < S_T_O_P_; i++)
#define fix(start, stop, increment) \
    for (s63_t i = (start), S_T_O_P_ = (stop), I_N_C_ = (increment); i < S_T_O_P_; i += I_N_C_)
#define fiN(iter_name, stop) \
    for (s63_t iter_name = 0, S_T_O_P_ = (stop); iter_name < S_T_O_P_; iter_name++)
#define fiX(iter_name, start, stop, increment) \
    for (s63_t iter_name = (start), S_T_O_P_ = (stop), I_N_C_ = (increment); \
    iter_name < S_T_O_P_; iter_name += I_N_C_)

#define PTHSPI(pointer_to_spinlock) pthread_spin_init(pointer_to_spinlock, 0);
#define PTHSPL(pointer_to_spinlock) pthread_spin_lock(pointer_to_spinlock);
#define PTHSPU(pointer_to_spinlock) pthread_spin_unlock(pointer_to_spinlock);
#define PTHSPD(pointer_to_spinlock) pthread_spin_destroy(pointer_to_spinlock);
#define PTHSPT pthread_spinlock_t
#define PTHTF(func_name) void* func_name(void *data)
#define PTHC(id_ptr, thread_func, data_ptr) \
if (pthread_create(id_ptr, NULL, thread_func, data_ptr) != 0)
#define PTHJ(id) if (pthread_join(id, NULL) != 0)

_Bool spawn_and_wait(u64_t T, void *d_arr, u64_t d_elem_sz, void* (*f)(void *));

#endif
