
/*

    Big_O:

    1e8 elems:
        - MSORT:           27.00 N
        - MSORT_MT:
            -   2 threads: 14.00 N
            -   4 threads:  7.75 N
            -   8 threads:  4.75 N
            -  16 threads:  3.30 N
            -  32 threads:  2.60 N
            -  64 threads:  2.30 N
            - 128 threads:  2.15 N
            - 256 threads:  2.05 N

    def big_O log_N
        log_N = Float log_N
        for i in 1..8
            return if i > log_N
            x = (log_N - i) / 2 ** i
            a = 1.0; i.times { x += a; a /= 2 }
            puts "\t- %3d: %13.10f" % [2 ** i, x]
        end
    end

*/

#ifndef ___MSORT_H
#define ___MSORT_H

#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MCMP(func_name) _Bool func_name(const void *a, const void *b)

typedef _Bool (* ___cmp_t) (const void *, const void *);

#define MSORT_DEC(func_name) \
void* func_name(void *Arr, size_t Arr_Len, void *Ext_Buff, ___cmp_t Is_Unordered);

#define MSORT_MT_DEC(func_name) \
void* func_name(size_t T, void *Arr, size_t Arr_Len, void *Ext_Buff);

#define ___I *x++ = *i++
#define ___J *x++ = *j++

#define MSORT(func_name, TYPE) /* TYPE: size_t for pointers */ \
void* func_name(void *const A, const size_t L, void *const E, ___cmp_t cmp) { \
    if (L < 2) return A; size_t z = 2, y = 1, l = L / 2, t = 0; \
    TYPE tmp, *i = A, *j = i  + 1, *I, *J, *a = A, *e = E, *x = a + L, *p; \
    for (; l--; i += 2, j += 2) { if (cmp(i, j)) { tmp = *j; *j = *i; *i = tmp; } } goto t; \
y2: for (; l--; I += 4, J += 4, i += 2, j += 2) { \
        if (cmp(i, j)) ___J; else ___I; \
        if (cmp(i, j)) { ___J; if (j == J) { ___I; ___I; continue; } } \
        else { ___I; if (i == I) { ___J; ___J; continue; } } \
        if (cmp(i, j)) { ___J; ___I; } else { ___I; ___J; } } goto t; \
y4: for (; l--; I += 8, J += 8, i += 4, j += 4) { \
        for (size_t k = 4; --k;) { if (cmp(i, j)) ___J; else ___I; } \
        for (;;) { if (cmp(i, j)) { ___J; if (j == J) { for (; i < I;) ___I; break; } } \
                   else { ___I; if (i == I) { for (; j < J;) ___J; break; } } } } goto t; \
y_:  for (; l--; I += z, J += z, i += y, j += y) { \
        for (size_t a, b, k = y;;) { for (; --k;) { if (cmp(i, j)) ___J; else ___I; } \
                                     a = I - i, b = J - j, k = a > b ? b : a; if (k < 4) break; } \
        for (;;) { if (cmp(i, j)) { ___J; if (j == J) { for (; i < I;) ___I; break; } } \
                   else { ___I; if (i == I) { for (; j < J;) ___J; break; } } } } \
t:  if (t == 0) { if (L % z == 0) goto n; else { \
        if (z == 2) { i = (I = j = (J = a + L) - 1) - 2; x = e + (L - 3); t = 3; } \
        else        { i = (I = x) - z; x = (j = (J = a + L) - y) - z; t = z + y; } } } \
    else { J = (j = p) + t; if (L % z < t) { t += y; } \
                            else { i = (I = x) - z; x = a + (i - e); t += z; } } p = x; \
    for (size_t a, b, k = y;;) { for (; --k;) { if (cmp(i, j)) ___J; else ___I; } \
                                 a = I - i, b = J - j, k = a > b ? b : a; if (k < 4) break; } \
    for (;;) { if (cmp(i, j)) { ___J; if (j == J) { for (; i < I;) ___I; break; } } \
               else { ___I; \
                   if (i == I) { if (x != j) { for (; j < J;) ___J; } else x = J; break; } } } \
n:  if (L / z == 1) return (x == a + L) ? a : e; if (z > 2) a = e, e = (e == E) ? A : E; \
    z = 2 * (y *= 2), l = L / z, J = (j = (I = (i = a) + y)) + y, x = e; if (L % z < t) l--; \
    if (y > 4) goto y_; else if (y == 2) goto y2; else goto y4; }

typedef struct ___mttfd_t { void *A, *E, *V; size_t L, X; } ___mttfd_t;

#define MSORT_MT(before, func_name, TYPE, msort_func_name, cmp) /* TYPE: size_t for pointers */ \
static void* func_name##_T1(void *data) { ___mttfd_t *d = data; \
    d->A = msort_func_name(d->A, d->L, d->E, cmp); return NULL; } \
static void* func_name##_T2(void *data) { ___mttfd_t *d = data; \
    TYPE *i = d->A, *I = i + d->L, \
         *j = (d->V == NULL) ? I : d->V, *J = j + (d->L + d->X), *x = d->E; \
    for (size_t a, b, k = d->L;;) { for (; --k;) { if (cmp(i, j)) ___J; else ___I; } \
                                    a = I - i, b = J - j, k = a > b ? b : a; if (k < 4) break; } \
    for (;;) { if (cmp(i, j)) { ___J; if (j == J) { for (; i < I;) ___I; break; } } \
               else { ___I; if (i == I) { if (x != j) { for (; j < J;) ___J; } break; } } } \
    return NULL; } \
before void* func_name(size_t T, void *const A, const size_t L, void *const E) { \
    if (T < 2 || T > 256 || (T & (T - 1))) return NULL; \
    if (L < 4096) return msort_func_name(A, L, E, cmp); \
    for (; L / T < 2048; T /= 2); ___mttfd_t d[T]; _Bool R = 0; pthread_t th[T]; \
    size_t t, P = T / 2, y = L / T, z = y * 2, X = L % T; TYPE *a = A, *e = E, *i = a, *j = e, *V; \
    for (t = 0; t < T; t++) { \
        d[t].A = i; i += y; d[t].E = j; j += y; d[t].L = (t < T - 1) ? y : y + X; \
        if (pthread_create(th + t, NULL, func_name##_T1, d + t) != 0) { R = 1; T = t; break; } } \
    for (t = 0; t < T; t++) if (pthread_join(th[t], NULL) != 0) R = 1; if (R) return NULL; \
    if (d[0].A == d[0].E) { a = E, e = A; V = (d[T - 1].A != d[T - 1].E) ? d[T - 1].A : NULL; } \
    else V = (d[T - 1].A == d[T - 1].E) ? d[T - 1].A : NULL; \
    for(i = a, j = e, T /= 2; T; T /= 2, y *= 2, z *= 2, a = e, e = (e == E) ? A : E, i = a, j = e) { \
        for (t = 0; t < T; t++) { \
            d[t].A = i; i += z; d[t].E = j; j += z; d[t].L = y; d[t].X = (t < T - 1) ? 0 : X; \
            d[t].V = (T == P && t == T - 1) ? V : NULL; \
            if (pthread_create(th + t, NULL, func_name##_T2, d + t) != 0) { R = 1; T = t; break; } } \
        for (t = 0; t < T; t++) if (pthread_join(th[t], NULL) != 0) R = 1; if (R) return NULL; } \
    return a; }

#endif
