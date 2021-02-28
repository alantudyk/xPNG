#include "until_fork.h"

static void* rsort_st(const size_t type, const _Bool asc,
                      void *const A, const size_t L, void *const E) {
                         
    size_t count[8][256] = {}, *m, r, offset[8][256]; if (L < 2) return A;
    
    s7_t  *sb;  u8_t      *b; s15_t *sw; u16_t *w, *we;
    s31_t *sd; u32_t *d, *de; s63_t *sq; u64_t *q, *qe;
    
    void *sa[] = { &&u1, &&s1, &&u2, &&s2, &&u4, &&s4, NULL, NULL, &&u8, &&s8 };
    void *ca[] = { NULL, NULL, &&c2, &&c3, &&c4, &&c5, NULL, NULL, &&c8, &&c9 };
    goto *sa[type];
    
s1: sb = A; m = count[0] + 128; fin(L) m[sb[i]]++;
    if (asc) fix(-128, 128, 1) { for (; m[i]--;) *sb++ = i; } else
           r_fix(-128, 128, 1) { for (; m[i]--;) *sb++ = i; } return A;
           
u1: b = A; fin(L) count[0][b[i]]++;
    if (asc) fin(256) { for (; count[0][i]--;) *b++ = i; } else
           r_fin(256) { for (; count[0][i]--;) *b++ = i; } return A;
           
s2: sw = A, w = A, we = E;
    fin(L) { w[i] = sw[i] + 0x8000; fiN(j, 2) count[j][((u8_t *)(w + i))[j]]++; } goto cx;
    
c3: fin(L) we[offset[0][((u8_t *)(w  + i))[0]]++] =  w[i];
    fin(L) sw[offset[1][((u8_t *)(we + i))[1]]++] = we[i] - 0x8000; return A;
    
u2: w = A, we = E; fin(L) { fiN(j, 2) count[j][((u8_t *)(w + i))[j]]++; } goto cx;

c2: fin(L) we[offset[0][((u8_t *)(w  + i))[0]]++] =  w[i];
    fin(L)  w[offset[1][((u8_t *)(we + i))[1]]++] = we[i]; return A;
    
s4: sd = A, d = A, de = E;
    fin(L) { d[i] = sd[i] + (1L << 31); fiN(j, 4) count[j][((u8_t *)(d + i))[j]]++; } goto cx;
    
c5: fiN(j, 3) { fin(L) de[offset[j][((u8_t *)(d + i))[j]]++] = d[i];
                    de = d, d = (d == A) ? E : A; }
    fin(L) sd[offset[3][((u8_t *)(d + i))[3]]++] = d[i] - (1L << 31); return A;
    
u4: d = A, de = E; fin(L) { fiN(j, 4) count[j][((u8_t *)(d + i))[j]]++; } goto cx;

c4: fiN(j, 4) { fin(L) de[offset[j][((u8_t *)(d + i))[j]]++] = d[i];
                    de = d, d = (d == A) ? E : A; } return A;

s8: sq = A, q = A, qe = E; s63_t x63 = (1UL << 63) - 1;
    fin(L) { q[i] = sq[i] < 0 ? (sq[i] + x63) + 1 : q[i] + (1UL << 63);
                fiN(j, 8) count[j][((u8_t *)(q +  i))[j]]++; } goto cx;

c9: fiN(j, 7) { fin(L) qe[offset[j][((u8_t *)(q + i))[j]]++] = q[i];
                    qe = q, q = (q == A) ? E : A; }
    fin(L) sq[offset[7][((u8_t *)(q + i))[7]]++] =
                     q[i] < (1UL << 63) ? (((s63_t)(q[i])) - x63) - 1
                                        : q[i] - (1UL << 63); return A;

u8: q = A, qe = E; fin(L) { fiN(j, 8) count[j][((u8_t *)(q + i))[j]]++; } goto cx;

c8: fiN(j, 8) { fin(L) qe[offset[j][((u8_t *)(q + i))[j]]++] = q[i];
                    qe = q, q = (q == A) ? E : A; } return A;

cx: fin(type & 14) {
    
        r = 0; if (asc) fiN(j, 256) offset[i][j] = r, r += count[i][j]; else
                      r_fiN(j, 256) offset[i][j] = r, r += count[i][j];
    
    }
    
    goto *ca[type];
}

typedef struct data_t { u64_t *count, *offset, L, type; void *A, *E; } data_t;

static PTHTF(rsort_T) {
    
    data_t *x = data;
    
    void *A = x->A, *E = x->E;
    u64_t type = x->type & 15, k = x->type >> 16, count[256], offset[256], *m, *n, L = x->L;
    _Bool s = (x->type >> 8) & 1;
    if (s) memset(count, 0, 2048); else { memcpy(offset, x->offset, 2048);
    if (type < 2) memcpy(count, x->count, 2048); }
    
    s7_t  *sb;  u8_t      *b; s15_t *sw; u16_t *w, *we; const s63_t x63 = (1UL << 63) - 1;
    s31_t *sd; u32_t *d, *de; s63_t *sq; u64_t *q, *qe;
    
    void *sa[] = { &&u1, &&s1, &&u2, &&s2, &&u4, &&s4, NULL, NULL, &&u8, &&s8 };
    void *ca[] = { &&c0, &&c1, &&c2, &&c3, &&c4, &&c5, NULL, NULL, &&c8, &&c9 };
    
    if (s) goto *sa[type]; else goto *ca[type];
    
s1: sb = A, m = count + 128; fin(L) m[sb[i]]++; goto r;
c1: sb = A, n = offset + 128, m = count + 128;
    fix(-128, 128, 1) { for (; m[i]--;) sb[n[i]++] = i; } goto r;
    
u1: b  = A; fin(L) count[b[i]]++; goto r;

c0: b  = A; fin(256) { for (; count[i]--;) b[offset[i]++] = i; } goto r;

s2: sw = A, w = A;
    fin(L) { if (!k) w[i] = sw[i] + 0x8000; count[((u8_t *)(w + i))[k]]++; } goto r;
    
c3: w = A; if (!k) { we = E; fin(L) we[offset[((u8_t *)(w  + i))[0]]++] = w[i]; }
    else { sw = E; fin(L) sw[offset[((u8_t *)(w + i))[1]]++] = w[i] - 0x8000; } goto r;
    
u2: w = A; fin(L) count[((u8_t *)(w + i))[k]]++; goto r;

c2: w = A, we = E; fin(L) we[offset[((u8_t *)(w  + i))[k]]++] = w[i]; goto r;

s4: sd = A, d = A;
    fin(L) { if (!k) d[i] = sd[i] + (1L << 31); count[((u8_t *)(d + i))[k]]++; } goto r;
    
c5: d = A; if (k < 3) { de = E; fin(L) de[offset[((u8_t *)(d + i))[k]]++] = d[i]; }
    else { sd = E; fin(L) sd[offset[((u8_t *)(d + i))[3]]++] = d[i] - (1L << 31); } goto r;
    
u4: d = A; fin(L) count[((u8_t *)(d + i))[k]]++; goto r;

c4: d = A, de = E; fin(L) de[offset[((u8_t *)(d + i))[k]]++] = d[i]; goto r;

s8: sq = A, q = A;
    fin(L) { if (!k) q[i] = sq[i] < 0 ? (sq[i] + x63) + 1 : q[i] + (1UL << 63);
                count[((u8_t *)(q + i))[k]]++; } goto r;
                
c9: q = A; if (k < 7) { qe = E; fin(L) qe[offset[((u8_t *)(q + i))[k]]++] = q[i]; }
    else { sq = E; fin(L) sq[offset[((u8_t *)(q + i))[7]]++] =
            q[i] < (1UL << 63) ? (((s63_t)(q[i])) - x63) - 1 : q[i] - (1UL << 63); } goto r;

u8: q = A; fin(L) count[((u8_t *)(q + i))[k]]++; goto r;

c8: q = A, qe = E; fin(L) qe[offset[((u8_t *)(q + i))[k]]++] = q[i]; goto r;

r:  if (s) memcpy(x->count, count, 2048); return NULL;
}

void* rsort(size_t T, const char *const Type, const char *const Ord,
            void *const A, const size_t L, void *const Ext_Buff) {
    
    _Bool asc; if (!strcmp(Ord,  "asc") || !strcmp(Ord, "<")) asc = 1; else
               if (!strcmp(Ord, "desc") || !strcmp(Ord, ">")) asc = 0; else ret NULL;
    
    const char *const ta[] = {
        "u8_t", "s7_t", "u16_t", "s15_t", "u32_t", "s31_t", NULL, NULL, "u64_t", "s63_t"
    };
    
    size_t _type = -1;
    for (; ++_type < 10; ) {
        if (ta[_type] == NULL) continue;
        if (!strcmp(Type,  ta[_type])) break;
    }
    if (_type == 10) ret NULL;
    const size_t type = _type, type_size = (_type &= 14) ? _type : 1;
    void *const E = (Ext_Buff == NULL && type_size > 1) ? malloc(L * type_size) : Ext_Buff;
    if (E == NULL && type_size > 1) ret NULL;
    
    if (T == 0) T = sysconf(_SC_NPROCESSORS_ONLN);
    for (; L / T < (u64_t)1e5 && T > 1; T--);
    if (T == 1) {
        void *r = rsort_st(type, asc, A, L, E);
        if (Ext_Buff == NULL) {
            if (r != A) { r = A; memcpy(A, E, type_size * L); }
            free(E);
        }
        ret r;
    }
    
    data_t d[T]; size_t (* count)[256], (* offset)[256], r, l = L / T;
    
    CALLOC(count, sizeof(size_t) * 256 * T * 2) ret NULL; offset = count + T;
    void *a = A, *e = E;
    
    fin(T) d[i].count = count[i], d[i].offset = offset[i], d[i].L = (i < T - 1) ? l : l + L % T;
    
    fiN(k, type_size) {
        
        fin(T) d[i].type = (k << 16) + 256 + type, d[i].A = a + i * l * type_size;
        if (spawn_and_wait(T, d, sizeof(data_t), rsort_T)) goto r;
        
        r = 0; if (asc) fiN(j, 256) fin(T) offset[i][j] = r, r += count[i][j]; else
                      r_fiN(j, 256) fin(T) offset[i][j] = r, r += count[i][j];
        
        fin(T)
            d[i].type = (k << 16) + type,
            d[i].A    = (type_size > 1) ? a + i * l * type_size : A,
            d[i].E    = e;
        if (spawn_and_wait(T, d, sizeof(data_t), rsort_T)) goto r;
        
        a = e, e = (e == E) ? A : E;
    }
    
r:  free(count); return a;

}
