#include <until_fork.h>

#if T_MAX < 1 || T_MAX > 256
#error T_MAX out of 1..256
#endif

typedef struct reg_t { u8_t *p; u64_t s; }      reg_t;
typedef struct  pm_t { u8_t *p; u64_t w, h, s; } pm_t;

#define F_SZ(filename,    s) if (f_size(filename,    s) != 0)
#define F_RD(filename, file) if (f_read(filename, file) != 0)
#define A (s63_t)3
#define BC(C, B) b->bc = (b->bc << (C)) + (B), b->l += (C)
#define BCF(C) ((b->bc >> (b->l -= (C))) & bitmask(C))
#define EqN == NULL

#define bitmask(n) ((1LU << (n)) - 1)
#define numBit(v) ((1 + (__builtin_clz(v) ^ 31)) & (-(v) >> 31))
#define pix_toU fin(3) { pix[i] = (s7_t)pix[i]; pix[i] = (pix[i] << 1) ^ (pix[i] >> 31); }
#define pix_toS fin(3) pix[i] = (pix[i] >> 1) ^ (-(pix[i] & 1));
#define FILL  if (b->l  < 32) { b->bc <<= 32, b->l += 32; if (b->_p < b->DP) b->bc += *(b->_p)++; }
#define FLUSH if (b->l >= 32) *(b->_p)++ = b->bc >> (b->l -= 32)
#define pix_copy(to, from) fin(3) to[i] = from[i];
#define p1x(operator) fin(3) pix[i] = pix[i] operator p[i - A];
#define p1y(operator) fin(3) pix[i] = pix[i] operator p[i - bpr];
#define p2a(operator) fin(3) pix[i] = pix[i] operator ((p[i - A] + p[i - bpr] + 1) >> 1);
#define p3a(operator) fin(3) pix[i] = \
    pix[i] operator ((((3 * p[i - A] + 3 * p[i - bpr]) - 2 * p[(i - bpr) - A]) + 2) >> 2)
#define numB(from) nl = from[0] | from[1] | from[2]; nl = numBit(nl);
#define ENC(pr) \
    pix_copy(pix, p); pr(-); pix_toU; numB(pix); \
    F[0][(pl << 4) + nl]++, pl = *cx[pl]++ = nl; \
    switch (nl) { \
        case  0: break; \
        case  1: F[1][*st[1]++ = (pix[0] << 2) | (pix[1] << 1) | pix[2]]++; break; \
        case  2: F[2][*st[2]++ = (pix[0] << 4) | (pix[1] << 2) | pix[2]]++; break; \
        default: fin(3) F[nl][*st[nl]++ = pix[i]]++; \
    } p += A;
#define ENC4 \
    pix_copy(pix, p); if (Y) p3a(-); else p2a(-); \
    if (G) pix[0] -= pix[1], pix[2] -= pix[1]; pix_toU; numB(pix); \
    F[0][(pl << 4) + nl]++, pl = *cx[pl]++ = nl; \
    switch (nl) { \
        case  0: break; \
        case  1: F[1][*st[1]++ = (pix[0] << 2) | (pix[1] << 1) | pix[2]]++; break; \
        case  2: F[2][*st[2]++ = (pix[0] << 4) | (pix[1] << 2) | pix[2]]++; break; \
        default: fin(3) F[nl][*st[nl]++ = pix[i]]++; \
    } p += A;


static u64_t f_size(const char *fn, u64_t *s) {
    struct stat st; if (stat(fn, &st) != 0) ret 1; *s = st.st_size; ret 0;
}

#define RET fclose(fp); ret
static u64_t f_read(const char *fn, reg_t *f) {
    FILE* fp = fopen(fn, "rb"); if (fp EqN) ret 3;
    F_SZ(fn, &(f->s)) { RET 3; }
    MALLOC(f->p, f->s) { RET 1; }
    if (fread(f->p, 1, f->s, fp) != f->s) { free(f->p); RET 3; }
    RET 0;
}

typedef struct task_t { u8_t *p, *f; u64_t w, h;                 } task_t;
typedef struct data_t { PTHSPT sp; task_t *t, *t_end; s63_t bpr; } data_t;

#define TILE (444 * 444)
static task_t* compute_props_of_each_tile(u64_t *N, const pm_t *pm) {
    
    task_t *r, *t; const u8_t *p = pm->p; const u64_t bpr = pm->w * 3;
    u64_t w0, w1, w, h0, h1, h, Nw, Nh;
    
    if (pm->s <= TILE * 3) Nw = Nh = 1, w0 = pm->w, h0 = pm->h; else {

        if (pm->w < 444) w1 = w = pm->w, h1 = h = TILE / w;
        else if (pm->h < 444) h1 = h = pm->h, w1 = w = TILE / h;
        else w1 = w = h1 = h = 444;

        u64_t _w = w / 2, _h = h / 2, r;

        Nw = pm->w / w, w0 = w + (r = pm->w % w); if (r > _w) Nw++, w0 = (w1 = w0 / 2) + (w0 & 1);
        Nh = pm->h / h, h0 = h + (r = pm->h % h); if (r > _h) Nh++, h0 = (h1 = h0 / 2) + (h0 & 1);

    }
    
    MALLOC(r = t, sizeof(task_t) * (*N = Nw * Nh)) ret t;
    
    fin(Nh) {
        fiN(j, Nw) {
            t->w = (j == 0) ? w0 : (j == 1 ? w1 : w);
            t->h = (i == 0) ? h0 : (i == 1 ? h1 : h);
            t->p = (u8_t *)p; p += t->w * 3;
            t++;
        }
        p = pm->p + h0 * bpr + (i > 0 ? h1 * bpr + h * bpr * (i - 1) : 0);
    }
    
    ret r;
}

#define STEP 4
#if STEP < 2
#error STEP < 2
#endif

#define pix_toU_(from) \
    fin(3) { from[i] = (s7_t)from[i]; from[i] = (from[i] << 1) ^ (from[i] >> 31); }
#define numB_(to, from) to = from[0] | from[1] | from[2]; to = numBit(to);
static void predict_predictor_rgb(_Bool *Y, _Bool *G, const task_t *t, const s63_t bpr) {
    
    u32_t F1 = 0, F2 = 0, F3 = 0, F4 = 0;
    const s63_t BPR = bpr * STEP, w = (t->w - t->w % STEP) * 3, W = BPR - w;
    const u8_t *p = t->p + (BPR - bpr) + (STEP - 1) * 3,
               *const P = p + BPR * (t->h / STEP),
               *L = p + w;
    
    for (; p != P; L += BPR, p += W) {
        for (; p != L; p += STEP * 3) {
            s31_t pix1[3], pix2[3], pix3[3], pix4[3], nl1, nl2, nl3, nl4;
            s31_t *pix[] = { pix1, pix2, pix3, pix4 };
            fin(3)
                pix2[i] = pix1[i] = p[i] - ((p[i - A] + p[i - bpr] + 1) >> 1),
                pix4[i] = pix3[i] = 
                    p[i] - ((((3 * p[i - A] + 3 * p[i - bpr]) - 2 * p[(i - bpr) - A]) + 2) >> 2);
            pix2[0] -= pix2[1], pix2[2] -= pix2[1]; pix4[0] -= pix4[1], pix4[2] -= pix4[1];
            fiN(j, 4) { s31_t *px = pix[j]; pix_toU_(px); }
            numB_(nl1, pix1); numB_(nl2, pix2); numB_(nl3, pix3); numB_(nl4, pix4);
            F1 += nl1, F2 += nl2, F3 += nl3, F4 += nl4;
        }
    }
    
    u64_t m = 0, r = F1;
    if (F2 < r) m = 1, r = F2;
    if (F3 < r) m = 2, r = F3;
    if (F4 < r) m = 3, r = F4;
    
    *Y = m & 2, *G = m & 1;
}

typedef struct bitstream_t { u64_t bc; u32_t l, *_p, *DP; } bitstream_t;
#define N_INIT u64_t N; task_t *t = compute_props_of_each_tile(&N, pm); \
    if (t == NULL) ret 1; if (T < 1 || T > 256) T = T_MAX; if (N < T) T = N; \
    data_t d = { .t = t, .t_end = t + N, .bpr = A * pm->w }; PTHSPI(&d.sp);
#define GET_TASK PTHSPL(&d->sp); \
    if (d->t == d->t_end) { PTHSPU(&d->sp); goto e; } t = d->t++; PTHSPU(&d->sp);

#define PROB_BITS  14
#define PROB_SCALE (1UL << PROB_BITS)
#define RANS64_L   (1UL << 31)

typedef struct Rans64EncSymbol {
    u64_t rcp_freq;
    u32_t freq, bias, cmpl_freq, rcp_shift;
} Rans64EncSymbol;

static void encode_stream(u32_t *F, const u64_t N, const u64_t nBit, u8_t **st,
                          const u64_t st_size, u8_t **_res, bitstream_t *b) {
    
    u32_t s = 0, cum[N + 1], *res = (u32_t *)*_res; cum[0] = 0;
    if (st_size == 0) { *_res = *st = (u8_t *)(--res), *res = 4; return; }
    fin(N) cum[i + 1] = cum[i] + F[i], s += (_Bool)F[i];
    if (s == 1) {
        *(--res) = st_size + (*(*st - 1) << 24),
       *_res = *st = (u8_t *)(--res), *res = 8 + (1 << 24); return;
    }
    fix(1, N + 1, 1) cum[i] = (cum[i] * PROB_SCALE) / st_size;
    fin(N) {
        if (F[i] && cum[i + 1] == cum[i]) {
            u32_t m = ~0u, best_steal, f;
            fiN(j, N) { f = cum[j + 1] - cum[j]; if (f > 1 && f < m) m = f, best_steal = j; }
            if (best_steal < i) fiX(j, best_steal + 1, i + 1, 1) cum[j]--;
            else fiX(j, i + 1, best_steal + 1, 1) cum[j]++;
        }
    }
    fin(N) F[i] = cum[i + 1] - cum[i];
    
    Rans64EncSymbol e[N]; fin(N) {
        
        e[i].freq = F[i];
        e[i].cmpl_freq = (1 << PROB_BITS) - F[i];
        
        if (F[i] < 2) {
            
            e[i].rcp_freq = ~0UL;
            e[i].rcp_shift = 0;
            e[i].bias = cum[i] + ((1 << PROB_BITS) - 1);
            
        } else {
            
            u32_t shift = 0; while (F[i] > (1u << shift)) shift++;
            u64_t x0, x1, t0, t1;
    
            x0 = F[i] - 1;
            x1 = 1UL << (shift + 31);
    
            t1 = x1 / F[i];
            x0 += (x1 % F[i]) << 32;
            t0 = x0 / F[i];
    
            e[i].rcp_freq = t0 + (t1 << 32);
            e[i].rcp_shift = shift - 1;
    
            e[i].bias = cum[i];
            
        }
    }
    
    const u8_t *x = *st, *const X = x - st_size;
    u64_t state0 = RANS64_L, state1 = RANS64_L;
    
    if (st_size & 1) {
        
        const Rans64EncSymbol *const sym0 = e + *(--x);
        
        u64_t q0 = (uint64_t)(((unsigned __int128)state0 * sym0->rcp_freq) >> 64);
        state0 += sym0->bias + (q0 >> sym0->rcp_shift) * sym0->cmpl_freq;
        
    }
    
    for (x -= 2; x >= X; x -= 2) {
        
        const Rans64EncSymbol *const sym1 = e + x[1];
        const Rans64EncSymbol *const sym0 = e + x[0];
        
        if (state1 >= ((RANS64_L >> PROB_BITS) << 32) * sym1->freq)
            *(--res) = (u32_t)state1, state1 >>= 32;
        if (state0 >= ((RANS64_L >> PROB_BITS) << 32) * sym0->freq)
            *(--res) = (u32_t)state0, state0 >>= 32;
        
        u64_t q1 = (uint64_t)(((unsigned __int128)state1 * sym1->rcp_freq) >> 64);
        state1 += sym1->bias + (q1 >> sym1->rcp_shift) * sym1->cmpl_freq;
        
        u64_t q0 = (uint64_t)(((unsigned __int128)state0 * sym0->rcp_freq) >> 64);
        state0 += sym0->bias + (q0 >> sym0->rcp_shift) * sym0->cmpl_freq;
        
    }
    
    *(u64_t *)(res -= 2) = state1, *(u64_t *)(res -= 2) = state0;
    
    u32_t tab_size = (N - s) + s * (PROB_BITS + 1);
    _Bool compr_tab = tab_size < N * PROB_BITS; tab_size = compr_tab ? tab_size : N * PROB_BITS;
    
    if (tab_size + 8 * (*_res - (u8_t *)res) >= nBit * st_size) {
        for (u8_t *X = *st, *x = X - st_size; x < X; x++) { BC(nBit, *x); FLUSH; }
        res = (u32_t *)*_res, *(--res) = st_size,
        *_res = *st = (u8_t *)(--res), *res = 8 + (2 << 24); return;
    }
    
    if (compr_tab) fin(N) { if (F[i]) BC(PROB_BITS + 1, F[i] + PROB_SCALE); else BC(1, 0); FLUSH; }
    else fin(N) { BC(PROB_BITS, F[i]); FLUSH; }
    *(--res) = st_size, *st = (u8_t *)(--res),
    *res = (*_res - (u8_t *)res) + ((3 + compr_tab) << 24), *_res = (u8_t *)res;
}

static PTHTF(enc_1_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u32_t F[9][256]; u8_t *f, *cx[9], *st[9], *xp[19]; MALLOC(xp[0], (u64_t)2e7) ret NULL;
    fix(1, 19, 1) xp[i] = xp[i - 1] + 500000 * (i < 12 ? 1 : (i < 18 ? 3 : 4));

t:  GET_TASK
    
    memset(F, 0, 4 * 9 * 256);
    memcpy(cx, xp, 8 * 9);
    memcpy(st + 1, xp + 9, 8 * 8); xp[17] = xp[18];
    bitstream_t b = { ._p = (u32_t *)(xp[18] + 4), .l = 24 };
    fin(3) b.bc = (b.bc << 8) + t->p[i];
    const u64_t W = bpr - t->w * 3; s31_t pix[3], nl; u32_t pl = 0;
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3; p += A;
    _Bool Y = 0, G = 0; if (t->w >= STEP && t->h >= STEP) predict_predictor_rgb(&Y, &G, t, bpr);
    
    for (; p != L;) { ENC(p1x); } p += W;
    for (; p != P;) { ENC(p1y); for (L += bpr; p != L;) { ENC4; } p += W; }
    
    fin(9) encode_stream(F[0] + i * 16, 9, 4, cx + i, cx[i] - xp[i], xp + 17, &b);
    fix(1, 9, 1) encode_stream(F[i], 1 << i * (i < 3 ? 3 : 1), i * (i < 3 ? 3 : 1),
                               st + i, st[i] - xp[i + 8], xp + 17, &b);
    
    if (b.l > 0) *(b._p)++ = b.bc << (32 - b.l);
    u64_t bsz = *(u32_t *)(xp[18]) = (u8_t *)(b._p) - xp[18], rsz = xp[18] - xp[17], tsz;

    if (bsz + rsz >= (tsz = t->w * t->h * 3)) {
        MALLOC(t->f = f, tsz + 4) goto e; *(u32_t *)f = tsz + 4; f += 4;
        for (p = t->p; p != P; p += bpr, f += t->w * 3) memcpy(f, p, t->w * 3); goto t;
    }
    
    MALLOC(t->f = f, bsz + rsz + 4) goto e; u32_t m = (1 << 26) + (Y << 25) + (G << 24);
    *(u32_t *)f = (bsz + rsz + 4) + m; f += 4; memcpy(f, xp[18], bsz); f += bsz;
    fin(9)       { memcpy(f, cx[i], tsz = *(u32_t *)(cx[i]) & bitmask(24)); f += tsz; }
    fix(1, 9, 1) { memcpy(f, st[i], tsz = *(u32_t *)(st[i]) & bitmask(24)); f += tsz; }
    
    goto t; e: free(xp[0]); return NULL;
}

#define G_p1x(to, operator) to = to operator p[0 - A];
#define G_p1y(to, operator) to = to operator p[0 - bpr];
#define G_p2a(to, operator) to = to operator ((p[0 - A] + p[0 - bpr] + 1) >> 1);
#define G_p3a(to, operator) to = \
    to operator ((((3 * p[0 - A] + 3 * p[0 - bpr]) - 2 * p[(0 - bpr) - A]) + 2) >> 2);
#define G_toU(from) from = (s7_t)from; from = (from << 1) ^ (from >> 31);
#define ENC_GRAY(pr) \
    if ((pix_0 = p[0]) != p[1] || p[1] != p[2]) ret 0; pr(pix_0, -); \
    G_toU(pix_0); fin(4) F[i][*st[i]++ = pix_0]++; p += A;

static _Bool is_grayscale(const s63_t bpr, task_t *t, u8_t *const *const xp) {
    
    if (t->p[0] != t->p[1] || t->p[1] != t->p[2]) ret 0; u8_t *res[4], *st[4], *f;
    memcpy(st, xp, sizeof(u8_t *) * 4); memcpy(res, xp + 5, sizeof(u8_t *) * 4);
    bitstream_t b[4] = {
        { ._p = (u32_t *)(xp[15] + 4), .l = 8, .bc = t->p[0] },
        { ._p = (u32_t *)(xp[16] + 4), .l = 8, .bc = t->p[0] },
        { ._p = (u32_t *)(xp[17] + 4), .l = 8, .bc = t->p[0] },
        { ._p = (u32_t *)(xp[18] + 4), .l = 8, .bc = t->p[0] }
    };
    const u64_t W = bpr - t->w * 3; u32_t F[4][256] = {};
    s31_t pix_0, pix_1, pix_2, pix_3;
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * A; p += A;
    
    for (; p != L;) { ENC_GRAY(G_p1x); } p += W;
    for (; p != P;) { ENC_GRAY(G_p1y); for (L += bpr; p != L; p += A) {
        if ((pix_0 = pix_1 = pix_2 = pix_3 = p[0]) != p[1] || p[1] != p[2]) ret 0;
        G_p1x(pix_0, -); G_toU(pix_0); F[0][*st[0]++ = pix_0]++;
        G_p1y(pix_1, -); G_toU(pix_1); F[1][*st[1]++ = pix_1]++;
        G_p2a(pix_2, -); G_toU(pix_2); F[2][*st[2]++ = pix_2]++;
        G_p3a(pix_3, -); G_toU(pix_3); F[3][*st[3]++ = pix_3]++;
    } p += W; }
    
    size_t bsz = 1e6, rsz = 1e6, tsz, m = 0;
    fin(4) {
        encode_stream(F[i], 256, 8, st + i, st[i] - xp[0 + i], res + i, b + i);
        if (b[i].l > 0) *(b[i]._p)++ = b[i].bc << (32 - b[i].l);
        size_t _bsz = *(u32_t *)(xp[15 + i]) = (u8_t *)(b[i]._p) - xp[15 + i], \
               _rsz = xp[5 + i] - res[i];
            if (_bsz + _rsz < bsz + rsz) m = i, bsz = _bsz, rsz = _rsz;
    }
    
    if (bsz + rsz >= (tsz = t->w * t->h)) {
        MALLOC(t->f = f, tsz + 4) ret 0; *(u32_t *)f = tsz + 4 + (5 << 27); f += 4;
        for (p = t->p, L = p + t->w * A; p != P; p += W, L += bpr)
            for (; p != L; p += A) *f++ = *p;
        ret 1;
    }
    
    MALLOC(t->f = f, bsz + rsz + 4) ret 0;
    *(u32_t *)f = (bsz + rsz + 4) + (2 << 28) + (m << 24); f += 4;
    memcpy(f, xp[15 + m], bsz); f += bsz; memcpy(f, res[m], rsz);
    ret 1;
}

static PTHTF(enc_2_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *f, *cx[9], *st[9], *xp[19]; MALLOC(xp[0], (u64_t)2e7) ret NULL;
    fix(1, 19, 1) xp[i] = xp[i - 1] + 500000 * (i < 12 ? 1 : (i < 18 ? 3 : 4));

t:  GET_TASK
    
    if (is_grayscale(bpr, t, xp)) goto t;
    
    u32_t F[9][256] = {};
    memcpy(cx, xp, 8 * 9); memcpy(st + 1, xp + 9, 8 * 8); xp[17] = xp[18];
    bitstream_t b = { ._p = (u32_t *)(xp[18] + 4), .l = 24 };
    fin(3) b.bc = (b.bc << 8) + t->p[i];
    const u64_t W = bpr - t->w * 3; s31_t pix[3], nl; u32_t pl = 0;
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3; p += A;
    _Bool Y = 0, G = 0; if (t->w >= STEP && t->h >= STEP) predict_predictor_rgb(&Y, &G, t, bpr);
    
    for (; p != L;) { ENC(p1x); } p += W;
    for (; p != P;) { ENC(p1y); for (L += bpr; p != L;) { ENC4; } p += W; }
    
    fin(9) encode_stream(F[0] + i * 16, 9, 4, cx + i, cx[i] - xp[i], xp + 17, &b);
    fix(1, 9, 1) encode_stream(F[i], 1 << i * (i < 3 ? 3 : 1), i * (i < 3 ? 3 : 1),
                               st + i, st[i] - xp[i + 8], xp + 17, &b);
    
    if (b.l > 0) *(b._p)++ = b.bc << (32 - b.l);
    u64_t bsz = *(u32_t *)(xp[18]) = (u8_t *)(b._p) - xp[18], rsz = xp[18] - xp[17], tsz;

    if (bsz + rsz >= (tsz = t->w * t->h * 3)) {
        MALLOC(t->f = f, tsz + 4) goto e; *(u32_t *)f = tsz + 4; f += 4;
        for (p = t->p; p != P; p += bpr, f += t->w * 3) memcpy(f, p, t->w * 3); goto t;
    }
    
    MALLOC(t->f = f, bsz + rsz + 4) goto e; u32_t m = (1 << 28) + (Y << 25) + (G << 24);
    *(u32_t *)f = (bsz + rsz + 4) + m; f += 4; memcpy(f, xp[18], bsz); f += bsz;
    fin(9)       { memcpy(f, cx[i], tsz = *(u32_t *)(cx[i]) & bitmask(24)); f += tsz; }
    fix(1, 9, 1) { memcpy(f, st[i], tsz = *(u32_t *)(st[i]) & bitmask(24)); f += tsz; }
    
    goto t; e: free(xp[0]); return NULL;
}

static _Bool ___encode(u64_t T, const u64_t mode, const pm_t *pm, const char *const filename) {
    
    TIME_PAIR; TIME_GET_START;
    
    if (pm->w > (1 << 24) || !pm->w || mode > 2 ||
        pm->h > (1 << 24) || !pm->h || pm->w * pm->h * 3 != pm->s || pm->p == NULL) ret 1;
    
    FILE *of = fopen(filename, "wb"); if (of EqN) ret 1;
    u32_t h[2] = { (pm->w - 1) | (mode << 24), pm->h - 1 };
    if (fwrite(h, 1, 8, of) != 8) ret 1;
    
    if (mode == 0) ret fwrite(pm->p, 1, pm->s, of) != pm->s || fclose(of);
    
    N_INIT; void* enc_th[] = { NULL, enc_1_th, enc_2_th };
    if (spawn_and_wait(T, &d, 0, enc_th[mode])) ret 1;
    
    TIME_GET_STOP; u64_t ns = TIME_DIFF_NS;
    pf("encode, %3d thread%c: %5lu MPx/s\n",
       (int)T_MAX, T_MAX > 1 ? 's' : ' ', (u64_t)((1e9 / ns) * (pm->s / 3e6)));
    
    u64_t x = 0;
    fin(N) {
        u32_t sz = *(u32_t *)(t[i].f) & bitmask(24);
        if (fwrite(t[i].f, 1, sz, of) != sz) ret 1;
        x += sz;
    }
    
    if (x >= pm->s) {
        if (fclose(of)) ret 1;
        of = fopen(filename, "wb"); if (of EqN) ret 1;
        h[0] &= bitmask(24);
        if (fwrite(h, 1, 8, of) != 8) ret 1;
        if (fwrite(pm->p, 1, pm->s, of) != pm->s) ret 1;
    }
    
    ret (_Bool)fclose(of);
}

_Bool png_from_pixmap(const u64_t mode, const pm_t *pm, const char *const filename) {
    
    ret ___encode(T_MAX, mode, pm, filename);
}

static void decode_stream(const u8_t *f, const u64_t N, const u64_t nBit, u8_t *x, bitstream_t *b) {
    
    const u32_t *res = (u32_t *)f,
    *const R = (u32_t *)(f + (*res & bitmask(24))), type = *res++ >> 24;
    
    switch (type) {
        case  1: memset(x, *res >> 24, *res & bitmask(24)); ret;
        case  2: fin(*res) { FILL; x[i] = BCF(nBit); } ret;
        case  3:
        case  4: break;
        default: ret;
    }
    
    u32_t st_size = *res++, F[N], cum[N + 1]; cum[0] = 0;
    if (type == 3) fin(N) { FILL; F[i] = BCF(PROB_BITS); }
    else fin(N) { FILL; if ((F[i] = BCF(1))) F[i] = BCF(PROB_BITS); }
    fin(N) cum[i + 1] = cum[i] + F[i];
    u8_t cum2sym[PROB_SCALE]; fin(N) memset(cum2sym + cum[i], i, F[i]);
    
    if (res + 4 > R) ret; u64_t state0 = *(u64_t *)res, state1 = *(u64_t *)(res + 2); res += 4;
    
    for (u8_t *const X = x + (st_size & ~1u); x < X; x += 2) {
        
        x[0] = cum2sym[state0 & bitmask(PROB_BITS)];
        x[1] = cum2sym[state1 & bitmask(PROB_BITS)];
    
        state0 = (F[x[0]] * (state0 >> PROB_BITS) + (state0 & bitmask(PROB_BITS))) - cum[x[0]];
        state1 = (F[x[1]] * (state1 >> PROB_BITS) + (state1 & bitmask(PROB_BITS))) - cum[x[1]];
    
        if (state0 < RANS64_L) { state0 = (state0 << 32) | *res; if (res < R) res++; }
        if (state1 < RANS64_L) { state1 = (state1 << 32) | *res; if (res < R) res++; }
        
    }
    
    if (st_size & 1) *x = cum2sym[state0 & bitmask(PROB_BITS)];
}

#define DEC(pr) \
    nl = *cx[nl]++; \
    switch (nl) { \
        case  0: fin(3) pix[i] = 0; break; \
        case  1: pix[2] = *st[1]++, \
                 pix[0] = (pix[2] >> 2), pix[1] = (pix[2] >> 1) & 1, pix[2] &= 1; break; \
        case  2: pix[2] = *st[2]++, \
                 pix[0] = (pix[2] >> 4), pix[1] = (pix[2] >> 2) & 3, pix[2] &= 3; break; \
        default: fin(3) pix[i] = *st[nl]++; \
    } pix_toS; pr(+); pix_copy(p, pix); p += A;
#define DEC4 \
    nl = *cx[nl]++; \
    switch (nl) { \
        case  0: fin(3) pix[i] = 0; break; \
        case  1: pix[2] = *st[1]++, \
                 pix[0] = (pix[2] >> 2), pix[1] = (pix[2] >> 1) & 1, pix[2] &= 1; break; \
        case  2: pix[2] = *st[2]++, \
                 pix[0] = (pix[2] >> 4), pix[1] = (pix[2] >> 2) & 3, pix[2] &= 3; break; \
        default: fin(3) pix[i] = *st[nl]++; \
    } pix_toS; \
    if (G) pix[0] += pix[1], pix[2] += pix[1]; if (Y) p3a(+); else p2a(+); \
    pix_copy(p, pix); p += A;


static PTHTF(dec_1_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *cx[9], *st[9], *xp[17]; MALLOC(xp[0], (u64_t)2e7) ret NULL;
    fix(1, 17, 1) xp[i] = xp[i - 1] + 500000 * (i < 12 ? 1 : 3);
    
t:  GET_TASK
    
    const u8_t *f = t->f; u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3;
    u32_t t_size = *(u32_t *)f & bitmask(24), m = f[3]; f += 4;
    if (!m) { for (; p != P; p += bpr, f += t->w * 3) memcpy(p, f, t->w * 3); goto t; }
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    fin(3) t->p[i] = (b.bc >> (b.l -= 8)) & 255;
    memcpy(cx, xp, 8 * 9); memcpy(st + 1, xp + 9, 8 * 8);
    
    fin(9) { decode_stream(f, 9, 4, cx[i], &b); f += *(u32_t *)f & bitmask(24); }
    fix(1, 9, 1) {
        decode_stream(f, 1 << i * (i < 3 ? 3 : 1), i * (i < 3 ? 3 : 1), st[i], &b);
        f += *(u32_t *)f & bitmask(24);
    }
    
    s31_t pix[3], nl = 0; const u64_t W = bpr - t->w * 3; p += A; _Bool Y = m & 2, G = m & 1;
    
    for (; p != L;) { DEC(p1x); } p += W;
    for (; p != P;) { DEC(p1y); for (L += bpr; p != L;) { DEC4; } p += W; }
    
    goto t;  e: free(xp[0]); return NULL;
}

#define DEC_GRAY(pr) \
    pix = *st++; pix = (pix >> 1) ^ (-(pix & 1)); p[0] = p[1] = p[2] = pr(pix, +); p += A;

static void if_grayscale(const s63_t bpr, task_t *t, u8_t *const *const xp) {
    
    u8_t *st = xp[0], *f = t->f + 4, m = t->f[3] & bitmask(2);
    
    u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * A;
    const u64_t W = bpr - t->w * A; s31_t pix;
    
    if (t->f[3] & (1 << 3)) {
        for (; p != P; p += W, L += bpr)
            for (; p != L; p += A) p[0] = p[1] = p[2] = *f++;
        ret;
    }
    
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    p[0] = p[1] = p[2] = (b.bc >> (b.l -= 8)) & 255; p += A;
    
    decode_stream(f, 256, 8, st, &b);
    
    for (; p != L;) { DEC_GRAY(G_p1x); } p += W;
    for (; p != P;) { DEC_GRAY(G_p1y); for (L += bpr; p != L; p += A) {
        pix = *st++; pix = (pix >> 1) ^ (-(pix & 1));
        switch (m) {
            case 0: G_p1x(pix, +); break;
            case 1: G_p1y(pix, +); break;
            case 2: G_p2a(pix, +); break;
            case 3: G_p3a(pix, +);
        }
        p[0] = p[1] = p[2] = pix;
    } p += W; }
    
}

static PTHTF(dec_2_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *cx[9], *st[9], *xp[17]; MALLOC(xp[0], (u64_t)2e7) ret NULL;
    fix(1, 17, 1) xp[i] = xp[i - 1] + 500000 * (i < 12 ? 1 : 3);
    
t:  GET_TASK
    
    const u8_t *f = t->f; u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3;
    u32_t t_size = *(u32_t *)f & bitmask(24), m = f[3]; f += 4;
    if (!m) { for (; p != P; p += bpr, f += t->w * 3) memcpy(p, f, t->w * 3); goto t; }
    if ((m >> 4) == 2) { if_grayscale(bpr, t, xp); goto t; }
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    fin(3) t->p[i] = (b.bc >> (b.l -= 8)) & 255;
    memcpy(cx, xp, 8 * 9); memcpy(st + 1, xp + 9, 8 * 8);
    
    fin(9) { decode_stream(f, 9, 4, cx[i], &b); f += *(u32_t *)f & bitmask(24); }
    fix(1, 9, 1) {
        decode_stream(f, 1 << i * (i < 3 ? 3 : 1), i * (i < 3 ? 3 : 1), st[i], &b);
        f += *(u32_t *)f & bitmask(24);
    }
    
    s31_t pix[3], nl = 0; const u64_t W = bpr - t->w * 3; p += A; _Bool Y = m & 2, G = m & 1;
    
    for (; p != L;) { DEC(p1x); } p += W;
    for (; p != P;) { DEC(p1y); for (L += bpr; p != L;) { DEC4; } p += W; }
    
    goto t;  e: free(xp[0]); return NULL;
}

static _Bool ___decode(u64_t T, const reg_t *r, pm_t *pm) {
    
    u64_t x = 8, mode; const u32_t *h = (const u32_t *)(r->p);
    pm->w = (h[0] & bitmask(24)) + 1, pm->h = (h[1] & bitmask(24)) + 1;
    if ((mode = h[0] >> 24) > 2) ret 1; pm->s = 3 * pm->w * pm->h;
    MALLOC(pm->p, pm->s) ret 1; if (!mode) ret !memcpy(pm->p, r->p + 8, pm->s);
    N_INIT; void* dec_th[] = { NULL, dec_1_th, dec_2_th };
    fin(N) t[i].f = r->p + x, x += *(u32_t *)t[i].f & bitmask(24);
    ret spawn_and_wait(T, &d, 0, dec_th[mode]);
}

static _Bool decode(const reg_t *r, pm_t *pm) {

    ret ___decode(T_MAX, r, pm);
}

_Bool png_from_jpg(const char *const jpg, const char *const png) {
    
    puts("Not Implemented.");
    
    ret 1;
}

#define PNG_COMPRESSION_TYPE_FAST 1
#define PNG_COMPRESSION_TYPE_SLOW 2
#define PNG_COMPRESSION_TYPE_EXJPEG 3
#define PNG_COMPRESSION_TYPE_UNCOMPRESSED 7

MAIN_ARGS {
    
    if (argc != 4 || argv[1][0] != '-' || strlen(argv[1]) != 2) goto h;
    pm_t pm; reg_t r; u64_t ns; TIME_PAIR;  
    
    switch (argv[1][1]) {
    
        case '0': puts("Not Implemented."); ret 1;
        case '1':
        case '2':
        case '7':
            
            F_RD(argv[2], &r) ret 1;
            if (sscanf((char *)(r.p + 3), "%lu %lu", &pm.w, &pm.h) != 2) ret 1;
            pm.s = 3 * pm.w * pm.h, pm.p = r.p + (r.s - pm.s);
            if (png_from_pixmap(argv[1][1] == '7' ? 0 : argv[1][1] - '0', &pm, argv[3])) ret 1;
            break;
        
        case '3': ret (int)png_from_jpg(argv[2], argv[3]);
        case 'd':
            
            F_RD(argv[2], &r) ret 1;
            TIME_DIFF_EXEC(if (decode(&r, &pm)) ret 1, NS, ns);
            pf("decode, %3d thread%c: %5lu MPx/s\n",
               (int)T_MAX, T_MAX > 1 ? 's' : ' ', (u64_t)((1e9 / ns) * (pm.s / 3e6)));
            char P6[100]; int l = sprintf(P6, "P6\n%lu %lu\n255\n", pm.w, pm.h);
            FILE *of = fopen(argv[3], "wb"); if (of EqN) ret 1;
            if (fwrite(P6, 1, l, of) != l || fwrite(pm.p, 1, pm.s, of) != pm.s || fclose(of)) ret 1;
            break;
        
        default: goto h;
    }
    
    ret 0; h: pf("\n"

    "encode: ./xpng -[0127] example.ppm  example.xpng\n"
    "        ./xpng -3      example.jpg  example.xpng\n"
    "decode: ./xpng -d      example.xpng example.ppm\n"

    "\n"); ret 1;
}
