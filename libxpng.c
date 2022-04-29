#include "xpng.h"

#if T_MAX < 1
#error
#endif

#define RGB  3
#define RGBA 4

#define BITSTREAM_WRITE(b, C, B) b.bc = (b.bc << (C)) | (B), b.l += (C)
#define BITSTREAM_READ(b, C) ((b.bc >> (b.l -= (C))) & BITMASK(C))
#define BITSTREAM_FLUSH(b) if (b.l >= 32) *(b._p)++ = b.bc >> (b.l -= 32)
#define BITSTREAM_END(b) if (b.l > 0) *(b._p)++ = b.bc << (32 - b.l)
#define BITSTREAM_FILL(b) if (b.l < 32) { b.bc <<= 32, b.l += 32; if (b._p < b.DP) b.bc += *(b._p)++; }

#define BC(C, B) b->bc = (b->bc << (C)) + (B), b->l += (C)
#define BCF(C) ((b->bc >> (b->l -= (C))) & BITMASK(C))
#define FILL  if (b->l  < 32) { b->bc <<= 32, b->l += 32; if (b->_p < b->DP) b->bc += *(b->_p)++; }
#define FLUSH if (b->l >= 32) *(b->_p)++ = b->bc >> (b->l -= 32)

#define numBit(v) ((1 + (__builtin_clz(v) ^ 31)) & (-(v) >> 31))
#define pix_toU fin(3) { pix[i] = (s7_t)pix[i]; pix[i] = (pix[i] << 1) ^ (pix[i] >> 31); }
#define pix_toS fin(3) pix[i] = (pix[i] >> 1) ^ (-(pix[i] & 1));
#define pix_copy(to, from) fin(3) to[i] = from[i];
#define p1x_(pxsz) p[i - pxsz]
#define p1x(operator) fin(3) pix[i] = pix[i] operator p1x_(RGB);
#define p1y(operator) fin(3) pix[i] = pix[i] operator p[i - bpr];
#define p2a_(pxsz) ((p[i - pxsz] + p[i - bpr] + 1) >> 1)
#define p2a(operator) fin(3) pix[i] = pix[i] operator p2a_(RGB);
#define p3a_(pxsz) ((((3 * p[i - pxsz] + 3 * p[i - bpr]) - 2 * p[(i - bpr) - pxsz]) + 2) >> 2)
#define p3a(operator) fin(3) pix[i] = pix[i] operator p3a_(RGB)
#define numB(from) nl = from[0] | from[1] | from[2]; nl = numBit(nl);

#define M1ENC_ \
    pix_toU; numB(pix); F[(pl << 4) + nl]++, pl = *cx[pl]++ = nl; \
    if (nl) { \
        const int n2 = nl << 1, n3 = n2 + nl; \
        BITSTREAM_WRITE(k, n3, (pix[0] << n2) | (pix[1] << nl) | pix[2]); \
        BITSTREAM_FLUSH(k); \
    } p += RGB;
#define M1ENC(pr) pix_copy(pix, p); pr(-); M1ENC_
#define M1ENC4 \
    pix_copy(pix, p); if (Y) p3a(-); else p2a(-); \
    if (G) pix[0] -= pix[1], pix[2] -= pix[1]; M1ENC_

#define ENC_ \
    pix_toU; numB(pix); F[0][(pl << 4) + nl]++, pl = *cx[pl]++ = nl; \
    switch (nl) { \
        case  0: break; \
        case  1: F[1][*st[1]++ = (pix[0] << 2) | (pix[1] << 1) | pix[2]]++; break; \
        case  2: F[2][*st[2]++ = (pix[0] << 4) | (pix[1] << 2) | pix[2]]++; break; \
        default: fin(3) F[nl][*st[nl]++ = pix[i]]++; \
    } p += RGB;
#define ENC(pr) pix_copy(pix, p); pr(-); ENC_
#define ENC4 \
    pix_copy(pix, p); if (Y) p3a(-); else p2a(-); \
    if (G) pix[0] -= pix[1], pix[2] -= pix[1]; ENC_

typedef struct task_t { u8_t *p, *f; u64_t w, h;                 } task_t;
typedef struct data_t { PTHSPT sp; task_t *t, *t_end; s63_t bpr; } data_t;

#define TILE (444 * 444)
__attribute__ ((noinline))
static task_t* compute_props_of_each_tile(u64_t *N, const xpng_t *pm) {
    
    task_t *r, *t; const u8_t *p = pm->p; const u64_t bpr = pm->w * RGB;
    u64_t w0, w1, w, h0, h1, h, Nw, Nh;
    
    if (pm->s <= TILE * RGB) Nw = Nh = 1, w0 = pm->w, h0 = pm->h; else {

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
            t->p = (u8_t *)p; p += t->w * RGB;
            t++;
        }
        p = pm->p + h0 * bpr + (i > 0 ? h1 * bpr + h * bpr * (i - 1) : 0);
    }
    
    ret r;
}

#define STEP 4
#define PP_RGBX(s) \
    fin(3) { pix##s[i] = (s7_t)pix##s[i]; pix##s[i] = (pix##s[i] << 1) ^ (pix##s[i] >> 31); } \
    n##s = pix##s[0] | pix##s[1] | pix##s[2]; F##s += numBit(n##s);
__attribute__ ((noinline))
static u64_t pp_rgbx(const task_t *t, const s63_t bpr, const s63_t PXSZ) {
    
    if (STEP < 2) exit(1);
    if (t->w < STEP || t->h < STEP) ret 0;
    
    u32_t F1 = 0, F2 = 0, F3 = 0, F4 = 0;
    const s63_t BPR = bpr * STEP, w = (t->w - t->w % STEP) * PXSZ, W = BPR - w;
    const u8_t *p = t->p + (BPR - bpr) + (STEP - 1) * PXSZ,
               *const P = p + BPR * (t->h / STEP),
               *L = p + w;
    
    s31_t pix1[3], pix2[3], pix3[3], pix4[3], n1, n2, n3, n4;
    
    if (PXSZ == 3) {
        
        for (; p != P; L += BPR, p += W) {
            for (; p != L; p += STEP * RGB) {
                fin(3)
                    pix2[i] = pix1[i] = p[i] - p2a_(RGB),
                    pix4[i] = pix3[i] = p[i] - p3a_(RGB);
                pix2[0] -= pix2[1], pix2[2] -= pix2[1];
                pix4[0] -= pix4[1], pix4[2] -= pix4[1];
                PP_RGBX(1); PP_RGBX(2); PP_RGBX(3); PP_RGBX(4);
            }
        }
    
    } else {
        
        for (; p != P; L += BPR, p += W) {
            for (; p != L; p += STEP * RGBA) {
                if (p[3] == 0) continue;
                fin(3)
                    pix2[i] = pix1[i] = p[i] - p2a_(RGBA),
                    pix4[i] = pix3[i] = p[i] - p3a_(RGBA);
                pix2[0] -= pix2[1], pix2[2] -= pix2[1];
                pix4[0] -= pix4[1], pix4[2] -= pix4[1];
                PP_RGBX(1); PP_RGBX(2); PP_RGBX(3); PP_RGBX(4);
            }
        }
    
    }
    
    u64_t m = 0, r = F1;
    
    if (F2 < r) m = 1, r = F2;
    if (F3 < r) m = 2, r = F3;
    if (F4 < r) m = 3, r = F4;
    
    ret m;
}

#undef PP_RGBX
#undef STEP

typedef struct bitstream_t { u64_t bc; u32_t l, *_p, *DP; } bitstream_t;
#define N_INIT u64_t N; task_t *t = compute_props_of_each_tile(&N, pm); \
    if (t == NULL) ret 1; if (T < 1 || T > 256) T = T_MAX; if (N < T) T = N; \
    data_t d = { .t = t, .t_end = t + N, .bpr = (RGB + pm->A) * pm->w }; PTHSPI(&d.sp);
#define GET_TASK PTHSPL(&d->sp); \
    if (d->t == d->t_end) { PTHSPU(&d->sp); goto e; } t = d->t++; PTHSPU(&d->sp);

#define RANS64_L (1UL << 31)

typedef struct Rans64EncSymbol {
    u64_t rcp_freq;
    u32_t freq, bias, cmpl_freq, rcp_shift;
} Rans64EncSymbol;

static void compress_block(u32_t *F, const u64_t N, u8_t **st, const u64_t st_size,
                           u8_t **_res, bitstream_t *b, const int PROB_BITS) {
    
    const u64_t PROB_SCALE = (1UL << PROB_BITS);
    const u64_t nBit = numBit((int)(N - 1));
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

static void decompress_block(const u8_t *f, const u64_t N, u8_t *x,
                             bitstream_t *b, const int PROB_BITS) {
    
    const u64_t PROB_SCALE = (1UL << PROB_BITS);
    const u64_t nBit = numBit((int)(N - 1));
    const u32_t *res = (u32_t *)f,
    *const R = (u32_t *)(f + (*res & BITMASK(24))), type = *res++ >> 24;
    
    switch (type) {
        case  1: memset(x, *res >> 24, *res & BITMASK(24)); ret;
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
        
        x[0] = cum2sym[state0 & BITMASK(PROB_BITS)];
        x[1] = cum2sym[state1 & BITMASK(PROB_BITS)];
    
        state0 = (F[x[0]] * (state0 >> PROB_BITS) + (state0 & BITMASK(PROB_BITS))) - cum[x[0]];
        state1 = (F[x[1]] * (state1 >> PROB_BITS) + (state1 & BITMASK(PROB_BITS))) - cum[x[1]];
    
        if (state0 < RANS64_L) { state0 = (state0 << 32) | *res; if (res < R) res++; }
        if (state1 < RANS64_L) { state1 = (state1 << 32) | *res; if (res < R) res++; }
        
    }
    
    if (st_size & 1) *x = cum2sym[state0 & BITMASK(PROB_BITS)];
}

#undef RANS64_L

static PTHTF(enc_1_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *xp[12]; MALLOC(xp[0], (u64_t)105e5) ret NULL;
    fix(1, 12, 1) xp[i] = xp[i - 1] + 500000 * (i < 10 ? 1 : 4); xp[9] = xp[10];

t:  GET_TASK
    
    u32_t F[256] = {};
    u8_t *f, *cx[10]; memcpy(cx, xp, 8 * 10);
    bitstream_t b = { ._p = (u32_t *)(xp[10] + 4), .l = 24 };
    fin(3) b.bc = (b.bc << 8) + t->p[i];
    bitstream_t k = { ._p = (u32_t *)(xp[11] + 4), };
    const u64_t W = bpr - t->w * RGB; s31_t pix[3], nl; u32_t pl = 0;
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * RGB; p += RGB;
    const u64_t pr = pp_rgbx(t, bpr, RGB); const _Bool Y = pr & 2, G = pr & 1;
    
    for (; p != L;) { M1ENC(p1x); } p += W;
    for (; p != P;) { M1ENC(p1y); for (L += bpr; p != L;) { M1ENC4; } p += W; }
    
    fin(9) compress_block(F + i * 16, 9, cx + i, cx[i] - xp[i], cx + 9, &b, 14);
    
    BITSTREAM_END(b);
    BITSTREAM_END(k);
    u64_t bsz = *(u32_t *)(xp[10]) = (u8_t *)(b._p) - xp[10],
          ksz = *(u32_t *)(xp[11]) = (u8_t *)(k._p) - xp[11],
          rsz = xp[9] - cx[9], tsz, fsz = rsz + bsz + ksz;

    if (fsz >= (tsz = t->w * t->h * RGB)) {
        MALLOC(t->f = f, tsz + 4) goto e; *(u32_t *)f = tsz + 4; f += 4;
        for (p = t->p; p != P; p += bpr, f += t->w * RGB) memcpy(f, p, t->w * RGB); goto t;
    }
    
    MALLOC(t->f = f, fsz + 4) goto e; u32_t m = (1 << 26) + (pr << 24);
    *(u32_t *)f = (fsz + 4) + m; f += 4;
    memcpy(f, xp[10], bsz); f += bsz;
    fin(9) { memcpy(f, cx[i], tsz = *(u32_t *)(cx[i]) & BITMASK(24)); f += tsz; }
    memcpy(f, xp[11], ksz);
    
    goto t; e: free(xp[0]); return NULL;
}

#define G_p1x(to, operator) to = to operator p[0 - RGB];
#define G_p1y(to, operator) to = to operator p[0 - bpr];
#define G_p2a(to, operator) to = to operator ((p[0 - RGB] + p[0 - bpr] + 1) >> 1);
#define G_p3a(to, operator) to = \
    to operator ((((3 * p[0 - RGB] + 3 * p[0 - bpr]) - 2 * p[(0 - bpr) - RGB]) + 2) >> 2);
#define G_toU(from) from = (s7_t)from; from = (from << 1) ^ (from >> 31);
#define ENC_GRAY(pr) \
    if ((pix_0 = p[0]) != p[1] || p[1] != p[2]) ret 0; pr(pix_0, -); \
    G_toU(pix_0); fin(4) F[i][*st[i]++ = pix_0]++; p += RGB;

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
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * RGB; p += RGB;
    
    for (; p != L;) { ENC_GRAY(G_p1x); } p += W;
    for (; p != P;) { ENC_GRAY(G_p1y); for (L += bpr; p != L; p += RGB) {
        if ((pix_0 = pix_1 = pix_2 = pix_3 = p[0]) != p[1] || p[1] != p[2]) ret 0;
        G_p1x(pix_0, -); G_toU(pix_0); F[0][*st[0]++ = pix_0]++;
        G_p1y(pix_1, -); G_toU(pix_1); F[1][*st[1]++ = pix_1]++;
        G_p2a(pix_2, -); G_toU(pix_2); F[2][*st[2]++ = pix_2]++;
        G_p3a(pix_3, -); G_toU(pix_3); F[3][*st[3]++ = pix_3]++;
    } p += W; }
    
    size_t bsz = 1e6, rsz = 1e6, tsz, m = 0;
    fin(4) {
        compress_block(F[i], 256, st + i, st[i] - xp[0 + i], res + i, b + i, 15);
        if (b[i].l > 0) *(b[i]._p)++ = b[i].bc << (32 - b[i].l);
        size_t _bsz = *(u32_t *)(xp[15 + i]) = (u8_t *)(b[i]._p) - xp[15 + i], \
               _rsz = xp[5 + i] - res[i];
            if (_bsz + _rsz < bsz + rsz) m = i, bsz = _bsz, rsz = _rsz;
    }
    
    if (bsz + rsz >= (tsz = t->w * t->h)) {
        MALLOC(t->f = f, tsz + 4) ret 0; *(u32_t *)f = tsz + 4 + (5 << 27); f += 4;
        for (p = t->p, L = p + t->w * RGB; p != P; p += W, L += bpr)
            for (; p != L; p += RGB) *f++ = *p;
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
    const u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3; p += RGB;
    const u64_t pr = pp_rgbx(t, bpr, 3); const _Bool Y = pr & 2, G = pr & 1;
    
    for (; p != L;) { ENC(p1x); } p += W;
    for (; p != P;) { ENC(p1y); for (L += bpr; p != L;) { ENC4; } p += W; }
    
    fin(9) compress_block(F[0] + i * 16, 9, cx + i, cx[i] - xp[i], xp + 17, &b, 14);
    fix(1, 9, 1) compress_block(F[i], 1 << i * (i < 3 ? 3 : 1),
                               st + i, st[i] - xp[i + 8], xp + 17, &b, 14);
    
    if (b.l > 0) *(b._p)++ = b.bc << (32 - b.l);
    u64_t bsz = *(u32_t *)(xp[18]) = (u8_t *)(b._p) - xp[18], rsz = xp[18] - xp[17], tsz;

    if (bsz + rsz >= (tsz = t->w * t->h * 3)) {
        MALLOC(t->f = f, tsz + 4) goto e; *(u32_t *)f = tsz + 4; f += 4;
        for (p = t->p; p != P; p += bpr, f += t->w * 3) memcpy(f, p, t->w * 3); goto t;
    }
    
    MALLOC(t->f = f, bsz + rsz + 4) goto e; u32_t m = (1 << 28) + (pr << 24);
    *(u32_t *)f = (bsz + rsz + 4) + m; f += 4; memcpy(f, xp[18], bsz); f += bsz;
    fin(9)       { memcpy(f, cx[i], tsz = *(u32_t *)(cx[i]) & BITMASK(24)); f += tsz; }
    fix(1, 9, 1) { memcpy(f, st[i], tsz = *(u32_t *)(st[i]) & BITMASK(24)); f += tsz; }
    
    goto t; e: free(xp[0]); return NULL;
}

__attribute__ ((noinline))
static _Bool normalize_RGBA(xpng_t *const d) {
    
    if (d->A == 0) ret 0;
    
    u32_t *p = (void *)(d->p), *const P = (void *)(d->p + d->s);
    _Bool a = 0, inv_rgb = 0;
    for (; p < P; p++) {
        if ((*p >> 24) == 0 && *p != 0) { inv_rgb = 1; break; };
        if ((*p >> 24) != 255) a = 1;
    }
    
    if (inv_rgb) {
        
        p = (void *)(d->p); MALLOC(d->p, d->s) ret 1;
        u32_t *p2 = (void *)(d->p);
        for (; p < P; p++, p2++)
            *p2 = (*p >> 24) ? *p : 0;
        
        ret 0;
    }
    
    if (a) ret 0;
    
    p = (void *)(d->p), d->s -= (d->s / 4), d->A = 0;
    MALLOC(d->p, d->s) ret 1;
    
    u8_t *p2 = d->p, *const P2 = p2 + (d->s - 3);
    for (; p2 < P2; p2 += 3, p++)
        *(u32_t *)p2 = *p;
        
    *(u16_t *)p2 = *p, p2[2] = *p >> 16;
    
    ret 0;
}

_Bool xpng_store_T(u64_t T, const u64_t mode, const xpng_t *const pm_, const char *const fn) {

#ifdef SYNC_IO

    TIME_PAIR; TIME_GET_START;
    
    if (pm_->w > (1 << 24) || !pm_->w || !(mode == 1 || mode == 2 || mode == 7) ||
        pm_->h > (1 << 24) || !pm_->h || pm_->w * pm_->h * (3 + pm_->A) != pm_->s ||
        pm_->p == NULL) ret 1;
    
    const xpng_t pmv = *pm_, *const pm = &pmv;
    if (normalize_RGBA((xpng_t *)pm)) ret 1;
    
    u32_t h[2] = { (pm->w - 1) | (mode << 24), pm->h - 1 };
    FILE *of = fopen(fn, "wb");
    if (of == NULL || fwrite(h, 1, 8, of) != 8) ret 1;
    
    if (mode == 7) ret fwrite(pm->p, 1, pm->s, of) != pm->s || fclose(of);
    
    N_INIT; if (spawn_and_wait(T, &d, 0, (void* []){ NULL, enc_1_th, enc_2_th }[mode])) ret 1;
    
    TIME_GET_STOP; u64_t ns = TIME_DIFF_NS;
    pf("encode, %3d thread%c: %5lu MPx/s\n",
       (int)T, T > 1 ? 's' : ' ', (u64_t)((1e9 / ns) * (pm->s / 3e6)));
    
    u64_t x = 0;
    fin(N) {
        u32_t sz = *(u32_t *)(t[i].f) & BITMASK(24);
        if (fwrite(t[i].f, 1, sz, of) != sz) ret 1;
        x += sz, free(t[i].f);
    }
    
    if (x >= pm->s) {
        ((u8_t *)h)[3] = XPNG_COMPRESSION_TYPE_UNCOMPRESSED;
        if (fclose(of)
            || (of = fopen(fn, "wb")) == NULL
            || fwrite(h, 1, 8, of) != 8
            || fwrite(pm->p, 1, pm->s, of) != pm->s) ret 1;
    }
    
e:  if (pm->p != pm_->p) free(pm->p);
    
    ret (_Bool)fclose(of);

#else
#error
#endif

}

_Bool xpng_store(const u64_t mode, const xpng_t *pm, const char *const filename) {
    
    ret xpng_store_T(T_MAX, mode, pm, filename);
}

#define M1DEC_ \
    if ((nl = *cx[nl]++)) { \
        BITSTREAM_FILL(k); \
        const int n2 = nl << 1, n3 = n2 + nl; \
        pix[2] = BITSTREAM_READ(k, n3); \
        pix[0] =  pix[2] >> n2, \
        pix[1] = (pix[2] >> nl) & BITMASK(nl), \
        pix[2] &= BITMASK(nl); \
    } else fin(3) pix[i] = 0; pix_toS;
#define M1DEC(pr) M1DEC_; pr(+); pix_copy(p, pix); p += RGB;
#define M1DEC4 M1DEC_; if (G) pix[0] += pix[1], pix[2] += pix[1]; \
                       if (Y) p3a(+); else p2a(+); pix_copy(p, pix); p += RGB;

static PTHTF(dec_1_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *cx[9], *xp[9]; MALLOC(xp[0], (u64_t)45e5) ret NULL;
    fix(1, 9, 1) xp[i] = xp[i - 1] + 500000;
    
t:  GET_TASK
    
    const u8_t *f = t->f;
    u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3;
    u32_t t_size = *(u32_t *)f & BITMASK(24), m = f[3]; f += 4;
    if (!m) { for (; p != P; p += bpr, f += t->w * RGB) memcpy(p, f, t->w * RGB); goto t; }
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    fin(3) t->p[i] = BITSTREAM_READ(b, 8);
    memcpy(cx, xp, 8 * 9);
    
    fin(9) { decompress_block(f, 9, cx[i], &b, 14); f += *(u32_t *)f & BITMASK(24); }
    
    bitstream_t k = { ._p = (u32_t *)(f + 4), .l = 32 };
    k.DP = (u32_t *)(f + *(u32_t *)f), k.bc = *(k._p)++;
    
    s31_t pix[3], nl = 0;
    const u64_t W = bpr - t->w * RGB; p += RGB;
    const _Bool Y = m & 2, G = m & 1;
    
    for (; p != L;) { M1DEC(p1x); } p += W;
    for (; p != P;) { M1DEC(p1y); for (L += bpr; p != L;) { M1DEC4; } p += W; }
    
    goto t;  e: free(xp[0]); return NULL;
}

#define DEC_GRAY(pr) \
    pix = *st++; pix = (pix >> 1) ^ (-(pix & 1)); p[0] = p[1] = p[2] = pr(pix, +); p += RGB;

static void if_grayscale(const s63_t bpr, task_t *t, u8_t *const *const xp) {
    
    u8_t *st = xp[0], *f = t->f + 4, m = t->f[3] & BITMASK(2);
    
    u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * RGB;
    const u64_t W = bpr - t->w * RGB; s31_t pix;
    
    if (t->f[3] & (1 << 3)) {
        for (; p != P; p += W, L += bpr)
            for (; p != L; p += RGB) p[0] = p[1] = p[2] = *f++;
        ret;
    }
    
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    p[0] = p[1] = p[2] = (b.bc >> (b.l -= 8)) & 255; p += RGB;
    
    decompress_block(f, 256, st, &b, 15);
    
    for (; p != L;) { DEC_GRAY(G_p1x); } p += W;
    for (; p != P;) { DEC_GRAY(G_p1y); for (L += bpr; p != L; p += RGB) {
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

#define DEC_ \
    nl = *cx[nl]++; \
    switch (nl) { \
        case  0: fin(3) pix[i] = 0; break; \
        case  1: pix[2] = *st[1]++, \
                 pix[0] = (pix[2] >> 2), pix[1] = (pix[2] >> 1) & 1, pix[2] &= 1; break; \
        case  2: pix[2] = *st[2]++, \
                 pix[0] = (pix[2] >> 4), pix[1] = (pix[2] >> 2) & 3, pix[2] &= 3; break; \
        default: fin(3) pix[i] = *st[nl]++; \
    } pix_toS;
#define DEC(pr) DEC_; pr(+); pix_copy(p, pix); p += RGB;
#define DEC4 DEC_; \
    if (G) pix[0] += pix[1], pix[2] += pix[1]; if (Y) p3a(+); else p2a(+); \
    pix_copy(p, pix); p += RGB;

static PTHTF(dec_2_th) {
    
    data_t *const d = data; const s63_t bpr = d->bpr; task_t *t;
    
    u8_t *cx[9], *st[9], *xp[17]; MALLOC(xp[0], (u64_t)2e7) ret NULL;
    fix(1, 17, 1) xp[i] = xp[i - 1] + 500000 * (i < 12 ? 1 : 3);
    
t:  GET_TASK
    
    const u8_t *f = t->f; u8_t *p = t->p, *const P = p + bpr * t->h, *L = p + t->w * 3;
    u32_t t_size = *(u32_t *)f & BITMASK(24), m = f[3]; f += 4;
    if (!m) { for (; p != P; p += bpr, f += t->w * 3) memcpy(p, f, t->w * 3); goto t; }
    if ((m >> 4) == 2) { if_grayscale(bpr, t, xp); goto t; }
    bitstream_t b = { ._p = (u32_t *)(f + 4), .l = 32 };
    b.DP = (u32_t *)(f += *(u32_t *)f), b.bc = *(b._p)++;
    fin(3) t->p[i] = (b.bc >> (b.l -= 8)) & 255;
    memcpy(cx, xp, 8 * 9); memcpy(st + 1, xp + 9, 8 * 8);
    
    fin(9) { decompress_block(f, 9, cx[i], &b, 14); f += *(u32_t *)f & BITMASK(24); }
    fix(1, 9, 1) {
        decompress_block(f, 1 << i * (i < 3 ? 3 : 1), st[i], &b, 14);
        f += *(u32_t *)f & BITMASK(24);
    }
    
    s31_t pix[3], nl = 0; const u64_t W = bpr - t->w * 3; p += RGB; _Bool Y = m & 2, G = m & 1;
    
    for (; p != L;) { DEC(p1x); } p += W;
    for (; p != P;) { DEC(p1y); for (L += bpr; p != L;) { DEC4; } p += W; }
    
    goto t;  e: free(xp[0]); return NULL;
}

_Bool xpng_load_T(u64_t T, const char *const xpng, xpng_t *pm) {

#ifdef SYNC_IO

    reg_t r; F_READ(xpng, &r) ret 1; TIME_PAIR; TIME_GET_START;
    
    u64_t x = 8; const u32_t *h = (const u32_t *)(r.p);
    pm->w = (h[0] & BITMASK(24)) + 1, pm->h = (h[1] & BITMASK(24)) + 1, pm->A = (h[1] >> 24) & 1;
    const u64_t mode = h[0] >> 24;
    unless (mode == 1 || mode == 2 || mode == 7) ret 1;
    pm->s = 3 * pm->w * pm->h;
    MALLOC(pm->p, pm->s) ret 1; if (mode == 7) ret !memcpy(pm->p, r.p + 8, pm->s);
    N_INIT; fin(N) t[i].f = r.p + x, x += *(u32_t *)t[i].f & BITMASK(24);
    if (spawn_and_wait(T, &d, 0, (void* []){ NULL, dec_1_th, dec_2_th }[mode])) ret 1;
    
    TIME_GET_STOP; u64_t ns = TIME_DIFF_NS;
    pf("decode, %3d thread%c: %5lu MPx/s\n",
        (int)T, T > 1 ? 's' : ' ', (u64_t)((1e9 / ns) * (pm->s / 3e6)));
    
    ret 0;
    
#else
#error
#endif

}

_Bool xpng_load(const char *const xpng, xpng_t *pm) {

    ret xpng_load_T(T_MAX, xpng, pm);
}

_Bool xpng_from_jpg_T(u64_t T, const char *const jpg, const char *const xpng) {
    
    puts("\nNot Implemented.\n");
    
    ret 1;
}

_Bool xpng_from_jpg(const char *const jpg, const char *const xpng) {
    
    ret xpng_from_jpg_T(T_MAX, jpg, xpng);
}
