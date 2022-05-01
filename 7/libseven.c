#include "seven.h"

_Bool store_7(const xpng_t *const pm, const char *const fn) {
    
    if (pm == nil  || pm->p == nil ||
        pm->w == 0 || pm->w > (1 << 24) ||
        pm->h == 0 || pm->h > (1 << 24) ||
        pm->s != pm->w * pm->h * (3 + pm->A)) ret 1;
    
    FILE *f;
    const u32_t h[2] = { (pm->w - 1) | (7 << 24), (pm->h - 1) | (pm->A << 24) };
    ret (f = fopen(fn, "wb")) == nil ||
        fwrite(h, 1, 8, f) != 8 ||
        fwrite(pm->p, 1, pm->s, f) != pm->s ||
        fclose(f);
}

_Bool load_7(const char *const fn, xpng_t *const pm) {
    
    const u64_t in_size; F_SIZE(fn, (u64_t *)&in_size) ret 1;
    if (in_size < 11) ret 1; u32_t h[2]; FILE *f;
    if ((f = fopen(fn, "rb")) == nil || fread(h, 1, 8, f) != 8) ret 1;
    
    *pm = (xpng_t){
        .w = (h[0] & BITMASK(24)) + 1,
        .h = (h[1] & BITMASK(24)) + 1,
        .A = (h[1] >> 24) & 1,
    }, pm->s = pm->w * pm->h * (3 + pm->A);
    
    if (pm->s + 8 != in_size || (h[0] >> 24) != 7) ret 1;
    MALLOC(pm->p, pm->s) ret 1;
    
    if (fread(pm->p, 1, pm->s, f) != pm->s) ret 1;
    
    ret fclose(f);
}
