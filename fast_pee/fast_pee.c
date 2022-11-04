#include "../7/seven.h"

#define PERR(CND, STR) if (CND) { pf("\n\t%s\n\tDecoding failed.\n\n", STR); exit(1); }

void decode(const reg_t *const r, xpng_t *const pm, u64_t *const T) {
    
    if (r->s < 25 || memcmp(r->p, "RIFF", 4) || memcmp(r->p + 8, "WEBP", 4)) exit(1);
    PERR(memcmp(r->p + 12, "VP8L", 4), "Bytes 12-15 must be equal to the string 'VP8L'.");
    
    u32_t *p = (u32_t *)(r->p + 21);
    PERR((*p >> 28) & 1, "Alpha channel not supported.");
    pm->w = (*p & BITMASK(14)) + 1;
    pm->h = ((*p >> 14) & BITMASK(14)) + 1;
    MALLOC(pm->p, pm->s = 3 * pm->w * pm->h) exit(1);
    
    *T = sysconf(_SC_NPROCESSORS_ONLN);
    
    sleep(1);
}

MAIN_ARGS {
    
    reg_t r; xpng_t pm; u64_t T, ns; TIME_PAIR;
    if (argc != 3 || f_read(argv[1], &r)) ret 1;
    
    TIME_DIFF_EXEC(decode(&r, &pm, &T), NS, ns);
    pf("\n\t%3d thread%c: %5lu MPx/s\n\n",
        (int)T, T > 1 ? 's' : ' ', (u64_t)((1e9 / ns) * (pm.s / 3e6)));
    
    FILE *o; char h[20]; int l = sprintf(h, "P6\n%lu %lu\n255\n", pm.w, pm.h);
    if ((o = fopen(argv[2], "wb")) == NULL ||
        fwrite(h, 1, l, o) != l || fwrite(pm.p, 1, pm.s, o) != pm.s) ret 1;
    
    ret 0;
}
