#include "../7/seven.h"

NOINLINE static void op_r90(xpng_t *const pm) {
    
    u8_t *const a = pm->p, *b; MALLOC(b, pm->s) exit(1);
    const u64_t z = pm->A + 3, bpr = z * pm->w, bw = z * pm->h;
    
    fin(pm->h)
        fiN(j, pm->w) {
            
            u8_t *x = (a + i * bpr) + j * z,
                 *y = (b + j * bw ) + (bw - (i + 1) * z);
            
            if (z == 3)
                *(u16_t *)y = *(u16_t *)x, y[2] = x[2];
            else
                *(u32_t *)y = *(u32_t *)x;
        }
    
    free(a); pm->p = b, pm->w = pm->h, pm->h = bpr / z;
    
}

NOINLINE static void op_r270(xpng_t *const pm) {
    
}

NOINLINE static void op_mv(xpng_t *const pm) {
    
    u64_t bpr = (pm->A + 3) * pm->w;
    u8_t *a = pm->p, *b = (a + pm->s) - bpr;
    
    while (a < b) {
        
        {
            u64_t t, *x = (void *)a, *y = (void *)b;
            fin((bpr & ~7UL) / 8) t = *x, *x++ = *y, *y++ = t;
        }
        
        {
            u8_t t, *x = a + (bpr & ~7UL), *y = b + (bpr & ~7UL);
            fin(bpr & 7) t = *x, *x++ = *y, *y++ = t;
        }
        
        a += bpr, b -= bpr;
        
    }
    
}

NOINLINE static void op_mh(xpng_t *const pm) {
    
    const u64_t z = pm->A + 3, bpr = z * pm->w;
    u8_t *a = pm->p, *const A = a + pm->s;
    
    while (a < A) {
        
        u8_t *x = a, *y = (a + bpr) - z; u32_t t;
        
        if (z == 3)
        
            while (x < y) 
                t = *(u16_t *)x, t |= x[2] << 16,
                *(u16_t *)x = *(u16_t *)y, x[2] = y[2],
                *(u16_t *)y = t, y[2] = t >> 16,
                x += z,
                y -= z;
            
        else
        
            while (x < y)
                t = *(u32_t *)x,
                *(u32_t *)x = *(u32_t *)y,
                *(u32_t *)y = t,
                x += z,
                y -= z;
        
        a += bpr;
        
    }
    
}

NOINLINE static void op_mvh(xpng_t *const pm) {
    
}

NOINLINE static void op_tl(xpng_t *const pm) {
    
}

NOINLINE static void op_tr(xpng_t *const pm) {
    
}

MAIN_ARGS {
    
    if (argc != 4) goto h; xpng_t pm; if (load_7(argv[2], &pm)) ret 1; u64_t o = 10;
    
    char *s[] = { "--r90", "--r270", "--mv", "--mh", "--mvh", "--tl", "--tr" };
    fin(7) if (!strcmp(argv[1], s[i])) { o = i; break; } if (o == 10) goto h;
    (typeof(op_r90)* []){ op_r90, op_r270, op_mv, op_mh, op_mvh, op_tl, op_tr }[o](&pm);
    
    ret (int)store_7(&pm, argv[3]);
    
h:  pf("\n\t./tool --(r90|r270|mv|mh|mvh|tl|tr) src.7 res.7\n\n"); ret 1;

}
