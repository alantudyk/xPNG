#include "seven.h"
#include <png.h>

NOINLINE static _Bool normalize_RGBA(xpng_t *const d) {
    
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

MAIN_ARGS {
    
    if (argc != 4) goto h; xpng_t pm; png_image img = { .version = PNG_IMAGE_VERSION };
    
    if (strcmp(argv[1], "--to_png") == 0) goto d;
    if (strcmp(argv[1],   "--to_7") != 0) goto h;
    
    png_image_begin_read_from_file(&img, argv[2]); if (img.warning_or_error > 1) ret 1;
    
    if (img.format & PNG_FORMAT_FLAG_LINEAR) ret 1;
    img.format = (img.format & PNG_FORMAT_FLAG_ALPHA) ? PNG_FORMAT_RGBA : PNG_FORMAT_RGB;
    
    pm = (xpng_t){
        .w = img.width,
        .h = img.height,
        .A = img.format == PNG_FORMAT_RGBA,
    };
    
    MALLOC(pm.p, pm.s = pm.w * pm.h * (3 + pm.A)) ret 1;
    
    png_image_finish_read(&img, nil, pm.p, 0, nil); if (img.warning_or_error > 1) ret 1;
    
    ret (int)(normalize_RGBA(&pm) || store_7(&pm, argv[3]));
    
d:  if (load_7(argv[2], &pm)) ret 1;
    
    img.width  = pm.w,
    img.height = pm.h,
    img.format = pm.A ? PNG_FORMAT_RGBA : PNG_FORMAT_RGB;
    
    png_image_write_to_file(&img, argv[3], 0, pm.p, 0, nil);
    
    ret (int)(img.warning_or_error > 1);
    
h:  pf("\n"

    "./seven --to_7   example.png example.7\n"
    "./seven --to_png example.7   example.png\n"

    "\n"); ret 1;
}
