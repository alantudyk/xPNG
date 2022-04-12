#include "7/seven.h"

MAIN_ARGS {
    
    if (argc != 4 || argv[1][0] != '-' || strlen(argv[1]) != 2) goto h;
    xpng_t pm; reg_t r;
    
    switch (argv[1][1]) {
    
        case '1':
        case '2':
        case '7':
            
            if (f_read(argv[2], &r) ||
                sscanf((char *)(r.p + 3), "%lu %lu", &pm.w, &pm.h) != 2) ret 1;
            pm.s = 3 * pm.w * pm.h, pm.p = r.p + (r.s - pm.s), pm.A = 0;
            ret (int)xpng_from_pixmap(argv[1][1] == '7' ? 0 : argv[1][1] - '0', &pm, argv[3]);
        
        case '3': ret (int)xpng_from_jpg(argv[2], argv[3]);
        case 'd':
            
            if (xpng_to_pixmap(argv[2], &pm)) ret 1;
            char P6[100]; int l = sprintf(P6, "P6\n%lu %lu\n255\n", pm.w, pm.h);
            FILE *of = fopen(argv[3], "wb"); if (of == NULL) ret 1;
            ret (int)(fwrite(P6, 1, l, of) != l || fwrite(pm.p, 1, pm.s, of) != pm.s || fclose(of));
        
    }
    
h:  pf("\n"

    "encode: ./xpng -[127] example.ppm  example.xpng\n"
    "        ./xpng -3     example.jpg  example.xpng\n"
    "decode: ./xpng -d     example.xpng example.ppm\n"

    "\n"); ret 1;
}
