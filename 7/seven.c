#include "seven.h"
#include <png.h>

MAIN_ARGS {
    
    if (argc != 4) goto h; xpng_t pm;
    
    if (strcmp(argv[1], "--to_png") == 0) goto d;
    if (strcmp(argv[1],   "--to_7") != 0) goto h;
    
    
    
    ret (int)store_7(&pm, argv[3]);
    
d:  if (load_7(argv[2], &pm)) ret 1;
    
    
    
    ret 0;
    
h:  pf("\n"

    "./seven --to_7   example.png example.7\n"
    "./seven --to_png example.7   example.png\n"

    "\n"); ret 1;
}
