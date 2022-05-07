#include "7/seven.h"

MAIN_ARGS {
    
    if (argc != 4 || argv[1][0] != '-' || strlen(argv[1]) != 2) goto h; xpng_t pm;
    
    switch (argv[1][1]) {
    
        case '1':
        case '2':
        case '7': ret (int)(load_7(argv[2], &pm) || xpng_store(argv[1][1] - '0', &pm, argv[3]));
        case '3': ret (int)xpng_from_jpg(argv[2], argv[3]);
        case 'd': ret (int)(xpng_load(argv[2], &pm) || store_7(&pm, argv[3]));
        
    }
    
h:  pf("\n"

    "encode: ./xpng -[127] example.7    example.xpng\n"
    "        ./xpng -3     example.jpg  example.xpng\n"
    "decode: ./xpng -d     example.xpng example.7\n"

    "\n"); ret 1;
}
