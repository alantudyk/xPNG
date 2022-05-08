#ifndef ___XPNG_H
#define ___XPNG_H

#include "until_fork/until_fork.h"

#define XPNG_COMPRESSION_TYPE_FAST 1
#define XPNG_COMPRESSION_TYPE_SLOW 2
#define XPNG_COMPRESSION_TYPE_EXJPEG 3
#define XPNG_COMPRESSION_TYPE_UNCOMPRESSED 7

typedef struct xpng_t { u8_t *p; u64_t w, h, s; _Bool A; } xpng_t;

CHECK _Bool xpng_store(u64_t mode, const xpng_t *pm, const char *xpng);
CHECK _Bool xpng_load(const char *xpng, xpng_t *pm);

CHECK _Bool xpng_from_jpg(const char *jpg, const char *xpng);

CHECK _Bool xpng_store_T(u64_t T, u64_t mode, const xpng_t *pm, const char *xpng);
CHECK _Bool xpng_load_T(u64_t T, const char *xpng, xpng_t *pm);

CHECK _Bool xpng_from_jpg_T(u64_t T, const char *jpg, const char *xpng);

#endif
