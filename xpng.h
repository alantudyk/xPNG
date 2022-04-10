#ifndef ___XPNG_H
#define ___XPNG_H

#include <until_fork.h>

#define XPNG_COMPRESSION_TYPE_FAST 1
#define XPNG_COMPRESSION_TYPE_SLOW 2
#define XPNG_COMPRESSION_TYPE_EXJPEG 3
#define XPNG_COMPRESSION_TYPE_UNCOMPRESSED 7

typedef struct xpng_t { u8_t *p; u64_t w, h, s; _Bool A; } xpng_t;

_Bool xpng_from_pixmap(const u64_t mode, const xpng_t *pm, const char *const xpng);
_Bool xpng_to_pixmap(const char *const xpng, xpng_t *pm);

_Bool xpng_from_jpg(const char *const jpg, const char *const xpng);

_Bool xpng_from_pixmap_T(u64_t T, const u64_t mode, const xpng_t *pm, const char *const xpng);
_Bool xpng_to_pixmap_T(u64_t T, const char *const xpng, xpng_t *pm);

_Bool xpng_from_jpg_T(u64_t T, const char *const jpg, const char *const xpng);

#endif
