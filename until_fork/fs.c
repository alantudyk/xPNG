#include "until_fork.h"

_Bool f_size(const char *fn, u64_t *s) {
    
    struct stat st;
    if (stat(fn, &st) != 0) ret 1;
    *s = st.st_size;
    
    ret 0;
}

#define RET fclose(fp); ret
_Bool f_read(const char *fn, reg_t *f) {
    FILE* fp = fopen(fn, "rb"); if (fp == NULL) ret 1;
    F_SIZE(fn, &(f->s)) { RET 1; }
    MALLOC(f->p, f->s) { RET 1; }
    if (fread(f->p, 1, f->s, fp) != f->s) { free(f->p); RET 1; }
    RET 0;
}

_Bool f_write(const char *fn, const reg_t *f) {
    
    FILE* fp;
    
    if ((fp = fopen(fn, "wb")) == NULL) ret 1;
    if (fwrite(f->p, 1, f->s, fp) != f->s) { RET 1; }
    ret (_Bool)fclose(fp);
}
