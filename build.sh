#!/bin/sh

cd until_fork && ./build.rb --gcc && cd .. &&

gcc -flto -O3 -s 7/libseven.c libxpng.c xpng.c -o xpng \
-D T_MAX=`getconf _NPROCESSORS_ONLN` \
-lpthread -Luntil_fork -lguntil_fork
