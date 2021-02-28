#!/bin/sh

cd until_fork && ./build.rb --gcc && cd .. &&

gcc -flto -O3 -s xpng.c -o xpng \
-D T_MAX=`getconf _NPROCESSORS_ONLN` \
-Iuntil_fork -lpthread -Luntil_fork -lguntil_fork
