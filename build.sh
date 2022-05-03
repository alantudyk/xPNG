#!/bin/sh

gcc -flto -O3 -s until_fork/*.c 7/libseven.c libxpng.c xpng.c -o xpng \
-D T_MAX=`getconf _NPROCESSORS_ONLN` -D SYNC_IO -lpthread
