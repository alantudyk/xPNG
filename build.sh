#!/bin/sh

cd until_fork &&
./build.rb --gcc &&
cd .. &&
gcc -flto -O3 -s xpng.c -o xpng -D T_MAX=2 -Iuntil_fork \
-lpthread -Luntil_fork -lguntil_fork
