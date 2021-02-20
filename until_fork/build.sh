#!/bin/sh

gcc -c -flto -O3 until_fork.c &&
gcc-ar cr libguntil_fork.a until_fork.o
