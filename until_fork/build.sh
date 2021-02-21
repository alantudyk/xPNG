#!/bin/sh

gcc -c -flto -O3 until_fork.c &&
gcc-ar cr libguntil_fork.a until_fork.o

# clang -c -flto -O3 until_fork.c &&
# llvm-ar cr libcuntil_fork.a until_fork.o
