#!/bin/sh

gcc -flto -O3 -s ../until_fork/*.c *.c -o seven -lpng
