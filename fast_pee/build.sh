#!/bin/sh

gcc -O3 -s -flto ../until_fork/*.c fast_pee.c -o fast_pee
