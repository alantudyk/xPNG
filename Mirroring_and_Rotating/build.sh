#!/bin/sh

gcc -O3 -s -flto ../until_fork/*.c ../7/libseven.c tool.c -o tool
