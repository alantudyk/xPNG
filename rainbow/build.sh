#!/bin/sh

gcc -O3 -s -flto ../until_fork/*.c rainbow.c -o rainbow
