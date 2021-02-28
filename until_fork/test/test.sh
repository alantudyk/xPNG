#!/bin/sh

cd .. && ./build.rb --all && cd test &&

gcc   -flto -O3 -s rsort_test.c -o rsort_gcc   -lpthread -L.. -lguntil_fork &&
clang -flto -O3 -s rsort_test.c -o rsort_clang -lpthread -L.. -lcuntil_fork &&

echo '\nTest rsort_gcc:\n'   && ./rsort_gcc &&
echo   'Test rsort_clang:\n' && ./rsort_clang
