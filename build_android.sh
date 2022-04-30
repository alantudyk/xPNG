#!/bin/sh

/home/"$USER"/android_ndk/android-ndk-r23b/toolchains/llvm/prebuilt/\
linux-x86_64/bin/aarch64-linux-android31-clang \
-flto -O3 -s until_fork/*.c 7/libseven.c libxpng.c xpng.c -o xpng_android \
-D T_MAX=0 -D SYNC_IO
