#!/bin/sh

./build.rb --all &&
sudo cp until_fork.h /usr/include/until_fork.h &&
sudo cp libguntil_fork.a /usr/lib/libguntil_fork.a &&
sudo cp libcuntil_fork.a /usr/lib/libcuntil_fork.a
