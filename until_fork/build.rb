#!/bin/ruby

CMD = '%s -c -flto -O3 *.c && ' +
      '%s-ar cr lib%suntil_fork.a *.o'

FMT = [['gcc', 'gcc', 'g'], ['clang', 'llvm', 'c']]

case $*[0]
    when '--gcc'  ; system(CMD % FMT[0]) || raise
    when '--clang'; system(CMD % FMT[1]) || raise
    when '--all'  ; FMT.each { system(CMD % _1) || raise }
    else            raise
end
