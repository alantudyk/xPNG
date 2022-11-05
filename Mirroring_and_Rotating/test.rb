#!/bin/ruby

system('./build.sh && cd .. && ./build.sh && cd 7 && ./build.sh && ' + 
       'cd ../Mirroring_and_Rotating') || exit(1)

Z = -> do
    system('../xpng -2 /tmp/x.7 /tmp/r.xpng > /dev/null') || raise
    File.size '/tmp/r.xpng'
end

f = $*[0] || raise
system("../7/seven --to_7 #{f} /tmp/s.7") || raise
puts

a = []; w = 0

['   0', '  90', ' 180', ' 270'].each do | r |

    ['    ', ' + v', ' + h'].each do | m |
    
        system(m[1] == ?+ ?
                   "./tool --m#{m[3]} /tmp/s.7 /tmp/x.7" :
                   'cp /tmp/s.7 /tmp/x.7') || raise
        a << [z = Z.call, "\t#{r}#{m}    =>    "]
        w = z if (z = z.to_s.size) > w
        
        break if [' 180', ' 270'].include? r
        
    end
    
    system('./tool --r90 /tmp/s.7 /tmp/s.7') || raise
    
end

a .sort .each { puts _1[1] + ("%#{w}d" % [_1[0]]) }

puts
