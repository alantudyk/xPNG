#!/bin/ruby

system('./build.sh && cd ../7 && ./build.sh && cd ../Mirroring_and_Rotating') || exit(1)

['pigz-logo', 'olaf'].each do | f |
    
    f = '../images/' + f + '.png'
    
    system("../7/seven --to_7 #{f} /tmp/s.7") || raise
    
    system('./tool --r90  /tmp/s.7 /tmp/x.7 && ' +
           './tool --r270 /tmp/x.7 /tmp/r.7 && cmp /tmp/s.7 /tmp/r.7') || raise
    
    ['mv', 'mh', 'mvh', 'tl', 'tr'].each do
        system("./tool --#{_1} /tmp/s.7 /tmp/x.7 && " +
               "./tool --#{_1} /tmp/x.7 /tmp/r.7 && cmp /tmp/s.7 /tmp/r.7") || raise
    end
    
end
