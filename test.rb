#!/bin/ruby

system('./build.sh && cd 7 && ./build.sh && cd ..') || exit(1)

OK     = -> do ["\x1b[32mOK\x1b[m", File.size('/tmp/res.xpng')] end
Failed = -> do ["\x1b[31mFailed\x1b[m", 'n/a'] end

Print_Result = -> o, success do

    fmt =  "\t\toption '-#{o}': %28s\n\t\tcompressed size: %16s\n\n"
    puts fmt % (success ? OK : Failed).call
    
end

p = $*[0] || 'images'
p += ?/ if p[-1] != ?/

puts

Dir.entries(p).sort.each do | f |
    
    f[/\.(pn|jp)g$/] ? puts("\n\tFile '#{f}':\n\n") : next
    
    if f[-2] == ?n
        
        system("7/seven --to_7 #{p + f} /tmp/src.7") ||
            (puts "convertion to *.7 failed\n\n"; next)
        
        (1..2).each do | o |
            
            Print_Result.call o, system(<<~CMD)
                ./xpng -#{o} /tmp/src.7 /tmp/res.xpng > /dev/null &&
                ./xpng -d /tmp/res.xpng /tmp/res.7 > /dev/null &&
                cmp /tmp/src.7 /tmp/res.7 > /dev/null
            CMD
            
        end
        
    else
        
        Print_Result.call 3, system("./xpng -3 #{p + f} /tmp/res.xpng > /dev/null")
        
    end
    
end

puts
