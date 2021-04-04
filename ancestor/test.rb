#!/bin/ruby

system('gcc -flto -O3 -s gray.c -o gray -D T_MAX=2 -lpthread') || raise

['-1', '-2'].each do | m |
    p = $*[0] || ?.
    p += ?/ if p[-1] != ?/
    sum = 0
    Dir.entries(p).each do | f |
        next unless f[/\.png$/]
        system("convert #{p + f} /tmp/src.ppm") || (puts 'convert ' + f; raise)
        system("./gray #{m} /tmp/src.ppm /tmp/res.gray > /dev/null") ||
            (puts 'g_c ' + f; raise)
        system('./gray -d /tmp/res.gray /tmp/res.ppm > /dev/null') ||
            (puts 'g_d ' + f; raise)
        sum += (z = File.size '/tmp/res.gray')
        (puts 'g_cmp ' + f; raise) unless system 'cmp /tmp/src.ppm /tmp/res.ppm'
        # puts '%11d %s' % [z, f]
    end
    puts '%11d %s' % [sum, m]
end
