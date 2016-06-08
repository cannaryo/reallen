#! /usr/bin/ruby
#
# filtMapped.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

Version="1.1.7"
banner = "Usage: filtMapped.rb [option] <input SAM file>\n+Filter aligned pair\n"

out_file="tmp.sam"
m_len=-1
f_rp=false

opt = OptionParser.new(banner)
opt.on("-l len","keep pair distance  > len (default: none)") {|v| m_len=v.to_i}
opt.on("-r", "Write reference position") {|v| f_rp=true}
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nfiltMapped.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

MaxLine=1000
sam.read_head
out.write_head(sam.head)

ids = Hash.new
rec = Array.new
c_d,c_n=0,0

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  
  for d in sam.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    next if( k==nil || (d.flag & 256 == 256))
    if(ids.key?(k))
      if(d.rname != "*" && ids[k].rname != "*")
        if((d.rname != ids[k].rname) || (d.pos-ids[k].pos).abs > m_len)
          if(f_rp)
            d.rnext = ids[k].rname
            ids[k].rnext = d.rname
            d.rnext = ids[k].rnext = "=" if(d.rname == ids[k].rname)
            ids[k].pnext = d.pos
            d.pnext = ids[k].pos
          end
          rec.push(d)
          rec.push(ids[k])
        end
      end
      ids.delete(k)      
    else
      ids[k] = d
    end
  end

  c_d += sam.size
  c_n += rec.size
  STDERR.printf("filtMapped > Write records: %d / %d\r", c_n, c_d)
  out.write_data(rec)
  rec.clear
  sam.clear  
end

STDERR.print("\n")
printf("%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
