#! /usr/bin/ruby
#
# pairTest.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

class MeasureDistance
  def initialize
    @use_file = false
    @hash = nil
  end

  def open(file)
    @hash = Hash.new
    File.open(file).each_line do |s|
      if(s =~  /^([^:]+(?::[^:]+){1,9})\s+(\d+)$/)
        @hash[$1] = $2.to_i
      end          
    end
    @use_file = true
  end

  def bp_position(dat)
    k,tag=dat.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    if(tag =~ /L[MSF]$/)
      if(dat.flag & 16 == 16)
        return dat.pos
      else
        return dat.pos + SAMDecoder.get_aligned_ref_length(dat) -1
      end
    elsif(tag =~ /R[MSF]$/)
      if(dat.flag & 16 == 16)
        return dat.pos + SAMDecoder.get_aligned_ref_length(dat) -1
      else
        return dat.pos
      end
    end
    return 0
  end

  def get_id(dat)
    k,tag=dat.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    return k
  end

  def check_pair(l_dat, r_dat, min_len)
    if(l_dat.rname != r_dat.rname)
      return false
    elsif((l_dat.flag & 16 == 16) != (r_dat.flag & 16 == 16))
      return false
    end        
    l_bp = bp_position(l_dat)
    r_bp = bp_position(r_dat)
    id = get_id(l_dat)
    if( @use_file && @hash.key?(id) )
      size = @hash[id] - (l_dat.seq.size + r_dat.seq.size)
      return false if((l_bp - r_bp).abs > min_len + size)      
    else
      return false if((l_bp - r_bp).abs > min_len)      
    end
    return true
  end

end

Version="1.2.0"
banner = "Usage: pairTest.rb [option] <input SAM file>\n+Filter aligned pair\n"

out_file="tmp.sam"
m_len=-1
f_ref_pos=false
len_file = ""

opt = OptionParser.new(banner)
opt.on("-l len", "keep pair distance  > len (default: none)") {|v| m_len=v.to_i}
opt.on("-r", "Write reference position") {|v| f_ref_pos=true}
opt.on("-o file", "output file (default: tmp.sam)") {|v| out_file=v}
opt.on("-f file", "use length file to measure correct distance") {|v| len_file = v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\npairTest.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

measure_dist = MeasureDistance.new
if(len_file != "")
  measure_dist.open(len_file)
end

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
        if( !measure_dist.check_pair(d, ids[k], m_len) )
          if(f_ref_pos)
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
  STDERR.printf("pairTest > Write records: %d / %d\r", c_n, c_d)
  out.write_data(rec)
  rec.clear
  sam.clear  
end

STDERR.print("\n")
printf("%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
