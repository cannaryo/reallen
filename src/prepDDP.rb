#! /usr/bin/ruby
#
# prepDDP.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

def calc_start(dat)
  k,tag=dat.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
  if(tag =~ /L[MSF]$/)
    if(dat.flag & 16 == 16)
      pos = dat.pos + SAMDecoder.get_aligned_ref_length(dat) -1
    else
      pos = dat.pos
    end
  elsif(tag =~ /R[MSF]$/)
    if(dat.flag & 16 == 16)
      pos = dat.pos
    else
      pos = dat.pos + SAMDecoder.get_aligned_ref_length(dat) -1
    end
  end
  return pos
end

def calc_dir(dat)
  k,tag=dat.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
  dir="*"
  if(tag =~ /L[MSF]$/)
    if(dat.flag & 16 == 16)
      dir = "-"
    else
      dir = "+"
    end
  elsif(tag =~ /R[MSF]$/)
    if(dat.flag & 16 == 16)
      dir = "+"
    else
      dir = "-"
    end
  end
  return dir
end

def calc_fr(dat)
  return "r" if(dat.flag & 16 == 16)
  return "f"
end

def write_data(rec, out)
  for i in rec
    out.print(i, "\n")
  end
end


Version="1.0.5"
banner = "Usage: prepDDP.rb [option] <input SAM file>\n+prepare for DDP alignment\n"

out_file="tmp.pos"
ref_index_file=nil
ref_len=200

opt = OptionParser.new(banner)
opt.on("-o out_pos", "output position file (default: tmp.pos)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nprepDDP.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

sam=SAMReader.new
sam.open(ARGV[0])
out=File.open(out_file,"w")

MaxLine=1000
sam.read_head

ids = Hash.new
rec = Array.new
c_d,c_n=0,0


while true
  STDOUT.flush
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  
  for d in sam.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    next if( k==nil || (d.flag & 256 == 256))
    if(ids.key?(k))
      if(d.rname != "*" && ids[k].rname != "*")
        if(tag =~ /L[MSF]$/)
          t1, t2 = d, ids[k]
        else
          t1, t2 = ids[k], d
        end
        rec_data=sprintf("%s,%s,%d,%s,%s,%s,%d,%s,%s",k, t1.rname, calc_start(t1), calc_fr(t1), calc_dir(t1), 
                         t2.rname, calc_start(t2), calc_fr(t2), calc_dir(t2))
        rec.push(rec_data)
      end
      ids.delete(k)
    else
      ids[k] = d
    end
  end

  c_d += sam.size
  c_n += rec.size
  printf("Write records: %d / %d\r", c_n, c_d)
  write_data(rec, out)
  rec.clear
  sam.clear  
end

printf("\n%d records in %d (%d pair) were witten in %s\n", c_n, c_d, c_d/2, out_file)

sam.close
out.close
