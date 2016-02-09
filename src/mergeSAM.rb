#! /usr/bin/ruby
#
# mergeSAM.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

Version="1.0.2"
banner = "Usage: mergeSAM.rb [option] <primary SAM file> <secondery SAM file>\n+Merge two SAM files\n"

out_file="tmp.sam"
keep_id=false

opt = OptionParser.new(banner)
opt.on("-i","keep reads with same ID") {|v| leep_id=true}
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 2)
printf("\nmergeSAM.rb %s\n",Version)

file_1 = ARGV[0]
sam_1=SAMReader.new
sam_1.open(file_1)
file_2 = ARGV[1]
sam_2=SAMReader.new
sam_2.open(file_2)

out=SAMWriter.new(out_file)

MaxLine=1000
sam_1.read_head
sam_2.read_head
out.write_head(sam_1.head)

ids = Hash.new
rec = Array.new
c_d,c_n=0,0
printf("Wrinting from %s\n", file_1)
while true
  STDOUT.flush
  MaxLine.times { sam_1.read_record }
  
  break if(sam_1.size==0)  
  for d in sam_1.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    rec.push(d)
    ids[k]=true
  end
  c_d += sam_1.size
  c_n += rec.size
  printf("Write records: %d / %d\r", c_n, c_d)
  out.write_data(rec)
  rec.clear
  sam_1.clear  
end

printf("\nWrinting from %s\n", file_2)
c2_d,c2_n=0,0
while true
  STDOUT.flush
  MaxLine.times { sam_2.read_record }
  
  break if(sam_2.size==0)  
  for d in sam_2.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    rec.push(d) unless(ids.key?(k) || keep_id)
  end
  c2_d += sam_2.size
  c2_n += rec.size
  printf("Write records: %d / %d\r", c2_n, c2_d)
  out.write_data(rec)
  rec.clear
  sam_2.clear  
end

printf("\n%d records in %d were witten in %s\n", c_n+c2_n, c_d+c2_d, out_file)

sam_1.close
sam_2.close
out.close
