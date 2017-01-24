#! /usr/bin/ruby
#
# sampleRead.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

Version="1.0.1"
banner = "Usage: sampleRead.rb [option] <input SAM file>\n+Sample read for test\n"

sample_rate = 0.5
rnd = Random.new
out_file="tmp.sam"

opt = OptionParser.new(banner)
opt.on("-p rate","Sampling Rate (default: 0.5)") { |v| sample_rate = v.to_f }
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nsampleRead.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

MaxLine=1000

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

sam.read_head
out.write_head(sam.head)

c_n,c_d=0,0

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  dat = Array.new
  for d in sam.data
    if( rnd.rand < sample_rate)
      dat.push(d)
    end
  end
  c_n += dat.size
  c_d += sam.size
  out.write_data(dat)
  STDERR.printf("sampleRead > Write records: %d / %d\r", c_n, c_d)
  sam.clear
end
STDERR.print("\n")
printf("%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
