#! /usr/bin/ruby
#
# extractLength.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

Version="1.0.0"
banner = "Usage: extractLength.rb [option] <input SAM file>\n+Extract ID and length of sequence\n"

out_file=""

opt = OptionParser.new(banner)
opt.on("-o file", "output file (default: STDOUT)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)

MaxLine=1000

file = ARGV[0]
sam=SAMReader.new
sam.open(file)

sam.read_head

if(out_file == "")
  output = STDOUT
else
  output = File.open(out_file, "w")
end

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)

  for d in sam.data
    next if( d.flag & 256 == 256 || d.flag & 2048 == 2048)
    output.printf("%s\t%d\n", d.qname, d.seq.size)
  end

  sam.clear
end

output.close
