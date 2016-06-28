#! /usr/bin/ruby
#
# bp2bed.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require "optparse"

Version="1.2.0"
banner = "Usage: bp2bed.rb [option] <input BP file>\nConvert BP file in BED format"

t_sc=0
len=30
out_file=""

opt = OptionParser.new(banner)
opt.on("-s score","score value. (1:Dir. & Cov. ; 2:Dist. Pos.)") {|v| t_sc=v.to_i }
opt.on("-l len","Length of region (default: 30)") {|v| len=v.to_i }
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)

file=ARGV[0]
input=File.open(file)

if(out_file=="")
  output=STDOUT
else
  output=File.open(out_file,"w")
end

id=rand(100000)

output.printf("track name=\"%s\" description=\"Converted by bp2bed.rb\" type=bedDetail\n", file)

while(l=input.gets)
  arr = l.split(",").map{|i| i.strip}
  return if arr.size != 15

  for i in [0,1]
    pp=arr[1+i*7,7]
    rr=arr[8-i*7,7]
    # arr[0]: ID
    # pp : [0] side, [1] chr, [2] pos, [3] dir, [4] ter, [5] cov, [6] len
    if(pp[4] == "p")
      p1 = pp[2].to_i - len
      p2 = pp[2].to_i
    else
      p1 = pp[2].to_i - 1
      p2 = pp[2].to_i + len - 1
    end
    score="."
    if(t_sc == 1)
      score = pp[0] + pp[3] + pp[4] + ":" + pp[5]
    elsif(t_sc == 2)
      score= rr[1] + ":" + rr[2] + ":" + rr[4]
    end
    output.printf("%s\t%d\t%d\t%s\t%s\n", pp[1], p1, p2, arr[0], score)
  end

end

input.close
output.close
