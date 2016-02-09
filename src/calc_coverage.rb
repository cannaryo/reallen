#! /usr/bin/ruby
#
# calc_coverage.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require "optparse"

Version="1.1.0"
banner = "Usage: calcCoverage.rb [option] <input BED file> <original BAM file>\n+calculate coverage at each position\n"

out_file="tmp.cov.csv"
tmp_file="cov_tmp.bed"
sam_comm="samtools"

opt = OptionParser.new(banner)
opt.on("-o file","output file (default: tmp.bed)") {|v| out_file=v}
opt.on("--tmp file", "temporary file (default: cov_tmp.bed)") {|v| tmp_file=v}
opt.on("--samtools command", "command for Samtools (default: samtools)") {|v| sam_comm=v}

opt.parse!(ARGV)

exit if(ARGV.size < 2)
printf("\ncalcCoverage.rb %s\n",Version)
printf("input=%s, BAM=%s\n", ARGV[0],ARGV[1])

if(File.exist?(tmp_file))
  printf("Cannot write %s\n",tmp_file)
  printf("File exist\n")
  exit
end

str=sam_comm+" bedcov "+ARGV[0]+" "+ARGV[1]+" > "+tmp_file
system(str)

data_all=Array.new
File.open(tmp_file).each_line do |s|
  d=s.split(" ")
  next if d.size < 6
  data_all.push([d[0], d[1].to_i, d[2].to_i, d[3], d[4],d[5].to_f/(d[2].to_i - d[1].to_i)])
end

data_all.sort!{ |a,b| a[0] == b[0] ? a[2] <=> b[2] : a[0] <=> b[0] }

File.open(out_file,"w") do |fp|
  fp.printf("Chromosome,Start,End,ID,Score,Coverage\n")
  data_all.each do |d|
    fp.printf("%s,%d,%d,%s,%s,%.2f\n", d[0],d[1],d[2],d[3],d[4],d[5])
  end
end

File.delete(tmp_file)
