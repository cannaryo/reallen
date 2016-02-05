#! /usr/bin/ruby

fp=File.open(ARGV[0])

data_all=Array.new
fp.each_line do |s|
  d=s.split(" ")
  next if d.size < 6
  data_all.push([d[0], d[1].to_i, d[2].to_i, d[3], d[4],d[5].to_f/(d[2].to_i - d[1].to_i)])
end

data_all.sort!{ |a,b| a[0] == b[0] ? a[2] <=> b[2] : a[0] <=> b[0] }

printf("Chromosome,Start,End,ID,Score,Coverage\n")

data_all.each do |d|
  printf("%s,%d,%d,%s,%s,%.2f\n", d[0],d[1],d[2],d[3],d[4],d[5])
end
