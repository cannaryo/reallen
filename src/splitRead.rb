#! /usr/bin/ruby
#
# splitRead.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require_relative "FastqReader.rb"
require "optparse"

Version="1.0.5"
banner = "Usage: splitRead.rb [option] <input SAM file>\n+Split reads by fixed length\n"

min_len = 100
sz = 35
out_file="tmp.fa"
is_fasta=false
no_split=false

opt = OptionParser.new(banner)
opt.on("-l min_len","minimum read length (default: 100)") {|v| min_len=v.to_i }
opt.on("-s size","split read size (default: 35)") {|v| sz=v.to_i }
opt.on("-o file","output file (default: tmp.fa)") {|v| out_file=v}
opt.on("--fasta","output in fasta format (default: false)") {|v| is_fasta=true}
opt.on("--no-split", "Do not split, just convert (default: false)") {|v| no_split=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nsplitRead.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
sam.read_head
fastq=FastqReader.new

MaxLine=1000
c_d,c_n=0,0

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)

  for d in sam.data
    len=d.seq.size
    next if(d.seq.size < min_len)
    id=d.qname
    id = "@" + id if(id[0] != "@")
    l_seq=d.seq[0,sz]
    l_qual=d.qual[0,sz]
    r_seq=d.seq[len-sz, sz]
    r_qual=d.qual[len-sz, sz]
    if(no_split)
      fastq.add_data(FastqField.new(id + ":" + d.rname, d.seq, "+", d.qual))
    else
      fastq.add_data(FastqField.new(id + ":LF",l_seq,"+",l_qual))
      fastq.add_data(FastqField.new(id + ":RF",r_seq,"+",r_qual))
    end
  end

  c_d += sam.size
  c_n = fastq.size / 2
  STDERR.printf("splitRead > Write records: %d pair / %d\r", c_n, c_d)
  sam.clear
end

output=File.open(out_file,"w")

if(is_fasta)
  for i in fastq.data
    output.print(">",i.id,"\n",i.seq,"\n") 
  end
else
  fastq.output(output)
end

STDERR.print("\n")
printf("%d sequences (%d pair) from %d data were witten in %s\n", fastq.data.size, fastq.data.size/2, c_d, out_file)

output.close
fastq.close
sam.close
