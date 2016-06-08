#! /usr/bin/ruby
#
# splitSoftClip.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require_relative "FastqReader.rb"
require "optparse"

Version="1.3.2"
banner = "Usage: splitSoftClip.rb [option] <input SAM file>\n+Split soft-clipped sequences\n"

min_len = 25
min_rlen = 25
out_file="tmp.fa"
ex_file=""
is_fasta=false

opt = OptionParser.new(banner)
opt.on("-l len","minimum length clipped read (default: 25)") {|v| min_len=v.to_i }
opt.on("-r len","minimum length remaining read (default: 25)") {|v| min_rlen=v.to_i }
opt.on("-e file","extract filtered data (default: None)") {|v| ex_file=v}
opt.on("-o file","output file (default: tmp.fa)") {|v| out_file=v}
opt.on("--fasta","output in fasta format (default: false)") {|v| is_fasta=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nsplitSoftClip.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
sam.read_head
fastq=FastqReader.new

if(ex_file != "")
  out_sam=SAMWriter.new(ex_file)
  out_sam.write_head(sam.head)
  rec=Array.new
end

MaxLine=1000
c_d,c_n=0,0

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  
  for d in sam.data
    len=d.seq.size
    comm=d.cigar.scan(/(\d+)([MIDNSHP=X])/)
    qual=SAMDecoder.get_quality(d)
    l_tag="LM"
    r_tag="RM"
    l_seq=l_qual=r_seq=r_qual=""
    sz_l=sz_r=0
    if(comm.first != nil && comm.first[1] == "S" )
      sz_l=comm.first[0].to_i
    end
    if(comm.last != nil && comm.last[1] == "S" )
      sz_r=comm.last[0].to_i
    end
    if(sz_l > sz_r)
      next if(sz_l < min_len || (len-sz_l) < min_rlen)
      sz=sz_l
      l_seq=d.seq[0,sz]
      l_qual=d.qual[0,sz]
      r_seq=d.seq[sz, len-sz]
      r_qual=d.qual[sz, len-sz]
      l_tag="LS"
    else
      next if(sz_r < min_len || (len-sz_r) < min_rlen)
      sz=sz_r
      l_seq=d.seq[0, len-sz]
      l_qual=d.qual[0,len-sz]
      r_seq=d.seq[len-sz, sz]
      r_qual=d.qual[len-sz, sz]
      r_tag="RS"
    end
    id=d.qname
    id = "@" + id if(id[0] != "@")
    d.qname="*"
    fastq.add_data(FastqField.new(id + ":" + l_tag, l_seq,"+",l_qual))
    fastq.add_data(FastqField.new(id + ":" + r_tag, r_seq,"+",r_qual))
  end

  if(ex_file != "")
    for d in sam.data
      next if d.qname == "*"    
      rec.push(d)
    end
    out_sam.write_data(rec)
    rec.clear
  end

  c_d += sam.size
  c_n = fastq.size / 2
  STDERR.printf("splitSoftClip > Write records: %d pair / %d\r", c_n, c_d)
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

if(ex_file != "")
  out_sam.close
  printf("%d recoeds were witten in %s\n", c_d - c_n, ex_file)
end

output.close
fastq.close
sam.close
