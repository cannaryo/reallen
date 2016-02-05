#! /usr/bin/ruby

require "SAMReader.rb"
require "FastqReader.rb"
require "optparse"

Version="1.3.1"
banner = "Usage: splitSoftClip.rb [option] <input SAM file>\n+Split soft-clipped sequences\n"

min_len = 25
min_rlen = 25
out_file=""
ex_file=""
is_fasta=false

opt = OptionParser.new(banner)
opt.on("-l len","minimum length clipped read (default: 25)") {|v| min_len=v.to_i }
opt.on("-r len","minimum length remaining read (default: 25)") {|v| min_rlen=v.to_i }
opt.on("-e file","extract filtered data (default: None)") {|v| ex_file=v}
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}
opt.on("--fasta","output in fasta format (default: false)") {|v| is_fasta=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nsplitSoftClip.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)

sam.read_all

fastq=FastqReader.new

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

if(out_file=="")
  output=STDOUT
else
  output=File.open(out_file,"w")
end

if(is_fasta)
  for i in fastq.data
    output.print(">",i.id,"\n",i.seq,"\n") 
  end
else
  fastq.output(output)
end

if(ex_file != "")
  out_sam=SAMWriter.new(ex_file)
  out_sam.write_head(sam.head)
  data=Array.new
  for d in sam.data
    next if d.qname == "*"    
    data.push(d)
  end
  out_sam.write_data(data)
  out_sam.close
end

printf("%d sequences (%d pair) from %d data were witten in %s\n", fastq.data.size, fastq.data.size/2, sam.data.size, out_file)

output.close
fastq.close
sam.close
