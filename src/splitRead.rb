#! /usr/bin/ruby

require "SAMReader.rb"
require "FastqReader.rb"
require "optparse"

Version="1.0.2"
banner = "Usage: splitRead.rb [option] <input SAM file>\n+Split reads by fixed length\n"

min_len = 100
sz = 35
out_file=""
is_fasta=false
no_split=false

opt = OptionParser.new(banner)
opt.on("-l min_len","minimum read length (default: 100)") {|v| min_len=v.to_i }
opt.on("-s size","split read size (default: 35)") {|v| sz=v.to_i }
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}
opt.on("--fasta","output in fasta format (default: false)") {|v| is_fasta=true}
opt.on("--no-split", "Do not split, just convert (default: false)") {|v| no_split=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nsplitRead.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)

sam.read_all

fastq=FastqReader.new

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

printf("%d sequences (%d pair) from %d data were witten in %s\n", fastq.data.size, fastq.data.size/2, sam.data.size, out_file)

output.close
fastq.close
sam.close
