#! /usr/bin/ruby

require "SAMReader.rb"
require "optparse"

Version="1.0.2"
banner = "Usage: sam2bed.rb [option] <input SAM file>\n"

t_sc=0
out_file=""
out_csv=false

opt = OptionParser.new(banner)
opt.on("-s score","score value. (1:MAPQ, 2:RNEXT+PNEXT, 3:TLEN, 4:AS/LEN)") {|v| t_sc=v.to_i }
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}
opt.on("--comma", "output by CSV file") {|v| out_csv=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)

file = ARGV[0]
sam=SAMReader.new
sam.open(file)

sam.read_all

if(out_file=="")
  output=STDOUT
else
  output=File.open(out_file,"w")
end

if(out_csv)
  output.printf("chr,start,end,name,score\n")
else
  output.printf("track name=\"%s\" description=\"Converted by sam2bed.rb\" type=bedDetail\n", file)
end

for d in sam.data
  case t_sc
  when 1
    score=d.mapq
  when 2
    score=d.rnext+":"+d.pnext.to_s
  when 3
    score=d.tlen
  when 4
    score=SAMDecoder.get_option(d,"AS").to_s + "/" + d.seq.size.to_s
  else
    score="."
  end
  len=SAMDecoder.get_aligned_ref_length(d)
  if(out_csv)
    output.printf("%s,%d,%d,%s,%s\n", d.rname, d.pos-1, d.pos+len-1, d.qname, score)
  else
    output.printf("%s\t%d\t%d\t%s\t%s\n", d.rname, d.pos-1, d.pos+len-1, d.qname, score)
  end
end

output.close
sam.close
