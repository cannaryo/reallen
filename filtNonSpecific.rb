#! /usr/bin/ruby

require "SAMReader.rb"
require "optparse"

Version="1.2.1"
banner = "Usage: filtNonSpecific.rb [option] <input SAM file>\n+Filter non-specific alignment\n"

out_file="tmp.sam"
min_score = 25
min_frac=0.0
f_sorted=true

opt = OptionParser.new(banner)
opt.on("-s score", "minimum difference with second alignment (default: 25)") {|v| min_score = v.to_i}
opt.on("-f fraction", "minimum alignment score fraction (default:0.0)") {|v| min_frac=v.to_f}
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}
opt.on("--not-sorted","The alignment is not sorted (slow)") {|v| f_sorted=false }

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nfiltNonSpecific.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

MaxLine=1000

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

sam.read_head
out.write_head(sam.head)
ids = Hash.new
id_tmp = ["DUMMY_ID",nil]
c_d=0
c_n=0

while true
  STDOUT.flush
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  for d in sam.data
    k=d.qname.scan(/([^:]+(?::[^:]+){1,9})$/).first.first
    next if k==nil
    if(id_tmp[0] == k)
      id_tmp[1] = nil if(id_tmp[1] != nil && (SAMDecoder.get_option(id_tmp[1],"AS") - SAMDecoder.get_option(d,"AS")).abs < min_score)
    elsif(!f_sorted && ids.key?(k))
      ids.delete(k) if((SAMDecoder.get_option(ids[k],"AS") - SAMDecoder.get_option(d,"AS")).abs < min_score)
    elsif(d.flag & 256 != 256)
      ids[id_tmp[0]] = id_tmp[1] if(id_tmp[1] != nil && SAMDecoder.get_option(id_tmp[1],"AS").to_f/id_tmp[1].seq.size >= min_frac)
      id_tmp = [k,d]
    end
  end
  c_d += sam.size  
  if(f_sorted)
    c_n+=ids.size
    printf("Write records: %d / %d\r", c_n, c_d)
    out.write_data(ids.values)
    ids.clear
  else
    printf("Check records: %d / %d\r", ids.size, c_d)
  end
  sam.clear  
end

ids[id_tmp[0]] = id_tmp[1] if(id_tmp[1] != nil && SAMDecoder.get_option(id_tmp[1],"AS").to_f/id_tmp[1].seq.size >= min_frac)
c_n+=ids.size

out.write_data(ids.values)
printf("\n%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
