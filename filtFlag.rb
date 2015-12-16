#! /usr/bin/ruby

require "SAMReader.rb"
require "optparse"

Version="1.0.4"
banner = "Usage: filtFlag.rb [option] <flag array> <input SAM file>\n+Filter reads by flag\n"

out_file="tmp.sam"
show_help=false

opt = OptionParser.new(banner)
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}
opt.on("--table", "show flag table and exit") {|v| show_help=true}

opt.parse!(ARGV)

if(ARGV.size < 2 || show_help)
  print("Usage: filtByFlag.rb [option] <flag array> <input SAM file>\n\n")
  print("Flag array should be argued 12 degit format\n")
  print("0: pass bit 0, 1: pass bit 1, .: pass both\n")
  print("\nflag table is as follows:\n")
  print("0x1 template having multiple segments in sequencing\n")
  print("0x2 each segment properly aligned according to the aligner\n")
  print("0x4 segment unmapped\n")
  print("0x8 next segment in the template unmapped\n")
  print("0x10 SEQ being reverse complemented\n")
  print("0x20 SEQ of the next segment in the template being reversed\n")
  print("0x40 the first segment in the template\n")
  print("0x80 the last segment in the template\n")
  print("0x100 secondary alignment\n")
  print("0x200 not passing quality controls\n")
  print("0x400 PCR or optical duplicate\n")
  print("0x800 supplementary alignment\n\n")
  exit
end

MaxLine=1000

flags=ARGV[0]
file = ARGV[1]

if(flags.size!=12)
  print("flag array must be specified 12 degit format\n")
  exit
end

sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

sam.read_head
out.write_head(sam.head)
rec = Array.new

bit_1,bit_2=0,0
c_n,c_d=0,0

for i in 0...12
  if(flags[i,1] == "1")
    bit_1 += 2**i
  elsif(flags[i,1] != "0")
    bit_2 += 2**i
  end
end

while true
  STDOUT.flush
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  for d in sam.data
    if((4095^(d.flag^bit_1))|bit_2 == 4095)
      rec.push(d)
    end
  end
  c_d += sam.size
  c_n += rec.size
  printf("Write records: %d / %d\r", c_n, c_d)
  out.write_data(rec)
  rec.clear
  sam.clear
end

printf("\n%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
