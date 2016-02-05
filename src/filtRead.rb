#! /usr/bin/ruby
# filtRead.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

def count_cigar(cigar)
  t=cigar.scan(/(\d+)([MIDNSHP=X])/) 
  cc_hash={"M"=>0, "I"=>0, "D"=>0, "N"=>0, "S"=>0, "H"=>0, "P"=>0, "="=>0, "X"=>0}
  for i in t
    cc_hash[i[1]] += i[0].to_i
  end
  return cc_hash
end

def filter_data(dat, m_d, m_i, m_id, m_s, m_m, umap, rev)
  out = Array.new
  for i in dat
    cc_hash=count_cigar(i.cigar)
    if(umap && i.rname == "*")
      out.push(i)
    elsif(m_s != nil && ( (!rev && m_s <= cc_hash["S"]) || (rev && m_s >= cc_hash["S"])))
      out.push(i)
    elsif(m_d != nil && ( (!rev && m_d <= cc_hash["D"]) || (rev && m_d >= cc_hash["D"])))
      out.push(i)
    elsif(m_i != nil && ( (!rev && m_i <= cc_hash["I"]) || (rev && m_i >= cc_hash["I"])))
      out.push(i)
    elsif(m_m != nil && ( (!rev && m_m <= cc_hash["M"]) || (rev && m_m >= cc_hash["M"])))
      out.push(i)
    elsif(m_id != nil && m_id <= (cc_hash["I"] + cc_hash["D"]) )
      out.push(i)
    end
  end
  return out
end

Version="1.7.1"
banner = "Usage: filtRead.rb [option] <input SAM file>\n+Filter reads by cigar\n"

min_del = nil
min_ins = nil
min_indel = nil
min_softclip = nil
min_match = nil
f_rev = false
unmapped = false
out_file="tmp.sam"

opt = OptionParser.new(banner)
opt.on("-d N","Keep deletion >= N") {|v| min_del=v.to_i }
opt.on("-i N","Keep insertion >= N") {|v| min_ins=v.to_i }
opt.on("--indel N","Keep indel >= N") {|v| min_indel=v.to_i }
opt.on("-s N","Keep softclip >= N") {|v| min_softclip=v.to_i }
opt.on("-m N","Keep match >= N") {|v| min_match=v.to_i }
opt.on("-r","Reverse inequality sign for [dism]") {|v| f_rev=true}
opt.on("-u","Keep unmapped read") { |v| unmapped = true }
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nfiltRead.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

MaxLine=1000

file = ARGV[0]
sam=SAMReader.new
sam.open(file)
out=SAMWriter.new(out_file)

sam.read_head
out.write_head(sam.head)

c_n,c_d=0,0

while true
  STDOUT.flush
  MaxLine.times { sam.read_record }
  break if(sam.size==0)

  dat=filter_data(sam.data, min_del, min_ins, min_indel, min_softclip, min_match, unmapped, f_rev )
  c_n += dat.size
  c_d += sam.size
  out.write_data(dat)
  printf("Write records: %d / %d\r", c_n, c_d)
  sam.clear
end
printf("\n%d records in %d were witten in %s\n", c_n, c_d, out_file)

sam.close
out.close
