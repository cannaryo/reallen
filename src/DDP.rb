#! /usr/bin/ruby
#
# DDP.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#


require_relative "SAMReader.rb"
require_relative "AlignDDP.rb"
require_relative "SeqTools.rb"
require "optparse"

def write_data(rec, out)
  for i in rec
    out.print(i, "\n")
  end
end


Version="1.1.5"
banner = "Usage: DDP.rb [option] <reference> <input SAM file> <position file>\n+DDP alignment\n"

out_file="tmp.sam"
ref_index_file=nil
ref_margin=20
ref_sparse=nil

opt = OptionParser.new(banner)
opt.on("-i ref_index", "Reference index file (to skip indexing)") {|v| ref_index_file=v}
opt.on("-o file","output file (default: tmp.sam)") {|v| out_file=v}
opt.on("--margin-length length", "margin of reference cut (default: 20)") {|v| ref_margin=v.to_f} 
opt.on("--reference-sparse N", "set sparse code each N charactor in reference") {|v| ref_sparse=v.to_f}
opt.parse!(ARGV)

exit if(ARGV.size < 3)
printf("\nDDP.rb %s\n",Version)
printf("reference=%s\ninput=%s\nposition=%s\n", ARGV[0], ARGV[1], ARGV[2])

ref_file=File.open(ARGV[0])
ref_index=Hash.new

if(ref_index_file==nil)
  printf("Preparing index...\n")
  while(d=ref_file.gets)
    if(d=~/^>(\w+)$/)
      k=$1
      ref_index[k]=ref_file.pos
      printf("%s: %d\n",k, ref_file.pos)
    end
  end
  ref_file.rewind
else
  File.open(ref_index_file).each_line do |l|
    if(l =~ /^(\w+)\t\d+\t(\d+)/)
      ref_index[$1] = $2.to_i
    end
  end
end

# test reference sparse
if(ref_sparse == nil)
  while(s=ref_file.gets)
    t_1 = ref_file.gets
    t_2 = ref_file.gets
    if(t_1.size == s.size && t_2.size == s.size)
      ref_sparse = s.size
      break
    end
  end
end

sam=SAMReader.new
sam.open(ARGV[1])

out=SAMWriter.new(out_file)
#out=File.open(out_file, "w")

pos_data=Hash.new
File.open(ARGV[2]).each_line do |i|
  tmp = i.strip.split(",")
  pos_data[tmp[0]] = tmp[1,8]
end

MaxLine=1000
sam.read_head
out.write_head(sam.head)
rec=Array.new
c_d,c_n=0,0

while true
  STDOUT.flush
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  
  for d in sam.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})$/).first
    next if( k==nil || (d.flag & 256 == 256))
    next unless( pos_data.key?(k) )
    refs=Array.new(2)
#    d.seq = SeqTools.reverse(d.seq) if(d.flag & 16 == 16)
#    rec.push(">"+pos_data[k].join(","))
#    rec.push(d.seq)
    for i in [0,1]
      p_d=pos_data[k][i*4,4]
      refs[i] = SeqTools.find_sequence(ref_file, ref_index, p_d[0], p_d[1].to_i, d.seq.size + ref_margin, p_d[3], ref_sparse)
      refs[i] = SeqTools.reverse(refs[i]) if(p_d[2] == "r")
#      rec.push(refs[i])
    end

    align_ddp=AlignDDP.new
    align_ddp.set_query(d.seq, refs[0], refs[1])

    align_ddp.init_matrix_1
    align_ddp.calc_matrix_1
    align_ddp.calc_trace_back_1
    first_score = align_ddp.fin_score

    align_ddp.init_matrix_2
    align_ddp.calc_matrix_2
    align_ddp.calc_trace_back_2
  
    si,sj = align_ddp.get_start_point
    lbi,lbj = align_ddp.get_left_break_point
    rbi,rbj = align_ddp.get_right_break_point
    ei,ej = align_ddp.get_end_point
    alignment = align_ddp.get_alignment
    cig = align_ddp.get_cigar(" ")
    fin_score = align_ddp.fin_score
    c_d += 2
    next if first_score <= fin_score
    cigs=cig.split(" ")
    next if(cigs.size<2)

#    p alignment, cigs, pos_data[k]
#    print(si," ",sj," ",lbi," ",lbj," ",rbi," ",rbj," ",ei," ",ej,"\n")

    p_d=pos_data[k][0,4]
    d_l=SamField.new(k+":LF", 0, p_d[0], p_d[1].to_i, d.mapq, cigs[0], "*", 0, 0, d.seq[sj...lbj], d.qual[sj...lbj], "")
    if(p_d[2] == "r")
      d_l.flag=16 # if(d.flag & 16 != 16)
      d_l.cigar = SeqTools.reverse_cigar(d_l.cigar)
      d_l.seq = SeqTools.reverse(d_l.seq)
      d_l.qual = d_l.qual.reverse
      d_l.pos = d_l.pos - lbi +1
    else
#      d_l.flag=16 if(d.flag & 16 == 16)
      d_l.pos = d_l.pos + si 
    end

    p_d=pos_data[k][4,4]
    d_r=SamField.new(k+":RF", 0, p_d[0], p_d[1].to_i, d.mapq, cigs[1], "*", 0, 0, d.seq[rbj...ej], d.qual[rbj...ej], "")
    if(p_d[2] == "r")
      d_r.flag=16 # if(d.flag & 16 != 16)
      d_r.cigar = SeqTools.reverse_cigar(d_r.cigar)
      d_r.seq = SeqTools.reverse(d_r.seq)
      d_r.qual = d_r.qual.reverse
      d_r.pos = d_r.pos + refs[1].size - ei
    else
#      d_r.flag=16 if(d.flag & 16 == 16)
      d_r.pos = d_r.pos - refs[1].size + rbi + 1
    end
    d_l.rnext=d_r.rname
    d_l.pnext=d_r.pos
    d_r.rnext=d_l.rname
    d_r.pnext=d_l.pos
    d_r.opt=d_l.opt="AS:i:"+align_ddp.fin_score.to_s

    rec.push(d_l)
    rec.push(d_r)

#    p d.seq.size, d_l.seq.size, d_r.seq.size

    # SamField = Struct.new("SamField", :qname, :flag, :rname, :pos, :mapq, :cigar, :rnext, :pnext, :tlen, :seq, :qual, :opt)
  end

#  c_d += sam.size
  c_n += rec.size
  printf("Write records: %d / %d\r", c_n, c_d)
#  write_data(rec, out)
  out.write_data(rec)
  rec.clear
  sam.clear  
end

printf("\n%d records (%d pair) in %d were witten in %s\n", c_n, c_n/2, c_d, out_file)

ref_file.close
sam.close
out.close
