#! /usr/bin/ruby
#
# filtBP.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require_relative "SAMReader.rb"
require "optparse"

BPField = Struct.new("BPField", :id, :chr, :pos, :dir, :len, :side, :ter)

# Chrom_order = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM"]

def calc_bp_field(dat)
  k,tag=dat.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
  rec=BPField.new(k, dat.rname, dat.pos, "*", 0, "*", "*")
  if(tag =~ /L[MSF]$/)
    if(dat.flag & 16 == 16)
      rec.dir="-"
      rec.ter="q"
    else
      rec.dir="+"
      rec.ter="p"
      rec.pos += SAMDecoder.get_aligned_ref_length(dat) -1
    end
    rec.side = "L"
  elsif(tag =~ /R[MSF]$/)
    if(dat.flag & 16 == 16)
      rec.dir="-"
      rec.ter="p"
      rec.pos += SAMDecoder.get_aligned_ref_length(dat) -1
    else
      rec.dir="+"
      rec.ter="q"
    end
    rec.side = "R"
  end
  rec.len = SAMDecoder.get_aligned_ref_length(dat)
  return rec
end

def bp_compare(dat1, dat2, len)
# c1,c2=Chrom_order.index(dat1.chr),Chrom_order.index(dat2.chr)
  c1,c2=dat1.chr,dat2.chr
  if(c1<c2)
    return 1
  elsif(c1>c2)
    return -1
  else
    if(dat1.pos + len < dat2.pos)
      return 1
    elsif(dat1.pos - len > dat2.pos)
      return -1
    else
      if(dat1.ter < dat2.ter)
        return 1
      elsif(dat1.ter > dat2.ter)
        return -1
      else
        return 0
      end
    end
  end
end

def bp_equal?(dat1, dat2, len)
  return (dat1.chr == dat2.chr) && (dat1.pos + len >= dat2.pos && dat1.pos - len <= dat2.pos) && (dat1.ter == dat2.ter)
end

def add_to_node(hash, d1, d2, len)
  bp1=calc_bp_field(d1)
  bp2=calc_bp_field(d2)
  if(bp_compare(bp1, bp2, 0) == -1)
    bp1, bp2 = bp2, bp1
    d1,d2=d2,d1
  end
  key=bp1.chr+bp1.ter+bp2.chr+bp2.ter
  new_data=true
  if(hash.key?(key))
    for nn in hash[key]
      if(nn.is_range?(bp1,bp2,len))
        nn.add_bp(bp1, bp2)
        nn.add_data(d1, d2)
        new_data=false
        break
      end
    end
    if(new_data)
      bp=BPNode.new
      bp.add_bp(bp1,bp2)
      bp.add_data(d1,d2)
      hash[key].push(bp)
      return 1
    end
  else
    hash[key]=Array.new
    bp=BPNode.new
    bp.add_bp(bp1,bp2)
    bp.add_data(d1,d2)
    hash[key].push(bp)
    return 1
  end
  return 0
end

def filter_by_cover(hash, cov)
  result=Array.new
  for arr in hash.values
    for i in arr
      if(i.cov>=cov)
        result.push(i)
      end
    end
  end
  return result
end

class BPNode
  attr_reader :cov

  def initialize
    @data_arr=Array.new
    @bp_arr1=Array.new
    @bp_arr2=Array.new
    @cov=0
    @bp1=nil
    @bp2=nil
  end

  def is_range?(bp1, bp2, len)
    return true if(bp_compare(bp1, @bp1, len) == 0 && bp_compare(bp2, @bp2, len) == 0)
  end

  def add_bp(bp1, bp2)
    @bp_arr1.push(bp1)
    @bp_arr2.push(bp2)
    @cov+=1
    recalc_center
  end

  def add_data(data1, data2)
    @data_arr.push(data1)
    @data_arr.push(data2)
  end

  def recalc_center
    sorted = @bp_arr1.sort_by {|i| i.pos}
    pp = @bp_arr1.index(sorted[(@bp_arr1.size-1)/2])
    @bp1=@bp_arr1[pp]
#    @bp_arr2 = @bp_arr2.sort_by {|i| i.pos}
    @bp2=@bp_arr2[pp]
  end

  def max_len(bp_arr)
    m=0
    for i in bp_arr
      m=i.len if m<i.len
    end
    return m
  end

  def get_data
    return @data_arr
  end

  def show_by_pos(out)
    k,tag=@data_arr.first.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    out.printf("%s, %s, %s, %d, %s, %s, %d, %d, %s, %s, %d, %s, %s, %d, %d\n", 
               k, @bp1.side, @bp1.chr, @bp1.pos, @bp1.dir, @bp1.ter, @cov, max_len(@bp_arr1), 
               @bp2.side, @bp2.chr, @bp2.pos, @bp2.dir, @bp2.ter, @cov, max_len(@bp_arr2))
  end

  def show_by_bed(out, len=1)
    k,tag=@data_arr.first.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    for dd in [@bp1, @bp2]
      if(dd.ter=="p")
        p1 = dd.pos - len
        p2 = dd.pos
      else
        p1 = dd.pos - 1
        p2 = dd.pos + len -1
      end
      score = dd.side + dd.dir + dd.ter + ":"  + (@cov).to_s
      out.printf("%s\t%d\t%d\t%s\t%s\n",dd.chr,p1,p2,k, score)
    end
  end
end


Version="1.2.8"
banner = "Usage: filtBP.rb [option] <input SAM file>\nFilter by break-point group\n"

opts = {"cov"=>2, "size"=>10, "sam"=>"", "bp"=>"tmp.bp", "bed"=>""}

opt = OptionParser.new(banner)
opt.on("-c cover", "minimum coverage (default: 2)") {|v| opts["cov"]=v.to_i }
opt.on("-r size", "BP redundancy (default: 10)") {|v| opts["size"]=v.to_i}
opt.on("-o file", "output SAM file (default: None)") {|v| opts["sam"]=v}
opt.on("-b file", "output BED file (default: None)") {|v| opts["bed"]=v}
opt.on("-p file", "output BP position file (default: tmp.bp)") {|v| opts["bp"]=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nfiltBP.rb %s\n",Version)
printf("input=%s\n", ARGV[0])

file = ARGV[0]
sam=SAMReader.new
sam.open(file)

bp_out=File.open(opts["bp"], "w")

MaxLine=1000
sam.read_head
rec=Array.new
nodes = Hash.new
ids = Hash.new
c_d,c_n=0,0

while true
  MaxLine.times { sam.read_record }
  break if(sam.size==0)
  for d in sam.data
    k,tag=d.qname.scan(/^([^:]+(?::[^:]+){1,9})[:\/]([^:]*)$/).first
    next if( k==nil || (d.flag & 256 == 256))
    if(ids.key?(k))
      if(d.rname != "*" && ids[k].rname != "*")
        c_n+=add_to_node(nodes, d, ids[k], opts["size"])
      end
      ids.delete(k)
    else
      ids[k] = d
    end
  end  
  c_d += sam.size
  STDERR.printf("filtBP > Find BP node: %d / %d sequences\r", c_n, c_d)
  sam.clear 
end

rec_bp=filter_by_cover(nodes, opts["cov"])

STDERR.print("\n")
printf("%d nodes in %d were written in %s\n",rec_bp.size, c_n, opts["bp"])

for d in rec_bp
  rec.concat(d.get_data)
  d.show_by_pos(bp_out)
end

unless(opts["bed"] =="")
  File.open(opts["bed"],"w") do |fp|
    fp.printf("track name=\"%s\" description=\"Filtered_BP_position\" type=bedDetail\n", file)
    for d in rec_bp
      d.show_by_bed(fp, opts["size"])
    end
  end
  printf("%d nodes in %d were written in %s\n",rec_bp.size, c_n, opts["bed"])
end

unless(opts["sam"] == "")
  out=SAMWriter.new(opts["sam"])
  out.write_head(sam.head)
  out.write_data(rec)
  out.close
  printf("%d sequences in %d were written in %s\n", rec.size, c_d, opts["sam"])
end

sam.close
bp_out.close
