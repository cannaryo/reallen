#! /usr/bin/ruby
#
# SAMReader
# ver.1.8.1
# This class is used to read SAM file
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

SamField = Struct.new("SamField", :qname, :flag, :rname, :pos, :mapq, :cigar, :rnext, :pnext, :tlen, :seq, :qual, :opt)

class SAMReader
  attr_reader :head
  attr_reader :size
  
  def initialize
    @file=nil
    @head=Array.new
    @data=Hash.new
    @data_dummy=Array.new
    @start=Time.now
    @size=0
  end

  def finalize
    self.close
  end
    
  def open(file)
    @file=File.open(file)
  end

  def read_head
    while(true)
      fp=@file.pos
      break unless(l=@file.gets)
      if(l=~/^@\w\w\s/)
        @head.push(l.chomp)
      elsif(l=~/^(?:\S+\s+){11}(\S.*)?$/)
        @file.pos=fp
        break
      end
    end
    return l
  end
  
  def read_record
    while(l=@file.gets)
      if(l=~/^(?:\S+\s+){11}(\S.*)?$/)
        opt=$1
        d=l.chomp.split(/\s/)
        n=d[2]
        @data_dummy.push(SamField.new(d[0],d[1].to_i,d[2],d[3].to_i,d[4].to_i,d[5],d[6],d[7].to_i,d[8].to_i,d[9],d[10],opt))
        @data[n]=Array.new if @data[n]==nil
        @data[n].push(@data_dummy.last)
        @size += 1
        break
      end
    end
    return l
  end

  def read_next
    while(l=@file.gets)
      if(l=~/^@\w\w\s/)
        @head.push(l.chomp)
      elsif(l=~/^(?:\S+\s+){11}(\S.*)?$/)
        opt=$1
        d=l.chomp.split(/\s/)
        n=d[2]
#        if @data[n]==nil
#          @data_dummy.push(Array.new)
#          @data[n] = @data_dummy.last
#        end
        @data_dummy.push(SamField.new(d[0],d[1].to_i,d[2],d[3].to_i,d[4].to_i,d[5],d[6],d[7].to_i,d[8].to_i,d[9],d[10],opt))
        @data[n]=Array.new if @data[n]==nil
        @data[n].push(@data_dummy.last)
        @size += 1
        break
      end
    end
    return l
  end
  
  def read_all
    while(read_next)
    end
  end

  # chr: rname, start: region pos, size: region size, dir: "reverse" or "forward", mapq: mapping quality
  def filter(chr="all", start=0, size=0, dir=nil, mapq=0)
    arr=Array.new
    for d in self.data(chr)
      if(d.pos <start || (size != 0 && start+size < d.pos))
        next
      elsif((dir=="forward" && (d.flag & 16 == 16)) || (dir=="reverse" && (d.flag & 16 == 0)))
        next
      elsif(d.mapq < mapq)
        next
      end
      arr.push(d)
    end
    return arr
  end
  
  def dump
    for i in @head
      print(i,"\n")
    end
    for i in self.data
      print(i.values.join("  "),"\n")
    end
  end
  
  def header
    return @header
  end

  def data(chr="all")
    if(chr=="all")
      return @data_dummy
    else
      if(@data.key?(chr))
        return @data[chr]
      end
    end
    return []
  end
  
  def clear
    @size=0
    @data.clear
    @data_dummy.clear
  end
  
  def close
    @file.close if @file != nil
  end
end

class SAMWriter
  def initialize(file)
    @file=File.open(file,"w")
    @is_head=false
  end

  def finalize
    self.close
  end

  def write_head(header)
    if @is_head
      print("Header has been already written in files: skip [write_head]\n")
      return
    end
    for d in header
      @file.print(d,"\n")
    end
    @is_head=true
  end

  # chr: rname, start: region pos, size: region size, dir: "reverse" or "forward", mapq: mapping quality
  def write_data(data)
    for d in data
      @file.printf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n",d.qname, d.flag, d.rname, d.pos, d.mapq, d.cigar, d.rnext, d.pnext, d.tlen, d.seq, d.qual, d.opt)
    end
  end

  def close
    @file.close if @file != nil
  end
end

class SAMDecoder
  # dat: struct SamField
  def self.align_sequence(dat)
    comm=dat.cigar.scan(/(\d+)([MIDNSHP=X])/)
    seq=dat.seq
    str=""
    for t in comm
      n=t[0].to_i
      case t[1]
      when "M"
        str+=seq[0,n]
        seq=seq[n..-1]
      when "I"
        seq=seq[n..-1]
      when "D"
        str+="-"*n
      when "N"
        str+="."*n
      when "S"
        seq=seq[n..-1]
      when "H"
      else
        print(" !! ",t,"\n")
      end
    end
    return str
  end

  # dat: struct SamField
  def self.get_aligned_length(dat)
    comm=dat.cigar.scan(/(\d+)([MIDNSHP=X])/)
    size=0
    for t in comm
      n=t[0].to_i
      case t[1]
      when "M"
        size+=n
      when "D","N","S","H"
      when "I"
        size+=n
      end
    end
    return size
  end

  # dat: struct SamField
  def self.get_aligned_ref_length(dat)
    comm=dat.cigar.scan(/(\d+)([MIDNSHP=X])/)
    size=0
    for t in comm
      n=t[0].to_i
      case t[1]
      when "M"
        size+=n
      when "D","N"
        size+=n
      when "I","S","H"
      end
    end
    return size
  end
  
  # dat: struct SamField, tag: string
  def self.get_option(dat,tag)
    dat.opt.split("\t").each do |t|
      c1,c2,c3=t.scan(/(\w+):([AifZHB]):(.+)/).first
      if(c1 != nil && c1 == tag)
        case c2
        when "f"
          c3 = c3.to_f
        when "i"
          c3 = c3.to_i
        end
        return c3
      end
    end
    return 0
  end

  # dat: struct SamField
  def self.get_quality(dat)
    return dat.qual.unpack("C*").map! {|i| i-33}
  end
  
  # arr: Array of struct SamField
  def self.find_min_position(arr)
    if(arr.first)
      pos = arr.first.pos
      end
    for i in arr
      pos = i.pos if pos > i.pos
    end
    return pos
  end
  
  # dat: struct SamField
  def self.get_aligned_quality(dat)
    comm=dat.cigar.scan(/(\d+)([MIDNSHP=X])/)
    qual=get_quality(dat)
    qt=Array.new
    for t in comm
      n=t[0].to_i
      case t[1]
      when "M"
        qt+=qual[0,n]
        qual=qual[n..-1]
      when "I"
        qual=qual[n..-1]
      when "D"
        qt+=[0]*n
      when "N"
        qt+=[0]*n
      when "S"
        qual=qual[n..-1]
      when "H"
      else
        print(" !! ",t,"\n")
      end
    end
    return qt
  end
  
end

