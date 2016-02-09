#! /usr/bin/ruby
#
# FastqReader.rb 
# ver.1.2.0
# This class is used to read/write fastq file
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

FastqField = Struct.new("FastqField", :id, :seq, :dummy, :qual)

class FastqReader
  attr_reader :data
  attr_reader :line

  def initialize
    @data=Array.new
    @line=0
    @file=nil
  end
  
  def open(file_name)
    @file=File.open(file_name, "r")
  end
  
  def read_next
    n=0
    dat=FastqField.new("","","","")
    while l=@file.gets
      @line+=1
      if(n == 2)
        dat.qual=l.strip
        @data.push(dat)
        break
      elsif(l[0,1] == "+")
        dat.dummy="+"
        n=2
      elsif(n == 1)
        dat.seq = dat.seq + l.strip
      elsif(l[0,1] == "@")        
        dat.id=l.strip
        n=1
      end
    end
    return l
  end
  
  def read_all
    while(read_next)
    end
  end
  
  def size
    return @data.size
  end

  def dump
    for i in @data
      print(i.values.join(" "),"\n")
    end
  end

  def output(out=STDOUT)
    for i in @data
      out.printf("%s\n%s\n%s\n%s\n",i.id, i.seq, i.dummy, i.qual)
    end
  end

  def add_data(dat)
    @data.push(dat)
  end

  def clear
    @data.clear
  end

  def close
    @file.close if(@file!=nil)      
  end

  def self.get_quality(dat)
    return dat.qual.unpack("C*").map! {|i| i-33}
  end

  def finalize
    self.close
  end

end
