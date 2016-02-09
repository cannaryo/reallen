#! /usr/bin/ruby
#
# GenomeAnnotation.rb
# ver 1.0.3
# This class is used for searching in annotation table quickly
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require 'csv'

class GenomeAnnotation
  attr_reader :header
  attr_reader :table
  attr_reader :key
  
  ColChr="Chromosome"
  ColStr="Start"
  ColEnd="End"

  def initialize
    @table=nil
    @header=nil
    @key=Hash.new
    @chr_p=Hash.new
    @chr_s=Hash.new
    @index=nil
  end

  def open(file)
    @header, *@table=CSV.read(file)
    @header.each_with_index { |d, i| @key[d] = i}
    return 0 if @table.size==0
    p_key=@table[0][@key[ColChr]]
    @chr_p[p_key]=0
    @table.each_with_index do |d, i|
      c=d[@key[ColChr]]
      unless(@chr_p.include?(c))
        @chr_p[c]=i
        @chr_s[p_key]=i-@chr_p[p_key]
        p_key=c
      end
    end
    @chr_s[p_key]=@table.size - @chr_p[p_key]
    return @table.size
  end

  def get_data_in_chr(chr)
    if(@chr_p[chr] != nil && @chr_s[chr] != nil)
      data=@table[@chr_p[chr], @chr_s[chr]]
      return data
    end
    return []
  end
  
  def find_index(chr,pos)
    data=self.get_data_in_chr(chr)
    st=@key[ColStr]
    en=@key[ColEnd]
    pp = (0...data.size).bsearch {|i| data[i][en].to_i >= pos}
    if(pp != nil && data[pp][st].to_i <= pos && pos <= data[pp][en].to_i)
      return pp + @chr_p[chr]
    else
      return nil
    end
  end
  
  def find(chr, pos, name)
    i=find_index(chr,pos)
    return get_data(i, name) if(i != nil)
    return nil
  end

  def get_data(idx, name)
    if(idx != nil && (k=@key[name]) != nil)
      return @table[idx][k]
    else
      return ""
    end
  end

end
