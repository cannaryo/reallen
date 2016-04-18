#! /usr/bin/ruby
#
# CSV reader
# ver.1.5.1
# This class is used to read CSV data
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require 'csv'

class CSVReader
  attr_reader :head
  attr_reader :data

  def initialize
    @file
    @head=Array.new
    @data=Array.new
    #@sep=","
  end
  
  def open(f)
    @file=File.open(f)
    if(s=@file.gets)
      #@sep="\"," if s=~/".+"\s*,/
      #@head=s.split(@sep).map! {|x| x.strip.delete("\"")}
      @head=CSV.parse_line(s)
    end
  end
  
  def read_next
    if(s=@file.gets)
      #t=s.split(@sep)
      #t.map! {|x| x.strip.delete("\"")}
      @data.push(CSV.parse_line(s))
    else
      nil
    end
  end
  
  def read_all
    while(read_next)
    end
  end
  
  def delete_data
    @data.clear
  end

  def clear
    @data.clear
  end
  
  def output(f=STDOUT)
    once=true
    @head.each do |i|
      f.print(",") unless once
      f.print("\"",i,"\"")
      once=false
    end
    f.print("\n")
    for d in @data
      once=true
      d.each do |i|
        f.print(",") unless once
        f.print("\"",i,"\"")
        once=false
      end
      f.print("\n")
    end
  end
  
  def key(val)
    @head.index(val)
  end
  
  def val(n, k)
    @data[n][key(k)]
  end
  
  def set_head(hed)
    @head=hed
  end
  
  def add_data(dat)
    @data.push(dat)
  end
  
  def close
    @data.clear
    @head.clear
    @file.close
  end
  
  def add_blank_col(name)
    @head.push(name)
    for d in @data
      d.push("")
    end
    return key(name)
  end
  
  def delete_col(name)
    k=key(name)
    @head.delete_at(k)
    for d in @data
      d.delete_at(k)
    end
  end
  
  def expand_row(n)
    if(n>@data.size)
      (n-@data.size).times do |i|
        dat=Array.new(@head.size,"")
        @data.push(dat)
      end
    end
  end
  
  def reduce_row(n)
    if(n<@data.size)
      for i in n...@data.size
        @data.delete_at(i)
      end
    end
  end
  
  def self.copy_col_data(src, s_k, dest, d_k)
    s_k=src.key(s_k) if s_k.instance_of?(String)
    d_k=dest.key(d_k) if d_k.instance_of?(String)
    return if(s_k == nil || d_k == nil)
    dest.data.size.times do |i|
      if(src.data.size > i)
        dest.data[i][d_k] = src.data[i][s_k]
      else
        dest.data[i][d_k] = ""
      end
    end
  end
end
