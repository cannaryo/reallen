#! /usr/bin/ruby
#
# SeqTools
# ver.1.0.1
# This class is used to various operation for sequence data
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

class SeqTools
  def self.invert(seq)
    r=""
    seq.each_char do |i|
      case i
      when "A"
        r+="T"
      when "T"
        r+="A"
      when "G"
        r+="C"
      when "C"
        r+="G"
      else
        r+="N"
      end
    end
    return r
  end

  def self.reverse(seq)
    return self.invert(seq).reverse
  end

  def self.reverse_cigar(cigar)
    comm=cigar.scan(/(\d+[MIDNSHP=X])/)
    buf=""
    if(comm != nil)
      for i in comm.reverse
        buf += i.first
      end
    end
    return buf
  end

  def self.find_sequence(fp, index, rname, pos, size, dir, ref_sparse=nil)
    sp_n=(pos-1)/(ref_sparse-1) unless ref_sparse==nil
    fp.seek(index[rname] + pos + sp_n -1)
    if(dir=="-")
      sz_ac=(pos-1)/(ref_sparse-1) - (pos-size-1)/(ref_sparse-1) + size
    fp.seek(-sz_ac+1, IO::SEEK_CUR) if(dir=="-")
    else
      sz_ac=(pos+size-1)/(ref_sparse-1) - (pos-1)/(ref_sparse-1) + size
    end
    return fp.read(sz_ac).delete("\n")
  end
end
