#! /usr/bin/ruby
#
# AlignDDP.rb 
# ver.1.1.0
# This class is used to perform 
# "Discontinuous Dynamic Programming" algorithm
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

class AlignDDP
  private
  def calc_score(f1,f2)
    return @score_match if(f1==f2)
    return @score_mismatch
  end
  
  def trace_back(ali,ei,ej)
    i=ei
    j=ej
    cg=""
    while true
      if(ali[i][j] == "M")
        cg += ali[i][j]
        i-=1
        j-=1
      elsif(ali[i][j] == "I")
        cg += ali[i][j]
        j-=1
      elsif(ali[i][j] == "D")
        cg += ali[i][j]
        i-=1
      elsif(ali[i][j] == "F")
        return cg.reverse,i,j
      elsif(ali[i][j] == "B")
        cg+="B"
        return cg.reverse,i,j
      else
        return "No path found",0,0
      end
    end
  end
  
  def find_end_point(mat,sx,sy)
    score=0
    ei=sx
    ej=sy
    for i in 0..sx
      for j in 0..sy
        if(mat[i][j] <= score)
          score = mat[i][j]
          ei=i
          ej=j
        end
      end
    end
    return ei,ej
  end

  public
  attr_accessor :score_match, :score_mismatch, :score_gap
  attr_accessor :score_no_startclip, :score_endclip, :score_jump
  attr_reader :fin_score, :is_bp

  def initialize()
    @score_match=-1
    @score_mismatch=3
    @score_gap=5
    @score_startclip=-2
    @score_endclip=-2
    @score_jump=5

    @mat, @ali=nil,nil
    @mat_2, @ali_2=nil,nil
    @seq, @ref_1, @ref_2=nil,nil,nil
    @jump_array=nil
    @fin_score=0
    @alignment=""
    @si, @sj, @ei, @ej=0,0,0,0
    @lbi, @lbj, @rbi,@rbj=0,0,0,0
    @is_bp=false
  end

  def set_query(query, ref_1, ref_2=nil)
    @seq=query
    @ref_1=ref_1
    @ref_2=ref_2
  end

  def init_matrix_1
    sy=@seq.size
    sx=@ref_1.size
    @mat=Array.new(sx+1).map{Array.new(sy+1, 0)}
    @ali=Array.new(sx+1).map{Array.new(sy+1,"F")}    
    for i in 0..sx
      @mat[i][0] = @score_startclip
    end
  end

  def calc_matrix_1
    return if(@mat==nil)
    sy=@seq.size
    sx=@ref_1.size
    for i in 1..sx
      for j in 1..sy
        s_i = @mat[i-1][j]+@score_gap
        s_j = @mat[i][j-1]+@score_gap
        s_ij = @mat[i-1][j-1]+calc_score(@ref_1[i-1,1], @seq[j-1,1])
        s_st = 0
        @mat[i][j]=[s_i,s_j,s_ij,s_st].min
        if(s_ij == @mat[i][j])
          @ali[i][j]="M"
        elsif(s_i == @mat[i][j])
          @ali[i][j]="D"
        elsif(s_j == @mat[i][j])
          @ali[i][j]="I"
        else
          @ali[i][j]="F"
        end
      end
    end
    for i in 0..sx
      @mat[i][sy] += @score_endclip
    end
  end

  def calc_trace_back_1
    return if(@mat==nil)
    sy=@seq.size
    sx=@ref_1.size
    @ei,@ej = find_end_point(@mat,sx,sy)
    @fin_score = @mat[@ei][@ej]
    @alignment, @si, @sj =trace_back(@ali,@ei,@ej)
  end

  def get_alignment
    return @alignment
  end

  def get_start_point
    return @si,@sj
  end

  def get_end_point
    return @ei,@ej
  end

  def get_left_break_point
    return @lbi, @lbj
  end
  
  def get_right_break_point
    return @rbi, @rbj
  end

  def init_matrix_2
    return if(@ref_2==nil || @mat==nil)
    sy=@seq.size
    sx=@ref_1.size
    @jump_arr=Array.new(sy)
    for j in 0..sy
      score,pi=100,0
      for i in 0..sx
        if(@mat[i][j] <= score)
          score = @mat[i][j]
          pi=i
        end
        @jump_arr[j]=[score + @score_jump, pi]
      end
    end
    sx=@ref_2.size
    @mat2=Array.new(sx+1).map{Array.new(sy+1,0)}
    @ali2=Array.new(sx+1).map{Array.new(sy+1,"F")}
    for i in 0..sx
      @mat2[i][0] = @score_startclip
    end
  end

  def calc_matrix_2
    return if(@mat==nil || @mat2==nil)
    sy=@seq.size
    sx=@ref_2.size
    for i in 1..sx
      for j in 1..sy
        s_i = @mat2[i-1][j]+@score_gap
        s_j = @mat2[i][j-1]+@score_gap
        s_ij = @mat2[i-1][j-1]+calc_score(@ref_2[i-1,1], @seq[j-1,1])
        s_jp = @jump_arr[j][0]
        s_st = 0
        @mat2[i][j]=[s_i,s_j,s_ij,s_jp,s_st].min
        if(s_ij == @mat2[i][j])
          @ali2[i][j]="M"
        elsif(s_i == @mat2[i][j])
          @ali2[i][j]="D"
        elsif(s_j == @mat2[i][j])
          @ali2[i][j]="I"
        elsif(s_jp == @mat2[i][j])
          @ali2[i][j]="B"
        else
          @ali2[i][j]="F"
        end
      end
    end
    for i in 0..sx
      @mat2[i][sy] = @mat2[i][sy] + @score_endclip
    end
    
    def calc_trace_back_2
      return if(@mat==nil || @mat2==nil)
      sy=@seq.size
      sx=@ref_2.size
      ei,ej = find_end_point(@mat2,sx,sy)
      score = @mat2[ei][ej]
      return if(@fin_score <= score)
      @fin_score=score
      align, bi, bj =trace_back(@ali2,ei,ej)
      @ei, @ej = ei, ej
      if(align[0,1] == "B")
        @rbi, @rbj=bi, bj
        @lbi, @lbj=@jump_arr[bj][1], bj
        @is_break=true
        tmp_s, @si, @sj=trace_back(@ali, @lbi, @lbj)
        @alignment=tmp_s+align
      else
        @si,@sj=bi,bj
      end
    end
  end

  def get_cigar(sep="")
    return "" if @alignment==""
    cigar=""
    c=""
    j=0
    for i in @alignment.chars
      if(c==i)
        j += 1
        next
      end
      if(c == "B")
        cigar += sep
      elsif(c != "")
        cigar += j.to_s+c
      end
      c=i
      j=1
    end
    cigar += j.to_s+c
    return cigar
  end

end
