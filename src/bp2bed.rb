#! /usr/bin/ruby

require "optparse"

def is_region(chr, pos, file_lines)
  for line in file_lines
    arr = line.split(" ").map{|i| i.strip}
    if(chr == arr[0] && arr[1].to_i <= pos && pos <= arr[2].to_i)
       return arr[3]
    end    
  end
  return "."
end

Version="1.1.0"
banner = "Usage: bp2bed.rb [option] <input BP file>\n"

t_sc=0
out_file=""
out_csv=false
bed_region=nil
e_mode=1

opt = OptionParser.new(banner)
opt.on("-s score","score value. (1:Gene Symbol, 2:Dist. Pos.)") {|v| t_sc=v.to_i }
opt.on("-m mode","extract mode (1:OR, 2:AND)") {|v| e_mode=v.to_i}
opt.on("-b file","region BED file") { |v| bed_region=v}
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}
opt.on("--comma", "output by CSV file") {|v| out_csv=true}

opt.parse!(ARGV)

exit if(ARGV.size < 1)

file=ARGV[0]
input=File.open(file)

if(bed_region != nil)
  bed_lines=Array.new
  File.open(bed_region).each_line do |s|
    unless(s =~ /^track name/)
      bed_lines.push(s)
    end
  end
end

if(out_file=="")
  output=STDOUT
else
  output=File.open(out_file,"w")
end

if(out_csv)
  output.printf("chr,start,end,name,score\n")
else
  output.printf("track name=\"%s\" description=\"Converted by sam2bed.rb\" type=bedDetail\n", file)
end

id=rand(100000)

while(l=input.gets)
  arr = l.split(",").map{|i| i.strip}
  return if arr.size < 8
  
  id += 1
  for i in [0,1]
    pp=arr[1+i*4,3]
    rr=arr[5-i*4,3]
    rname=pp[0]
    qname=file.crypt("ID") + ":" + id.to_s + ":" + arr[0] + ":" + pp[2]
    score="."
    if(bed_region != nil)
      gene1=is_region(pp[0], pp[1].to_i, bed_lines)
      gene2=is_region(rr[0], rr[1].to_i, bed_lines) 
      if((e_mode == 1 && (gene1 !="." || gene2 != ".")) || (e_mode == 2 && (gene1 != "." && gene2 != ".")))
        is_ext = true
        score = gene1 if(t_sc==1)
      else
        is_ext=false
      end
    else
      is_ext=true
    end
    if(t_sc==2)
      score=rr[0]+":"+rr[1]+":"+rr[2]
    end
    len = 30
    if(pp[2] == "+")
      pos=pp[1].to_i - len
      pos2=pp[1].to_i
    else
      pos=pp[1].to_i
      pos2=pp[1].to_i + len
    end

    if(is_ext)
      if(out_csv)
        output.printf("%s,%d,%d,%s,%s\n", rname, pos, pos2, qname, score)
      else
        output.printf("%s\t%d\t%d\t%s\t%s\n", rname, pos, pos2, qname, score)
      end
    end
  end
end

input.close
output.close
