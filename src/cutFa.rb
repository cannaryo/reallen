#! /usr/bin/ruby

require "SeqTools.rb"
require "optparse"

Version="1.2.0"
banner = "Usage: cutFa.rb [option] <in.fa> <region file>\n+Cut sequence\n"

opts={"qname"=>"seq","show_format"=>false,"index"=>nil, "sparse"=>nil,"out"=>"tmp.fa"}

opt = OptionParser.new(banner)
opt.on("-i ref_index", "Reference index file (to skip indexing)") {|v| opts["index"]=v}
opt.on("-q qname", "qname root (default: seq)") {|v| opts["qname"]=v} 
opt.on("-o file", "ouput file (default: tmp.fa)") {|v| opts["out"]=v} 
opt.on("--reference-sparse N", "set sparse code each N charactor in reference") {|v| opts["sparse"]=v.to_f}
opt.on("--file-format","show file format and exit") {|v| opts["show_format"] = true}
opt.parse!(ARGV)

if(ARGV.size < 2 || opts["show_format"])
  printf("CutFa.rb cut sequences in fasta file and make merged one\n")
  printf("[-h] shows details information\n")
  printf("----------------------\n")
  printf("region file format is:\n")
  printf("qname,start,end, +/-\n")
  printf("qname,start,end, +/-\n")
  printf("*\n")
  printf("qname,start,end, +/-\n...\n")
  printf("----------------------\n")
  printf("# each line is merged to one sequence\n")
  printf("# '*' is separator (>qname)\n")
  exit
end


ref_file=File.open(ARGV[0])
region_file=File.open(ARGV[1])
ref_index=Hash.new

if(opts["index"]==nil)
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
  File.open(opts["index"]).each_line do |l|
    if(l =~ /^(\w+)\t\d+\t(\d+)/)
      ref_index[$1] = $2.to_i
    end
  end
end

# test reference sparse
if(opts["sparse"] == nil)
  while(s=ref_file.gets)
    t_1 = ref_file.gets
    t_2 = ref_file.gets
    if(t_1.size == s.size && t_2.size == s.size)
      opts["sparse"] = s.size
      break
    end
  end
end

cc=1
out=File.open(opts["out"],"w")
out.print(">",opts["qname"],cc,"\n")
ref = ""
for s in region_file
  if(s[0,1]=="*")
    printf("%s%d : %d bp\n",opts["qname"],cc,ref.size)
    ms=opts["sparse"]-1
    for i in 1..(ref.size/ms)
      ref.insert(i*ms+i-1,"\n")
    end
    out.print(ref,"\n")
    ref = ""
    cc += 1
    out.print(">",opts["qname"],cc,"\n")
  else    
    dd=s.split(" ").map(&:strip)
    next if(dd.size != 4)
    printf("chr=%s, start=%s, end=%s, strand=%s\n",*dd[0,4])
    sz =dd[2].to_i - dd[1].to_i + 1
    r_tmp = SeqTools.find_sequence(ref_file, ref_index, dd[0], dd[1].to_i, sz, "+", opts["sparse"])
    r_tmp=SeqTools.reverse(r_tmp) if(dd[3] == "-")
    ref += r_tmp
  end
end

printf("%s%d : %d bp\n",opts["qname"],cc,ref.size)
ms=opts["sparse"]-1
for i in 1..(ref.size/ms)
  ref.insert(i*ms+i-1,"\n")
end
out.print(ref,"\n")

printf("%d sequences was witten in %s\n", cc, opts["out"])

ref_file.close
region_file.close
out.close
