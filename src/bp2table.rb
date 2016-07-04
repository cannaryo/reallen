#! /usr/bin/ruby
#
# bp2table.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require "optparse"
require_relative "CSVReader.rb"
require_relative "GenomeAnnotation.rb"

# input: str, int, str, str, int, str, int
def find_sv_type(chr1, pos1, ter1, chr2, pos2, ter2, size)
  if(chr1 != chr2)
    return "CTX"
  elsif(size > 0 && (pos1-pos2).abs - 1 > size)
    return "ITX"
  elsif(ter1 == "p" && ter2 == "q")
    if(pos1 < pos2)
      return "DEL"
    else
      return "DUP"
    end
  elsif(ter1 == "q" && ter2 == "p")
    if(pos1 > pos2)
      return "DEL"
    else
      return "DUP"
    end
  elsif((ter1 == "p" || ter1 == "q") && ter1 == ter2)
    return "INV"
  end
  return "Unknown"
end

# input: str, int, Array
def is_region(chr, pos, bed_data)
  for dd in bed_data
    if(chr == dd[0] && dd[1].to_i <= pos && pos <= dd[2].to_i)
       return dd[3]
    end    
  end
  return "."
end

# input: Hash
def make_csv_data(opts)
  csv=CSVReader.new
  cols=["index", "ID","chrom","pos","term","dir","read_side","cov","max_len"]
  cols.push("orig_cov","cov_ratio") if(opts["cov"] != "")
  cols.push("region_name") if(opts["bed"] != "")
  cols.push("gene","geneID","strand") if(opts["gene"] != "")
  cols.push("exon") if(opts["exon"] != "")
  cols.push("cytoband") if(opts["cyto"] != "")
  cols.push("sv_type", "sv_size") if(opts["type"] != -1)
  cols.push("dist_chr","dist_pos","dist_ter") if(opts["dist"])
  cols.push("filter")
  for s in cols
    csv.add_blank_col(s)
  end
  return csv
end

# input: str, Hash
def set_default_files(dir, op)
  return if dir==""
  op["gene"]= dir+"/hg19_gene_region.csv"
  op["exon"]= dir+"/hg19_exon_region.csv"
  op["cyto"]= dir+"/cytoband.csv"
end

# input: CSV, str, Array
def csv_push_main(csv, id, pp)
  # pp : side, chr, pos, dir, ter, cov, len
  dat=Array.new(csv.head.size, 0)
  dat[csv.key("ID")] = id
  arr_order=["read_side","chrom","pos","dir","term","cov","max_len"]
  for i in 0...arr_order.size 
    k=csv.key(arr_order[i])
    dat[k] = pp[i]
  end
  csv.add_data(dat)
  return dat
end

# input: Array, Hash
def load_annotation_tables(keys, opts)
  tables = Hash.new
  for tag in keys
    if(opts[tag] != "")
      STDERR.printf("Load annotation data from %s .. ", opts[tag])
      tables[tag] = GenomeAnnotation.new
      sz = tables[tag].open(opts[tag])    
      STDERR.printf("load %d lines\n", sz)
    end
  end
  return tables
end

# input: CSV, str, str
def find_in_csv(csv, id, side)
  k_id,k_s=csv.key("ID"),csv.key("read_side")
  for d in csv.data
    return d if(d[k_id] ==id && d[k_s]==side)
  end
  return nil
end

# input: int
def sep_comma(num)
  return num.to_s.gsub(/(\d)(?=(\d{3})+(?!\d))/, '\1,')
end

Version="1.3.3"
banner = "Usage: bp2table.rb [option] <input BP file>\n+Output BP file in table format\n"

opts = {"cov"=>"","bed"=>"","mode"=>"OR","out"=>"","cyto"=>"", "gene"=>"", "exon"=>"","type"=> -1,
  "dist"=>false,"blast"=>false, "ncov"=>0, "rcov"=>0.0, "keep"=>false}
out_file = ""
annotation_dir = ""

opt = OptionParser.new(banner)
opt.on("-c file", "output original coverage (specify output of 'samtools bedcov')") { |v| opts["cov"]=v} 
opt.on("-d","output distination position") {|v| opts["dist"]=true  }
opt.on("-g file", "output gene symbol (specify annotation CSV)") { |v| opts["gene"]=v} 
opt.on("-e file", "detect exon region (specify annotation CSV)") { |v| opts["exon"]=v} 
opt.on("-i file","output cytoband (specify annotation CSV)") { |v| opts["cyto"]=v}
opt.on("-o file","output file (default: STDOUT)") {|v| out_file=v}
opt.on("-t [size]", "Estimate SV type (size: maximum SV size for call)") {|v| opts["type"]=v.to_i}
opt.on("-b file","filter by region BED file") {|v| opts["bed"] = v}
opt.on("-m", "mark filter-records instead of removing") {|v| opts["keep"] = true}
opt.on("--and-mode","region filter by AND mode (default: OR mode)") {|v| opts["mode"]="AND"}
opt.on("--annotation dir","use default annotation CSV files in dir(instead of -gei)") {|v| annotation_dir=v} 
opt.on("--blast", "output by BLAST format") {|v| opts["blast"]=true}
opt.on("--cov-filter number", "filter by minimum coverage") {|v| opts["ncov"]=v.to_i}
opt.on("--cov-rate-filter fraction", "filter by coverage rate") {|v| opts["rcov"]=v.to_f}

opt.parse!(ARGV)

exit if(ARGV.size < 1)
printf("\nbp2table.rb %s\n",Version) if(out_file != "")
printf("input=%s\n", ARGV[0]) if(out_file != "")

file=ARGV[0]
input=File.open(file)

if(opts["bed"] != "")
  bed_lines=Array.new
  File.open(opts["bed"]).each_line do |s|
    if(s =~ /^\S+\s+\d+\s+\d+\s+\S+\s+\S+/)
      dat = s.split(" ").map{|i| i.strip}
      bed_lines.push(dat)
    end
  end
end

if(out_file=="")
  output=STDOUT
else
  output=File.open(out_file,"w")
end

set_default_files(annotation_dir, opts)
csv=make_csv_data(opts)

annot_keys=["gene","exon","cyto","cov"]
annot_tables=load_annotation_tables(annot_keys, opts)

index=1
while(l=input.gets)
  arr = l.split(",").map{|i| i.strip}
  exit if arr.size != 15
  
  filter = "PASS"
  reg_names = [".", "."]
  for i in [0,1]
    pp=arr[1+i*7,7]
    rr=arr[8-i*7,7]
    # arr[0]: ID
    # pp : [0] side, [1] chr, [2] pos, [3] dir, [4] ter, [5] cov, [6] len
    if(opts["bed"] != "")
      gene1=is_region(pp[1], pp[2].to_i, bed_lines)
      gene2=is_region(rr[1], rr[2].to_i, bed_lines) 
      unless((opts["mode"] == "OR" && (gene1 !="." || gene2 != ".")) || (opts["mode"] == "AND" && (gene1 != "." && gene2 != ".")))
        filter = "OuterRegion"
      end
      reg_names[i] = gene1
    end
    if(filter == "PASS" && pp[5].to_i < opts["ncov"])
      filter = "LowCoverage"
    end
    if(filter == "PASS" && opts["cov"] != "")
      cov = annot_tables["cov"].find(pp[1],pp[2].to_i,"Coverage").to_f
      if(cov > 0.0 && pp[5].to_f/cov < opts["rcov"])
        filter = "LowCoverageRate"
      end
    end
  end

  if(filter == "PASS" || opts["keep"])
    for i in [0,1]
      pp=arr[1+i*7,7]
      rr=arr[8-i*7,7]
      # arr[0]: ID
      # pp : [0] side, [1] chr, [2] pos, [3] dir, [4] ter, [5] cov, [6] len
      annot_index = Hash.new
      for tag in annot_keys
        if(annot_tables[tag] != nil)
          annot_index[tag] = annot_tables[tag].find_index(pp[1],pp[2].to_i)
        end
      end

      dat=csv_push_main(csv, arr[0], pp)
      for k in 0...csv.head.size
        case csv.head[k]
        when "index"
          dat[k] = index
        when "region_name"
          dat[k] = reg_names[i]
        when "dist_chr"
          dat[k] = rr[1]
        when "dist_pos"
          dat[k] = rr[2]
        when "dist_ter"
          dat[k] = rr[4]
        when "cytoband"
          dat[k] = annot_tables["cyto"].get_data(annot_index["cyto"],"Cytoband")
        when "gene"
          dat[k] = annot_tables["gene"].get_data(annot_index["gene"], "Name")
        when "geneID"
          dat[k] = annot_tables["gene"].get_data(annot_index["gene"], "ID")
        when "strand"
          dat[k] = annot_tables["gene"].get_data(annot_index["gene"], "Strand")
        when "exon"
          dat[k] = annot_tables["exon"].get_data(annot_index["exon"], "Name")
        when "orig_cov"
          dat[k] = annot_tables["cov"].get_data(annot_index["cov"], "Coverage")
        when "cov_ratio"
          if((cov = annot_tables["cov"].get_data(annot_index["cov"], "Coverage").to_f) == 0.0)
            dat[k] = "NA"
          else
            dat[k] = pp[5].to_f/cov
          end
        when "sv_type"
          dat[k] = find_sv_type(pp[1], pp[2].to_i, pp[4], rr[1], rr[2].to_i, rr[4], opts["type"])
        when "sv_size"
          if(pp[1] == rr[1])
            dat[k] = (pp[2].to_i - rr[2].to_i).abs
          else
            dat[k] = "-"
          end
        when "filter"
          dat[k] = filter
        end
      end

    end   
    index += 1
    # end data push 
  end
end


if(opts["blast"])
  if(opts["cyto"]=="")
    STDERR.printf("! Blast format needs cytoband for correct output. Instead, rname is used\n")
    k_c=csv.key("chrom")
  else
    k_c=csv.key("cytoband")
  end
  k_i,k_d,k_p,k_s=csv.key("ID"),csv.key("dir"),csv.key("pos"),csv.key("read_side")
  for d in csv.data
    if(d[k_s]=="L")
      rr=find_in_csv(csv, d[k_i], "R")
      output.printf("%s(%s)(%s)::%s(%s)(%s)\n",d[k_c],d[k_d],sep_comma(d[k_p]),rr[k_c],rr[k_d],sep_comma(rr[k_p]))
    end
  end
  printf("%d records were witten in %s\n", csv.data.size/2, out_file) if(out_file != "")
else
  csv.output(output)
  printf("%d records were witten in %s\n", csv.data.size, out_file) if(out_file != "")
end

input.close
output.close
