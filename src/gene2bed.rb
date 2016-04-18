#! /usr/bin/ruby
#
# gene2bed.rb
#
# Copyright (c) 2016 - Ryo Kanno
#
# This software is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
#

require "optparse"
require_relative "CSVReader.rb"

Version="1.0.1"
banner = "Usage: gene2bed.rb [option] <Gene List>\n+Extract region from gene symbol\n"

opts={"list"=>"", "out"=>"", "annot"=>"annotation/hg19_gene_region.csv"}

opt = OptionParser.new(banner)
opt.on("-l file","specify gene symbols by list file") {|v| opts["list"]=v }
opt.on("-g file","annotation file for gene symbol (default: hg19_gene_region.csv)") {|v| opts["annot"]=v}
opt.on("-o file","output file (default: STDOUT)") {|v| opts["out"]=v}

opt.parse!(ARGV)

exit if(ARGV.size < 1 && opts["list"]=="")

if(opts["list"] != "")
  gene_list=File.open(opts["list"]).each_line.map(&:strip) + ARGV
  track_name=opts["list"]
else
  gene_list=ARGV
  track_name="gene_list"
end

if(opts["out"]=="")
  output=STDOUT
else
  output=File.open(opts["out"],"w")
end

csv=CSVReader.new
csv.open(opts["annot"])

printf("Load annotation data from %s .. ", opts["annot"])
csv.read_all
printf("load %d records\n", csv.data.size)

k1=csv.key("Chromosome")
k2=csv.key("Start")
k3=csv.key("End")
k4=csv.key("Name")
k5=csv.key("ID")
rec=Array.new
c_n=0

for d in csv.data
  if(gene_list.include?(d[k4]))
    rec.push(d)
  end
  c_n += 1
  printf("Find %d records in %d \r", rec.size, c_n)
end

print("\n")
output.printf("track name=\"%s\" description=\"Converted_by_gene2.rb\" type=bedDetail\n",track_name)
for d in rec
    output.printf("%s\t%s\t%s\t%s\t%s\n",d[k1],d[k2],d[k3],d[k4],d[k5])
end

printf("%d records were witten in %s\n", rec.size, opts["out"]) if(opts["out"] != "")

csv.close
output.close
