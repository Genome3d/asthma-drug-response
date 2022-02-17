#!/usr/bin/env Rscript
library("optparse")

option_list = list(
    make_option(c("-i", "--infile"), type="character", default=NULL,
    		help="", metavar="character")
    make_option(c("-out", "--out"), type="character", default=NULL,
    		help="", metavar="character")
		);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

print(opt$infile)
print(opt$infile[1])
print(opt$infile[2])
