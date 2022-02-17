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




gwas <- read_tsv(opt$infile) %>% 
mutate(CHR_POS = as.character(CHR_POS)) %>% 
  separate_rows( # Split single entries that have multiple SNPs into distinct rows
    CHR_ID, CHR_POS, MAPPED_GENE, `STRONGEST SNP-RISK ALLELE`, SNPS, CONTEXT,sep=";") %>% 
  mutate(
    CHR_ID = str_trim(CHR_ID), 
    CHR_POS = str_trim(CHR_POS), 
    MAPPED_GENE = str_trim(MAPPED_GENE), 
    `STRONGEST SNP-RISK ALLELE` = str_trim(`STRONGEST SNP-RISK ALLELE`), 
    SNPS = str_trim(SNPS), 
    CONTEXT = str_trim(CONTEXT),
    CONTEXT = replace_na(CONTEXT, "Not annotated")
  ) %>% 
 distinct(SNPS, .keep_all = TRUE)


gwas %>% 
  filter(grepl("response", tolower(`MAPPED_TRAIT`)) & !is.na(`OR or BETA`)) %>% 
  select(MAPPED_GENE, `STRONGEST SNP-RISK ALLELE`, `INITIAL SAMPLE SIZE`,
  		      MAPPED_TRAIT, PUBMEDID, `OR or BETA`, `95% CI (TEXT)`, SNPS) #%>%
  #write_tsv(opt$out)
