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


resolve_genes_range <- function(df) {
  # Some SNPs are mapped to a range of genes instead of a single gene.
  # This function identifies all genes within that range.
  res.df <- tibble()
  gene.ref <- read_tsv("../data/gene_reference.bed",
                       col_names = c("chrom", "start", "end", "gene", "gencode_id"))
    for (i in 1:length(df$Gene)) {
      first_gene <- ""
      second_gene <- ""
      if (str_detect(df$Gene[i], " - ")) {
        first_gene <- str_split(df$Gene[i], " - ", simplify = T)[[1]]
        second_gene <- str_split(df$Gene[i], " - ", simplify = T)[[2]]
      }else if (str_detect(df$Gene[i], ", ")) {
        first_gene <- str_split(df$Gene[i], ", ", simplify = T)[[1]]
        second_gene <- str_split(df$Gene[i], ", ", simplify = T)[[2]]
      }
      first_gene_df <- gene.ref %>% filter(gene == first_gene)
      second_gene_df <- gene.ref %>% filter(gene == second_gene)
      if (dim(first_gene_df)[1] == 0  | dim(second_gene_df)[1] == 0) {
        to.res <- df[i,] %>% mutate(
          Gene = NULL,
          n = 2) %>%
          uncount(n) %>%
          add_column(Gene = c(first_gene, second_gene))
        res.df <- res.df %>% bind_rows(to.res)
      } else {
        positions <- c(first_gene_df$start, first_gene_df$end, second_gene_df$start, second_gene_df$end)
        overlapping.genes <- gene.ref %>%
        filter(chrom == first_gene_df$chrom &
               end >= min(positions) &
               start <= max(positions)) %>% distinct(gene)
        to.res <- df[i,] %>% mutate(
          Gene = NULL,
          n = length(overlapping.genes$gene)) %>%
          uncount(n) %>%
          add_column(Gene = overlapping.genes$gene)
        res.df <- res.df %>% bind_rows(to.res)
        }
  }
  return(res.df %>% distinct())
}


read_tsv(opt$infile) %>%
  resolve_genes_range(drugs %>% filter(grepl("-", Gene) | grepl(", ", Gene))) %>%
  bind_rows(drugs %>% filter(!grepl("-", Gene) & !grepl(",", Gene))) %>%
  filter(!is.na(SNP)) %>%
   mutate(
    SNP = case_when(
      SNP =="rs2412222" ~ "rs321081", # Merged rsID
      TRUE ~ as.character(SNP))) %>%
  write_tsv(opt$out)
  #write_tsv("../drugs/drug_response_combined_curated_resolved_gene_range.txt")
