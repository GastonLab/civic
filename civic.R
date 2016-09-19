library(jsonlite)
library(magrittr)

civic_df <- fromJSON("https://civic.genome.wustl.edu/api/genes/?count=20000")

all_civic_variants <- NULL
for(i in 1:nrow(civic_df$records)){
  gene_entry <- civic_df$records[i,]
  if(nrow(gene_entry$variants[[1]]) > 0){
    df <- gene_entry$variants[[1]] %>% dplyr::select(name,id)
    all_civic_variants <- rbind(all_civic_variants,df)
  }
}

all_civic <- NULL
for(i in 1:nrow(all_civic_variants)){
  try({
    variant_data <- fromJSON(paste0("https://civic.genome.wustl.edu/api/variants/",all_civic_variants[i,]$id))

    api_url <- paste0("https://civic.genome.wustl.edu/api/variants/",all_civic_variants[i,]$id)
    api_variant_id <- all_civic_variants[i,]$id
    genesymbol <- variant_data$entrez_name
    entrezgene <- variant_data$entrez_id
    variant_name <- variant_data$name
    variant_description <- variant_data$description
    chromosome <- variant_data$coordinates$chromosome
    chr_start <- variant_data$coordinates$start
    chr_stop <- variant_data$coordinates$stop
    refbase <- variant_data$coordinates$reference_bases
    altbase <- variant_data$coordinates$variant_bases
    refbuild <- variant_data$coordinates$reference_build

    if(variant_description == ''){
      variant_description <- NA
    }

    if(is.null(chromosome)){
      chromosome <- NA
    }
    if(is.null(chr_start)){
      chr_start <- NA
    }
    if(is.null(chr_stop)){
      chr_stop <- NA
    }
    if(is.null(refbase)){
      refbase <- NA
    }
    if(is.null(altbase)){
      altbase <- NA
    }
    if(is.null(refbuild)){
      refbuild <- NA
    }

    variant_description <- stringr::str_replace_all(variant_description,'\\\n',' ')
    all_variant_groups <- NA
    variant_groups <- character()
    if(!is.null(nrow(variant_data$variant_groups))){
      for(k in 1:nrow(variant_data$variant_groups)){
        variant_groups <- c(variant_groups,variant_data$variant_groups[k,]$name)
      }
    }
    if(length(variant_groups) > 0){
      all_variant_groups <- paste(variant_groups,collapse=", ")
    }
    for(j in 1:nrow(variant_data$evidence_items)){

      variant_df <- data.frame('genesymbol' = genesymbol, 'entrezgene' = entrezgene, 'variant_name' = variant_name, 'variant_description' = variant_description, 'variant_groups' = all_variant_groups, 'api_url' = api_url, 'api_variant_id' = api_variant_id, 'chromosome' = chromosome, 'chr_start' = chr_start, 'chr_stop' = chr_stop, 'refbase' = refbase, 'altbase' = altbase, 'refbuild' = refbuild, stringsAsFactors = F, row.names = NULL)
      

      evidence_item <- variant_data$evidence_items[j,]
      variant_df$evidence_type <- evidence_item$evidence_type
      if(evidence_item$evidence_level == 'A'){
        variant_df$evidence_level <- 'A: Validated'
      }
      if(evidence_item$evidence_level == 'B'){
        variant_df$evidence_level <- 'B: Clinical evidence'
      }
      if(evidence_item$evidence_level == 'C'){
        variant_df$evidence_level <- 'C: Case study'
      }
      if(evidence_item$evidence_level == 'D'){
        variant_df$evidence_level <- 'D: Preclinical evidence'
      }
      if(evidence_item$evidence_level == 'E'){
        variant_df$evidence_level <- 'E: Indirect evidence'
      }
      
      variant_df$evidence_description <- stringr::str_replace_all(evidence_item$description,'\\\n',' ')
      variant_df$evidence_description <- stringr::str_replace_all(evidence_item$description,'\\\r','')
      
      variant_df$evidence_direction <- evidence_item$evidence_direction
      variant_df$clinical_significance <- evidence_item$clinical_significance
      variant_df$variant_origin <- evidence_item$variant_origin
      variant_df$variant_hgvs <- evidence_item$variant_hgvs
      variant_df$status <- evidence_item$status
      variant_df$type <- evidence_item$type
      variant_df$status <- evidence_item$status
      variant_df$pubmed_id <- evidence_item$source$pubmed_id
      variant_df$pubmed_html_link <- paste0("<a href='",evidence_item$source$source_url,"' target='_blank'>",evidence_item$source$citation,"</a>")
      variant_df$disease_name <- evidence_item$disease$name
      variant_df$disease_ontology_id <- evidence_item$disease$doid
      variant_df$disease_url <- evidence_item$disease$url
      variant_df$drug_names <- paste(evidence_item$drugs[[1]]$name, collapse = ", ")
      variant_df$drug_interaction_type <- evidence_item$drug_interaction_type
      variant_df$rating <- evidence_item$rating
      all_civic <- rbind(all_civic, variant_df)

    }
    cat(genesymbol,'\n')
  })
}


all_civic[!is.na(all_civic$variant_origin) & all_civic$variant_origin == "N/A",]$variant_origin <- NA
all_civic[!is.na(all_civic$variant_hgvs) & all_civic$variant_hgvs == "N/A",]$variant_hgvs <- NA
all_civic[!is.na(all_civic$variant_hgvs) & all_civic$variant_hgvs == "",]$variant_hgvs <- NA
all_civic[!is.na(all_civic$clinical_significance) & all_civic$clinical_significance == "NA",]$clinical_significance <- NA

all_civic[!is.na(all_civic$altbase) & all_civic$altbase == "",]$altbase <- NA
all_civic[!is.na(all_civic$refbase) & all_civic$refbase == "",]$refbase <- NA
all_civic[!is.na(all_civic$chr_start) & all_civic$chr_start == "",]$chr_start <- NA
all_civic[!is.na(all_civic$chr_stop) & all_civic$chr_stop == "",]$chr_start <- NA

all_civic[all_civic$drug_names == "",]$drug_names <- NA
all_civic[!is.na(all_civic$clinical_significance) & all_civic$clinical_significance == "N/A",]$clinical_significance <- NA

save(all_civic,file="all_civic.rda")

civic_mutation_variants <- all_civic %>% dplyr::select(genesymbol,entrezgene,variant_name,variant_description,evidence_type,evidence_direction,evidence_level,evidence_description,clinical_significance,variant_origin,pubmed_id,pubmed_html_link,disease_url,disease_ontology_id,disease_name,drug_names,rating,drug_interaction_type,variant_groups,api_url,chromosome,chr_start,chr_stop,refbase,altbase) %>% dplyr::distinct() %>% dplyr::filter(!stringr::str_detect(variant_name,"PHOSPHORYLATION|MISLOCALIZATION|POLYMORPHISM|REARRANGEMENT|DUPLICATION|TRANSLOCATION|FUSION|HYPERMETHYLATION|METHYLATION|SERUM|LOH|Polymorphism|N-TERMINAL"))

civic_mutation_variants <- civic_mutation_variants %>% dplyr::filter(!(stringr::str_detect(variant_name,"-") & !stringr::str_detect(variant_name,"^(DEL|EXON)")))

civic_mutation_variants$alteration_type <- 'MUT'
civic_mutation_variants[stringr::str_detect(civic_mutation_variants$variant_name,"EXPRESSION"),]$alteration_type <- 'EXP'
civic_mutation_variants[stringr::str_detect(civic_mutation_variants$variant_name,"AMPLIFICATION|LOSS|COPY"),]$alteration_type <- 'CNA'

civic_mutation_variants <- civic_mutation_variants %>% dplyr::filter(!is.na(chromosome) & !is.na(chr_start) & !is.na(chr_stop))
civic_mutation_variants[civic_mutation_variants$chr_stop == 167729630,]$chr_stop <- 162729630


# get_pmid_data <- function(pmid){
#   res <- RISmed::EUtilsSummary(pmid, type="esearch", db="pubmed")
#   year <- RISmed::YearPubmed(RISmed::EUtilsGet(res))
#   first_author <- paste(RISmed::Author(RISmed::EUtilsGet(res))[[1]][1,]$LastName," et al.",sep="")
#   journal <- RISmed::ISOAbbreviation(RISmed::EUtilsGet(res))
#   link_url <- paste("<a href='https://www.ncbi.nlm.nih.gov/pubmed/",pmid,"' target='_blank'>",paste(first_author,year,journal,sep=", "),"</a>",sep="")
#   
#   return(link_url)
#   
# }
# 
# load('pmids.rda')
# 
# unique_pmids2 <- data.frame('pmid' = unique(civic_mutation_variants$pubmed_id), stringsAsFactors = F)
# tmp <- dplyr::anti_join(unique_pmids2,unique_pmids,by=c("pmid"))
# 
# i <- 1
# tmp$pubmed_html_link <- NA
# while(i <= nrow(tmp)){
#   tmp[i,"pubmed_html_link"] <- get_pmid_data(tmp[i,]$pmid)
#   cat(i,'\n')
#   i <- i + 1
# }
# 
# unique_pmids <- rbind(unique_pmids,tmp)
# save(unique_pmids,file="pmids.rda")
# 
# civic_mutation_variants <- dplyr::left_join(civic_mutation_variants,unique_pmids,by=c("pubmed_id" = "pmid"))



civic_mutation_variants$chr_start <- as.numeric(civic_mutation_variants$chr_start)
civic_mutation_variants$chr_start <- civic_mutation_variants$chr_start - 1
civic_mutation_variants$region_identifier <- paste("CIVIC",civic_mutation_variants$chromosome,civic_mutation_variants$chr_start,civic_mutation_variants$chr_stop,sep="_")
write.table(dplyr::filter(civic_mutation_variants,!is.na(chromosome)),file="civic.biomarkers.tsv",sep="\t",row.names=F,quote=F,col.names=T)

civic_bed_regions <- civic_mutation_variants %>% dplyr::select(chromosome,chr_start,chr_stop) %>% dplyr::filter(!is.na(chromosome) & chromosome != "") %>% dplyr::distinct() %>% dplyr::arrange(chromosome,chr_start,chr_stop)
chrOrder <- c(as.character(c(1:22)),"X","Y")
civic_bed_regions$chromosome <- factor(civic_bed_regions$chromosome, levels=chrOrder)
civic_bed_regions <- civic_bed_regions[order(civic_bed_regions$chromosome),]
civic_bed_regions$chr_start <- as.numeric(civic_bed_regions$chr_start)
civic_bed_regions$chr_stop <- as.numeric(civic_bed_regions$chr_stop)

civic_bed_regions_sorted <- NULL
for(chrom in chrOrder){
  if(nrow(civic_bed_regions[civic_bed_regions$chromosome == chrom,]) > 0){
    chrom_regions <- civic_bed_regions[civic_bed_regions$chromosome == chrom,]
    chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(chr_start, chr_stop)),]
    civic_bed_regions_sorted <- rbind(civic_bed_regions_sorted, chrom_regions_sorted)
  }
  cat(chrom,'\n')
}


civic_bed_regions_sorted$region_identifier <- paste("CIVIC",civic_bed_regions_sorted$chromosome,civic_bed_regions_sorted$chr_start,civic_bed_regions_sorted$chr_stop,sep="_")
write.table(civic_bed_regions_sorted,file="civic.bed",sep="\t",row.names=F,quote=F,col.names=F)

#bgzip civic.bed
#tabix -p bed civic.bed.gz


