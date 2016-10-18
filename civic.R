library(jsonlite)
library(magrittr)

civic_df <- fromJSON("https://civic.genome.wustl.edu/api/genes/?count=20000")

civic_variant_names <- NULL
for(i in 1:nrow(civic_df$records)){
  gene_entry <- civic_df$records[i,]
  if(nrow(gene_entry$variants[[1]]) > 0){
    df <- gene_entry$variants[[1]] %>% dplyr::select(name,id)
    civic_variant_names <- rbind(civic_variant_names,df)
  }
}

all_civic <- NULL
for(i in 1:nrow(civic_variant_names)){
  try({
    variant_data <- fromJSON(paste0("https://civic.genome.wustl.edu/api/variants/",civic_variant_names[i,]$id))

    api_url <- paste0("https://civic.genome.wustl.edu/api/variants/",civic_variant_names[i,]$id)
    api_variant_id <- civic_variant_names[i,]$id
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
all_civic[!is.na(all_civic$clinical_significance) & all_civic$clinical_significance == "NA",]$clinical_significance <- NA

all_civic[!is.na(all_civic$altbase) & all_civic$altbase == "",]$altbase <- NA
all_civic[!is.na(all_civic$refbase) & all_civic$refbase == "",]$refbase <- NA
all_civic[!is.na(all_civic$chr_start) & all_civic$chr_start == "",]$chr_start <- NA
all_civic[!is.na(all_civic$chr_stop) & all_civic$chr_stop == "",]$chr_start <- NA

all_civic[all_civic$drug_names == "",]$drug_names <- NA
all_civic[!is.na(all_civic$clinical_significance) & all_civic$clinical_significance == "N/A",]$clinical_significance <- NA

all_civic$evidence_description <- stringr::str_replace_all(all_civic$evidence_description,"\\\n","")

save(all_civic,file="all_civic.rda")

load(file="all_civic.rda")
civic_all_mut_cna_exp <- all_civic %>% dplyr::select(genesymbol,entrezgene,variant_name,variant_description,evidence_type,evidence_direction,evidence_level,evidence_description,clinical_significance,variant_origin,pubmed_id,pubmed_html_link,disease_url,disease_ontology_id,disease_name,drug_names,rating,drug_interaction_type,variant_groups,api_url,chromosome,chr_start,chr_stop,refbase,altbase) %>% dplyr::distinct() %>% dplyr::filter(!stringr::str_detect(variant_name,"PHOSPHORYLATION|MISLOCALIZATION|POLYMORPHISM|REARRANGEMENT|DUPLICATION|TRANSLOCATION|FUSION|HYPERMETHYLATION|METHYLATION|SERUM|LOH|Polymorphism|N-TERMINAL"))

civic_all_mut_cna_exp$variant_name <- stringr::str_trim(civic_all_mut_cna_exp$variant_name,side="both")

civic_all_mut_cna_exp <- civic_all_mut_cna_exp %>% dplyr::filter(!(stringr::str_detect(variant_name,"-") & !stringr::str_detect(variant_name,"^(DEL|EXON)")))

civic_all_mut_cna_exp$alteration_type <- 'MUT'
civic_all_mut_cna_exp[stringr::str_detect(civic_all_mut_cna_exp$variant_name,"EXPRESSION"),]$alteration_type <- 'EXP'
civic_all_mut_cna_exp[stringr::str_detect(civic_all_mut_cna_exp$variant_name,"AMPLIFICATION|LOSS|COPY"),]$alteration_type <- 'CNA'
civic_all_mut_cna_exp$civic_id <- paste0('CIVIC_',rep(1:nrow(civic_all_mut_cna_exp)))

## PREPARE INPUT FOR TRANSVAR - PROTEIN TO GENOME MAPPING TOOL
protein_variants <- civic_all_mut_cna_exp[stringr::str_detect(civic_all_mut_cna_exp$variant_name,"^[A-Z]{1}[0-9]{1,}[a-zA-Z]{0,}(\\*[0-9]{0,}){0,}$"),] %>% dplyr::select(genesymbol,variant_name) %>% dplyr::distinct()
protein_variants$id <- paste(protein_variants$genesymbol,paste0('p.',protein_variants$variant_name),sep=":")
protein_variants$id <- stringr::str_replace(protein_variants$id,"FS$","fs")
protein_variants$id <- stringr::str_replace(protein_variants$id,"DEL","del")
protein_variants$id <- stringr::str_replace(protein_variants$id,"INS","ins")
protein_variants <- dplyr::rename(protein_variants, transvar_id = id)
civic_all_mut_cna_exp <- dplyr::left_join(civic_all_mut_cna_exp,protein_variants)
transvar_df <- data.frame('id' = protein_variants$transvar_id)
write.table(transvar_df, file="civic.transvar.input.v3.tsv",sep="\t",quote=F,col.names = F,row.names = F)

## RUN TRANSVAR IN VCFANNO DOCKER, PARSE OUTPUT HERE

transvar_output_df <- as.data.frame(readr::read_tsv(file="civic.transvar.output.v3.tsv",col_names = F))
transvar_output <- as.data.frame(stringr::str_split_fixed(transvar_output_df$X5,"/",n = 3))
transvar_output$V4 <- transvar_output_df$X1
transvar_output$V5 <- stringr::str_replace(transvar_output_df$X2," \\(protein_coding\\)","")
transvar_output <- transvar_output %>% dplyr::filter(!stringr::str_detect(V1,"coordinates") & stringr::str_detect(V1,">|del[A-Z]{1,}ins[A-Z]{1,}")) %>% dplyr::distinct() %>% dplyr::select(V1,V4,V5)
colnames(transvar_output) <- c('gdna','transvar_id','transcript_id')
#transvar_output <- as.data.frame(transvar_output %>% dplyr::group_by(gdna,transvar_id) %>% dplyr::summarise(transvar_ensembl_transcript_id = paste(transcript_id, collapse=",")))


transvar_output_part2 <- dplyr::select(transvar_output_df, X1,X2,X7)
transvar_output_part2$X2 <- stringr::str_replace(transvar_output_part2$X2," \\(protein_coding\\)","")
transvar_output_part2$X3 <- stringr::str_replace(stringr::str_match(transvar_output_df$X7,"candidate_snv_variants=chr[0-9]{1,}:g\\.[0-9]{1,}(A|G|C|T)>(A|G|C|T)")[,1],"candidate_snv_variants=","")
transvar_output_part2 <- dplyr::select(transvar_output_part2, -X7) %>% dplyr::filter(X1 != 'input' & !is.na(X3))
colnames(transvar_output_part2) <- c('transvar_id','transcript_id','gdna')
transvar_output <- as.data.frame(rbind(transvar_output, transvar_output_part2) %>% dplyr::group_by(gdna,transvar_id) %>% dplyr::summarise(transvar_ensembl_transcript_id = paste(transcript_id, collapse=",")))



##insertion/deletions
transvar_output_indels <- transvar_output %>% dplyr::filter(stringr::str_detect(gdna,"del|ins"))
transvar_output_indels$chr_start <- stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna,"chr.{1,}:g\\.",''),"_",n=2)[,1]
alleles <- stringr::str_split_fixed(stringr::str_replace(transvar_output_indels$gdna,'chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}_[0-9]{1,}del',''),'ins',2)
transvar_output_indels$refbase <- alleles[,1]
transvar_output_indels$altbase <- alleles[,2]

transvar_output_indels$chr_stop <- as.integer(transvar_output_indels$chr_start) + nchar(transvar_output_indels$altbase) - 1
transvar_output_indels$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_indels$gdna,':',2)[,1],'chr','')

##snvs
transvar_output_snvs <- transvar_output %>% dplyr::filter(!stringr::str_detect(gdna,"del|ins"))
transvar_output_snvs$chr_start <- stringr::str_replace(stringr::str_replace(transvar_output_snvs$gdna,'chr([0-9]{1,}|X|Y):g\\.',''),'(A|G|C|T)>(A|G|C|T)$','')
transvar_output_snvs$chr_stop <- transvar_output_snvs$chr_start
alleles <- stringr::str_split_fixed(stringr::str_replace(transvar_output_snvs$gdna,'chr([0-9]{1,}|X|Y):g\\.[0-9]{1,}',''),'>',2)
transvar_output_snvs$refbase <- alleles[,1]
transvar_output_snvs$altbase <- alleles[,2]
transvar_output_snvs$chromosome <- stringr::str_replace(stringr::str_split_fixed(transvar_output_snvs$gdna,':',2)[,1],'chr','')

all_transvar <- rbind(transvar_output_indels,transvar_output_snvs)


## ENTRIES MAPPED TO VARIANT ALLELES (USING TRANSVAR)
civic_vcf_mapped_transvar <- as.data.frame(dplyr::left_join(dplyr::select(dplyr::filter(civic_all_mut_cna_exp,!is.na(transvar_id)),civic_id,transvar_id),all_transvar) %>% dplyr::filter(!is.na(gdna)) %>% dplyr::select(chromosome,civic_id,chr_start,refbase,altbase) %>% dplyr::group_by(chromosome,chr_start,refbase,altbase) %>% dplyr::summarise(CIVIC_ID = paste(civic_id, collapse = ",")))
civic_vcf <- as.data.frame(dplyr::left_join(dplyr::select(dplyr::filter(civic_all_mut_cna_exp,!is.na(transvar_id)),civic_id,transvar_id),all_transvar) %>% dplyr::filter(!is.na(gdna)) %>% dplyr::select(civic_id) %>% dplyr::distinct())
civic_vcf_mapped_transvar <- dplyr::rename(civic_vcf_mapped_transvar, chrom = chromosome, pos = chr_start, ref = refbase, alt = altbase)

## ENTRIES MAPPED WITH VARIANT ALLELES FROM CIVIC (NOT MAPPED USING TRANSVAR)
civic_vcf_mapped_original <- dplyr::anti_join(civic_all_mut_cna_exp, civic_vcf) %>% dplyr::filter(alteration_type == 'MUT') %>% dplyr::filter(!is.na(chromosome) & !is.na(chr_start) & !is.na(refbase) & !is.na(altbase) & altbase != '-' & refbase != '-' & !stringr::str_detect(altbase,'/')) %>% dplyr::select(chromosome,chr_start,refbase,altbase,civic_id)
civic_vcf_mapped_original <- dplyr::rename(civic_vcf_mapped_original, chrom = chromosome, pos = chr_start, ref = refbase, alt = altbase)
tmp <- dplyr::select(civic_vcf_mapped_original, civic_id) %>% dplyr::distinct()
tmp2 <- as.data.frame(civic_vcf_mapped_original %>% dplyr::group_by(chrom,pos,ref,alt) %>% dplyr::summarise(CIVIC_ID = paste(civic_id, collapse=",")))

civic_vcf <- rbind(tmp, civic_vcf) %>% dplyr::distinct()
civic_vcf_mapped <- rbind(civic_vcf_mapped_transvar, tmp2)

civic_vcf <- dplyr::left_join(dplyr::select(civic_vcf,civic_id),dplyr::select(civic_all_mut_cna_exp,alteration_type,genesymbol,civic_id))
civic_vcf$mapping_category <- 'exact'
civic_vcf$mapping_rank <- 1
civic_vcf$civic_exon <- NA
civic_vcf$civic_codon <- NA
civic_vcf$civic_consequence <- NA

civic_bed <- dplyr::anti_join(civic_all_mut_cna_exp, civic_vcf) %>% dplyr::filter(!is.na(chromosome) & !is.na(chr_start) & !is.na(chr_stop) & (is.na(refbase) | is.na(altbase) | altbase == '-' | refbase == '-')) %>% dplyr::select(genesymbol,variant_name,alteration_type,transvar_id,civic_id,chromosome,chr_start,chr_stop)

civic_bed_gene_mutations <- civic_bed %>% dplyr::filter(alteration_type == 'MUT' & (variant_name == 'MUTATION'|  variant_name == 'DELETION' | variant_name == 'FRAMESHIFT MUTATION' | variant_name == 'FRAMESHIFT TRUNCATION' | variant_name == 'DELETERIOUS MUTATION' | variant_name == 'TRUNCATING MUTATION'))
civic_bed_gene_mutations$mapping_category <- 'gene'
civic_bed_gene_mutations$civic_codon <- NA
civic_bed_gene_mutations$civic_consequence <- NA
civic_bed_gene_mutations[civic_bed_gene_mutations$variant_name == 'DELETION',]$civic_consequence <- 'inframe_deletion'
civic_bed_gene_mutations[civic_bed_gene_mutations$variant_name == 'TRUNCATING MUTATION',]$civic_consequence <- 'stop_gained'
civic_bed_gene_mutations[stringr::str_detect(civic_bed_gene_mutations$variant_name,'FRAMESHIFT '),]$civic_consequence <- 'frameshift_variant'

civic_bed_gene_mutations$civic_exon <- NA
civic_bed_gene_mutations$mapping_rank <- 4
civic_bed_gene_mutations <- civic_bed_gene_mutations %>% dplyr::select(-c(transvar_id,variant_name))

civic_bed_gene_cna <- civic_bed %>% dplyr::filter(alteration_type == 'CNA' & stringr::str_detect(variant_name,"LOSS|AMPLIFICATION"))
civic_bed_gene_cna$mapping_category <- 'gene'
civic_bed_gene_cna$civic_codon <- NA
civic_bed_gene_cna$civic_consequence <- NA
civic_bed_gene_cna[stringr::str_detect(civic_bed_gene_cna$variant_name,"LOSS"),]$civic_consequence <- 'loss'
civic_bed_gene_cna[stringr::str_detect(civic_bed_gene_cna$variant_name,"AMPLIFICATION"),]$civic_consequence <- 'gain'
civic_bed_gene_cna$civic_exon <- NA
civic_bed_gene_cna$mapping_rank <- 1
civic_bed_gene_cna <- civic_bed_gene_cna %>% dplyr::select(-c(transvar_id,variant_name))

civic_bed_codon <- civic_bed %>% dplyr::filter(alteration_type == 'MUT' & stringr::str_detect(variant_name,"^[A-Z]{1}[0-9]{1,}((FS|fs)(\\*[0-9]{1,}){0,}){0,}$"))
civic_bed_codon$mapping_category <- 'codon'
civic_bed_codon$civic_codon <- as.integer(stringr::str_replace_all(civic_bed_codon$variant_name,"(^[A-Za-z])|(FS$)|fs\\*[0-9]{1,}$",""))
civic_bed_codon$civic_consequence <- NA
civic_bed_codon[stringr::str_detect(civic_bed_codon$variant_name,"(FS|fs\\*[0-9])$"),]$civic_consequence <- 'frameshift_variant'
civic_bed_codon$civic_exon <- NA
civic_bed_codon$mapping_rank <- 2
civic_bed_codon <- civic_bed_codon %>% dplyr::select(-c(transvar_id,variant_name))

civic_bed_exon <- civic_bed %>% dplyr::filter(alteration_type == 'MUT' & stringr::str_detect(variant_name,"EXON "))
civic_bed_exon$mapping_category <- 'exon'
civic_bed_exon$civic_codon <- NA
civic_bed_exon$civic_consequence <- NA
civic_bed_exon[stringr::str_detect(civic_bed_exon$variant_name," DELETION"),]$civic_consequence <- 'inframe_deletion'
civic_bed_exon[stringr::str_detect(civic_bed_exon$variant_name," FRAMESHIFT"),]$civic_consequence <- 'frameshift_variant'

tmp <- dplyr::select(civic_bed_exon,civic_id,variant_name)
tmp$variant_name <- stringr::str_trim(stringr::str_replace_all(tmp$variant_name,"EXON |([A-Z]| )*$",""),side="both")
single_exon <- tmp[!stringr::str_detect(tmp$variant_name,"-"),]
single_exon <- dplyr::rename(single_exon, civic_exon = variant_name)
civic_bed_exon_single <- dplyr::left_join(dplyr::filter(single_exon, stringr::str_detect(civic_exon,"^[0-9]{1,}$")),civic_bed_exon)

multiple_exons <- tmp[stringr::str_detect(tmp$variant_name,"-"),]

multiple_exons_expanded <- NULL
for(n in 1:nrow(multiple_exons)){
  start <- as.integer(stringr::str_split_fixed(multiple_exons[n,]$variant_name,"-",2)[,1])
  stop <- as.integer(stringr::str_split_fixed(multiple_exons[n,]$variant_name,"-",2)[,2])
  for(j in start:stop){
    multiple_exons_expanded <- rbind(multiple_exons_expanded, data.frame('civic_id' = multiple_exons[n,]$civic_id, 'civic_exon' = as.character(j), stringsAsFactors = F))
  }
}

civic_bed_exon_multiple <- dplyr::left_join(multiple_exons_expanded,civic_bed_exon)

civic_bed_exon <- rbind(dplyr::select(civic_bed_exon_multiple,-c(transvar_id,variant_name)),dplyr::select(civic_bed_exon_single,-c(transvar_id,variant_name)))
civic_bed_exon$mapping_rank <- 3

civic_bed_nonexact <- civic_bed %>% dplyr::filter(alteration_type != 'EXP' & alteration_type != 'CNA' & !stringr::str_detect(variant_name,"EXON ") & !stringr::str_detect(variant_name,"^[A-Z]{1}[0-9]{1,}((FS|fs)(\\*[0-9]{1,}){0,}){0,}$") & !stringr::str_detect(variant_name,'MUTATION|DELETION|DELETERIOUS|TRUNCATING|FRAMESHIFT|VIII|ITD|WILD')) %>% dplyr::select(-c(variant_name,transvar_id))

civic_bed_nonexact$mapping_category <- 'gene'
civic_bed_nonexact$civic_codon <- NA
civic_bed_nonexact$civic_consequence <- NA
civic_bed_nonexact$civic_exon <- NA
civic_bed_nonexact$mapping_rank <- 5

civic_bed <- rbind(civic_bed_gene_mutations,civic_bed_exon,civic_bed_codon,civic_bed_gene_cna,civic_bed_nonexact)
civic_bed_mapped <- civic_bed %>% dplyr::select(chromosome,chr_start,chr_stop,civic_id)
civic_bed <- civic_bed %>% dplyr::select(-c(chromosome,chr_start,chr_stop))
civic_bed_mapped$chr_start <- as.integer(civic_bed_mapped$chr_start) - 1

civic_biomarkers1 <- dplyr::inner_join(dplyr::select(civic_all_mut_cna_exp,-c(chromosome,chr_start,chr_stop,refbase,altbase)),civic_vcf)
civic_biomarkers2 <- dplyr::inner_join(dplyr::select(civic_all_mut_cna_exp,-c(chromosome,chr_start,chr_stop,refbase,altbase)),civic_bed)
civic_biomarkers <- rbind(civic_biomarkers1,civic_biomarkers2)
civic_biomarkers_rest <- dplyr::anti_join(dplyr::select(civic_all_mut_cna_exp,-c(chromosome,chr_start,chr_stop,refbase,altbase)),civic_biomarkers)
tags <- c('mapping_category','mapping_rank','civic_exon','civic_codon','civic_consequence')
for(t in tags){
  civic_biomarkers_rest[t] <- NA
}
civic_biomarkers_all <- rbind(civic_biomarkers,civic_biomarkers_rest)
                                     
#civic_all_mut_cna_exp$chr_start <- as.numeric(civic_all_mut_cna_exp$chr_start)
#civic_all_mut_cna_exp$chr_start <- civic_all_mut_cna_exp$chr_start - 1
#civic_all_mut_cna_exp$region_identifier <- paste("CIVIC",civic_all_mut_cna_exp$chromosome,civic_all_mut_cna_exp$chr_start,civic_all_mut_cna_exp$chr_stop,sep="_")
write.table(civic_biomarkers_all,file="civic.biomarkers.tsv",sep="\t",row.names=F,quote=F,col.names=T)

chrOrder <- c(as.character(c(1:22)),"X","Y")
civic_bed_mapped$chromosome <- factor(civic_bed_mapped$chromosome, levels=chrOrder)
civic_bed_mapped <- civic_bed_mapped[order(civic_bed_mapped$chromosome),]
civic_bed_mapped$chr_start <- as.numeric(civic_bed_mapped$chr_start)
civic_bed_mapped$chr_stop <- as.numeric(civic_bed_mapped$chr_stop)

civic_bed_mapped_sorted <- NULL
for(chrom in chrOrder){
  if(nrow(civic_bed_mapped[civic_bed_mapped$chromosome == chrom,]) > 0){
    chrom_regions <- civic_bed_mapped[civic_bed_mapped$chromosome == chrom,]
    chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(chr_start, chr_stop)),]
    civic_bed_mapped_sorted <- rbind(civic_bed_mapped_sorted, chrom_regions_sorted)
  }
  cat(chrom,'\n')
}

write.table(civic_bed_mapped_sorted,file="civic.bed",sep="\t",row.names=F,quote=F,col.names=F)

#bgzip civic.bed
#tabix -p bed civic.bed.gz

header_lines <- c("##fileformat=VCFv4.2","##SOURCE_CIVIC=2016_10_19","##INFO=<ID=CIVIC_ID,Number=.,Type=String,Description=\"Identifier for evidence item of clinical biomarker in CiVIC database\">","#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
write(header_lines,file="civic.vcf",sep="\n")

civic_vcf_mapped$CIVIC_ID <- paste0('CIVIC_ID=',civic_vcf_mapped$CIVIC_ID)
civic_vcf_mapped$QUAL <- '.'
civic_vcf_mapped$FILTER <- 'PASS'
civic_vcf_mapped$ID <- '.'
civic_vcf_mapped <-  dplyr::rename(civic_vcf_mapped, CHROM = chrom, POS = pos, REF = ref, ALT = alt, INFO = CIVIC_ID)

civic_vcf_mapped <- civic_vcf_mapped[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]

write.table(civic_vcf_mapped, file="civic_vcfcontent.tsv",sep="\t",col.names = F,quote=F, row.names = F)

system("cat civic_vcfcontent.tsv | tail -n +2 | egrep -v \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 >> civic.vcf")
system("cat civic_vcfcontent.tsv | tail -n +2 | egrep \"^[XYM]\" | sort -k1,1n -k2,2n -k4,4 -k5,5 >> civic.vcf")
system("bgzip -c civic.v2.vcf > civic.vcf.gz")
system("tabix -p vcf civic.vcf.gz")

