library(QuasR)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)

### 1. example_NRF1pair_amplicon.bam & example_NRF1pair_amplicon.bam.bai

# Raw reads (*.fq) have been aligned using the Bioconductor package QuasR.
prj <- QuasR::qAlign(sampleFile=sampleSheet,
              genome="BSgenome.Mmusculus.UCSC.mm10",
              aligner = "Rbowtie",
              paired="fr",
              bisulfite="undir",
              projectName="prj",
              alignmentParameter = "-e 70 -X 1000 -k 2 --best -strata",
              alignmentsDir="./")
# Each resulting bam file has been subset for the reads covering the example region of interest as follows
# system('samtools view -hb $file "chr6:88105500-88107000" > $Subset_file.bam')
# Multiple bam file subsets were merged as follows
# system('samtools merge NRF1pair.bam $(ls Subset_file*.bam)')
# the resulting file was indexed as follows
# system('samtools index NRF1pair.bam')

### 2. EnrichmentRegions_mm10.rds
# Coordinates provided by capture probes manufacturer

### 3. ReducedRefMat.rds
# Raw reads (*.fq) have been aligned using the Bioconductor package QuasR.
QuasRprj <- QuasR::qAlign(sampleFile=sampleSheet,
              genome="BSgenome.Mmusculus.UCSC.mm10",
              aligner = "Rbowtie",
              paired="fr",
              bisulfite="undir",
              projectName="prj",
              alignmentParameter = "-e 70 -X 1000 -k 2 --best -strata",
              alignmentsDir="./")
# Replicates were merged and deduplicated as follows
# system('java -Xmx4g -jar picard.jar MergeSamFiles I=Replicate1.bam I=Replicate2.bam O=Merged.bam')
# system('java -Xmx4g -jar picard.jar MarkDuplicates INPUT=Merged.bam OUTPUT=Merged_deduplicated.bam METRICS_FILE=Deduplication_stats.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true')
# system('samtools index -@ 16 Merged_deduplicated.bam')
# Bulk-level methylation was retrieved using the Bioconductor package QuasR.
Methylation = QuasR::qMeth(QuasRprj[grep(sample, Samples)], mode="allC", range, collapseBySample = TRUE, keepZero = TRUE)
MethylationMatrix = as.matrix(GenomicRanges::elementMetadata(Methylation)[grepl("_M", colnames(GenomicRanges::elementMetadata(Methylation)))])
colnames(MethylationMatrix) = c("SMF_MM_ES_NO", "SMF_MM_NP_NO", "SMF_MM_TKO_DE")
saveRDS(MethylationMatrix, "ReferenceMethylation.rds")

### 4. AllCreduced.rds
Contexts = c("DGCHN","CGH","NWCGW","GCH","GCG")
lapply(1:21, function(chromosome){

  lapply(Contexts, function(context){

    Cytosines = Biostrings::matchPattern(context, BSgenome.Mmusculus.UCSC.mm10[[chromosome]], fixed = FALSE)
    Cytosines_GR = GenomicRanges::resize(GRanges(seqnames(BSgenome.Mmusculus.UCSC.mm10)[chromosome], Cytosines@ranges), width = 1, fix = ifelse(context == "CGH", "start","center"))
    Cytosines_GR$type = context
    Cytosines_GR

  }) -> Cytosines_GR_allContexts
  names(Cytosines_GR_allContexts) = Contexts
  ov1 = findOverlaps(Cytosines_GR_allContexts$DGCHN, Cytosines_GR_allContexts$GCH)
  ov2 = findOverlaps(Cytosines_GR_allContexts$GCG, Cytosines_GR_allContexts$GCH)
  ov3 = findOverlaps(Cytosines_GR_allContexts$NWCGW, Cytosines_GR_allContexts$CGH)
  ov4 = findOverlaps(Cytosines_GR_allContexts$GCG, Cytosines_GR_allContexts$CGH)
  Cytosines_GR_allContexts$GCH = Cytosines_GR_allContexts$GCH[-unique(c(subjectHits(ov1), subjectHits(ov2)))]
  Cytosines_GR_allContexts$CGH = Cytosines_GR_allContexts$CGH[-unique(c(subjectHits(ov3), subjectHits(ov4)))]
  Cytosines_GR_allContexts_Reduced = sort(Reduce(c, Cytosines_GR_allContexts))
  Cytosines_GR_allContexts_Reduced

}) -> Cytosines_GR_allContexts_Reduced_allChromosomes
AllCreduced = Cytosines_GR_allContexts_Reduced = sort(Reduce(c, Cytosines_GR_allContexts_Reduced_allChromosomes))
saveRDS(AllCreduced, "AllCs.rds")





