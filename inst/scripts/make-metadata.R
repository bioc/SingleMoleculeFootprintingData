# ?ExperimentHubData::makeExperimentHubMetadata()

# SpeciesList = AnnotationHubData::getSpeciesList()
# grep("^Mus ", SpeciesList, ignore.case = TRUE, value = TRUE)

# TaxIDs = GenomeInfoDb::loadTaxonomyDb()
# TaxIDs[which(TaxIDs$tax_id == 10090),]

Metadata_DF =
  data.frame(Title = c("NRF1pair.bam",
                       "NRF1pair.bam.bai",
                       "EnrichmentRegions_mm10.rds",
                       "ReferenceMethylation.rds",
                       "AllCs.rds"),
             Description = c("Bam file contanining reads covering example NRF1 pair binding locus used for SingleMoleculeFootprinting vignette",
                             "Bam index file to Bam file used as example data in SingleMoleculeFootprinting vignette",
                             "GRanges obj of mouse genomic regions enriched for SMF signal in genome-wide capture experiments. Can be used to compute bait capture efficiency",
                             "Reference matrix of genome-wide bulk SMF(%) values for published experiments in mouse cell lines",
                             "GRanges obj referencing the genomic context cytosines for mm10"),
             BiocVersion = c("3.13", "3.13", "3.13", "3.13", "3.13"),
             Genome = c("mm10", "mm10", "mm10", "mm10", "mm10"),
             SourceType = c("FASTA", "FASTA", "FASTA", "FASTA", "FASTA"),
             SourceUrl = c("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/004/ERR4078784/ERR4078784_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/004/ERR4078784/ERR4078784_2.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/005/ERR4078785/ERR4078785_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/005/ERR4078785/ERR4078785_2.fastq.gz",
                           "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/004/ERR4078784/ERR4078784_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/004/ERR4078784/ERR4078784_2.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/005/ERR4078785/ERR4078785_1.fastq.gz,ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR407/005/ERR4078785/ERR4078785_2.fastq.gz",
                           "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/*",
                           "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/*",
                           "ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/*"),
             SourceVersion = c(1, 1, 1, 1, 1),
             Species = c("Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus"),
             TaxonomyId = c("10090", "10090", "10090", "10090", "10090"),
             Coordinate_1_based = c(FALSE, FALSE, FALSE, FALSE, FALSE),
             DataProvider = c("Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>"),
             Maintainer = c("Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>"),
             RDataClass = c("character", "character", "GRanges", "matrix", "GRanges"),
             DispatchClass = c("FilePath", "FilePath", "Rds", "Rds", "Rds"),
             Location_Prefix = c("", "", "", "", ""), #<------------
             RDataPath = c("SingleMoleculeFootprintingData/example_NRF1pair_amplicon.bam",
                           "SingleMoleculeFootprintingData/example_NRF1pair_amplicon.bam.bai",
                           "SingleMoleculeFootprintingData/EnrichmentRegions_mm10.rds",
                           "SingleMoleculeFootprintingData/ReducedRefMat.rds",
                           "SingleMoleculeFootprintingData/AllCreduced.rds"), #<------------
             Tags = c("", "", "", "", ""))

#final data are written out with
write.csv(Metadata_DF, file = "inst/extdata/metadata.csv", row.names=FALSE)
