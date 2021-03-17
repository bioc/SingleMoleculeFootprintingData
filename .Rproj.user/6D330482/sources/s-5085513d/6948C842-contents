# ?ExperimentHubData::makeExperimentHubMetadata()

# SpeciesList = AnnotationHubData::getSpeciesList()
# grep("^Mus ", SpeciesList, ignore.case = TRUE, value = TRUE)

# TaxIDs = GenomeInfoDb::loadTaxonomyDb()
# TaxIDs[which(TaxIDs$tax_id == 10090),]

Metadata_DF =
  data.frame(Title = c("example_NRF1pair_amplicon.bam",
                       "example_NRF1pair_amplicon.bam.bai",
                       "EnrichmentRegions_mm10.rds",
                       "ReducedRefMat.rds",
                       "AllCreduced.rds"),
             Description = c("Bam file contanining reads covering example NRF1 pair binding locus used for SingleMoleculeFootprinting vignette",
                             "Bam index file to Bam file used as example data in SingleMoleculeFootprinting vignette",
                             "GRanges obj of mouse genomic regions enriched for SMF signal in genome-wide capture experiments. Can be used to compute bait capture efficiency",
                             "Reference matrix of genome-wide bulk SMF(%) values for published experiments in mouse cell lines",
                             "GRanges obj referencing the genomic context cytosines for mm10"),
             BiocVersion = c("3.13", "3.13", "3.13", "3.13", "3.13"),
             Genome = c("mm10", "mm10", "mm10", "mm10", "mm10"),
             SourceType = c("FASTA", "FASTA", "FASTA", "FASTA", "FASTA"),
             SourceUrl = c("", "", "", "", ""), #<------------
             SourceVersion = c(1, 1, 1, 1, 1),
             Species = c("Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus", "Mus musculus"),
             TaxonomyId = c("10090", "10090", "10090", "10090", "10090"),
             Coordinate_1_based = c(FALSE, FALSE, FALSE, FALSE, FALSE),
             DataProvider = c("Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>", "Arnaud Krebs <arnaud.krebs@embl.de>"),
             Maintainer = c("Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>", "Guido Barzaghi <guido.barzaghi@embl.de>"),
             RDataClass = c("character", "character", "GRanges", "matrix", "GRanges"),
             DispatchClass = c("FilePath", "FilePath", "Rds", "Rds", "Rds"),
             Location_Prefix = c("", "", "", "", ""), #<------------
             RDataPath = c("", "", "", "", ""), #<------------
             Tags = c("", "", "", "", "")) #<------------

#final data are written out with
write.csv(Metadata_DF, file = "inst/extdata/metadata.csv", row.names=FALSE)

