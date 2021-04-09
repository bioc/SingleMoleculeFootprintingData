\name{NRF1pair.bam}
\alias{NRF1pair.bam}
\alias{NRF1pair.bam.bai}
\alias{EnrichmentRegions_mm10.rds}
\alias{ReferenceMethylation.rds}
\alias{AllCs.rds}

\title{SingleMoleculeFootprintingData}

\description{
This Data package contains r objcets necessary to run some of the functions from the SingleMoleculeFootprinting package.
SingleMoleculeFootprinting is an R package providing functions to analyze Single Molecule Footprinting (SMF) data.
}

\usage{
NRF1pair.bam(metadata = FALSE)
NRF1pair.bam.bai(metadata = FALSE)
EnrichmentRegions_mm10.rds(metadata = FALSE)
ReferenceMethylation.rds(metadata = FALSE)
AllCs.rds(metadata = FALSE)
}

\arguments{
  \item{metadata}{FALSE (default) returns data. TRUE returns metadata}
}

\value{
Returns respectively:
    NRF1pair.bam - Bam file contanining reads covering example NRF1 pair binding locus used for SingleMoleculeFootprinting vignette.
    NRF1pair.bam.bai - Bam index file to Bam file used as example data in SingleMoleculeFootprinting vignette
    EnrichmentRegions_mm10.rds - GRanges obj of mouse genomic regions enriched for SMF signal in genome-wide capture experiments. Can be used to compute bait capture efficiency
    ReferenceMethylation.rds - Reference matrix of genome-wide bulk SMF values for published experiments in mouse cell lines
    AllCs.rds - GRanges obj referencing the genomic context cytosines for mm10
}

\examples{
NRF1pair.bam(metadata = TRUE)
NRF1pair.bam.bai(metadata = TRUE)
EnrichmentRegions_mm10.rds(metadata = TRUE)
ReferenceMethylation.rds(metadata = TRUE)
AllCs.rds(metadata = TRUE)
}
