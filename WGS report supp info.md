# Cancer Whole Genome Research Report Supplementary information

## Clinical Significance of Variants Based on AMP / ASCO / CAP Guidelines
Adapted from [PMID: 27993330](https://pubmed.ncbi.nlm.nih.gov/27993330/)
### Strong Significance
#### Tier 1A
- Biomarker predicts response or resistance to an FDA or EMA approved therapy, according to drug label or professional guidelines for this diagnosis
- Biomarker included in professional guidelines is prognostic or diagnostic for this diagnosis

#### Tier 1B	
- Biomarker predicts response or resistance to a therapy for this diagnosis based on well powered studies
- Biomarker is prognostic or diagnostic for this diagnosis based on well-powered studies

### Potential Significance
#### Tier 2C
- Biomarker is associated with response or resistance to an FDA or EMA approved therapy, according to drug label or professional guidelines but only for different diagnosis
- Biomarker is an inclusion criterion for an active clinical trial
- Biomarker is prognostic or diagnostic based on multiple small studies

#### Tier 2D
- Biomarker shows plausible response or resistance based on case or preclinical studies
- Biomarker may assist in disease diagnosis or prognosis based on small studies

### Uncertain Significance
#### Tier 3
- Biomarker has uncertain clinical significance and not known to be likely benign or benign

## WGS: Test Methodology & Limitations:
Whole Genome Sequencing was performed on DNA extracted from fresh frozen tissue and matched blood. Libraries were prepared using the Illumina TruSeq Nano library method using 200ng of DNA. Indexed libraries were pooled and sequenced to a depth of 50x for the normal and 100x for tumour using paired 150bp reads on the Illumina Novaseq 6000 platform. Sequence reads were aligned to the build 38 of the genome human reference genome using DRAGEN v3.9. Variants were detected by the DRAGEN mutation caller (https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html). Quality control, variant filtering, virus detection and mutational signature detection are performed using the UMCCR post-processing framework, [umccrise v2.2.0](https://github.com/umccr/umccrise). All variants were annotated using the personalized cancer genome reporter (https://github.com/sigven/pcgr). SNV/InDels in the report were classified according to a five-tiered structure, similar to proposed recommendations [PMID: 24768039](https://pubmed.ncbi.nlm.nih.gov/24768039/), also adopting the MLVD framework for description of clinically relevant cancer variants. Whole Genome Sequencing at 100x sequencing depth, this assay has >95% sensitivity for SNPs and >90% sensitivity for InDels in coding sequences

## Identification of hotspot variants at low allele frequency:
Small nucleotide variants are reported only if they occur at a variant allele frequency ≥10%, except for variants occurring in mutation hotspot regions. Pathogenic hotspot variants in the coding region of the genome (sourced from CIVIC, OncoKb and CGI) are analysed using Somatic Alterations in Genome (SAGE) [PMID: 31645765](https://pubmed.ncbi.nlm.nih.gov/31645765/).  SAGE searches for each hotspot in the tumour BAM files directly and reports evidence of hotspot SNVs and small in frame indels. Variants in these regions are filtered by variant read base and mapping quality, exact variant match, allele frequency, tumour allelic depth and whether detected in the germline.

## Copy number and Structural variant analysis:
Whole genome analysis of Copy number was performed using the PURPLE package. Structural variants were detected using the DRAGEN SV caller. 

## Known genomic targets mediating tumour immune checkpoint inhibitor resistance:
Based on [PMID: 30842677](https://pubmed.ncbi.nlm.nih.gov/30842677/), the following genes are curated for mutations in tumour samples that display at least one of the following features: Microsatellite instability, high tumour mutational burden (>10 mutations/Mb) and ultramutation due to hotspot mutations in the proofreading domains of POLD1 or POLE. B2M, IFNGR1, IFNGR2, JAK1, JAK2, STAT genes, IRF1, SOCS1, PIAS4, CD274, CDK4, SPOP, CMTM4, CMTM6, STK11.
