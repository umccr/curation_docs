
# UMCCR Somatic Whole Genome Sequencing Report Supplementary information

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

## Test Methodology & Limitations:
### Whole Genome Sequencing (WGS)
Whole genome sequencing was performed on DNA extracted from fresh or formalin-fixed paraffin-embedded (FFPE) tissue and matched blood. Libraries were prepared using the Illumina TruSeq Nano library method using 200ng of DNA. Indexed libraries were pooled and sequenced to a depth of 40x for the normal and 80x (fresh) or 120x (FFPE) for tumour using paired end 150bp reads on the Illumina Novaseq 6000 platform. 

Secondary analysis of WGS data was performed using the Illumina [DRAGEN Bio-IT Platform v4.4](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) and [Oncoanalyser v2.2](https://github.com/nf-core/oncoanalyser). Sequence reads were aligned using [DRAGEN](https://github.com/umccr/cwl-ica/releases/tag/dragen-wgts-dna-pipeline%2F4.4.4__20250617224149) to the GRCh38 human reference genome with mapping, variant calling, and quality control completed using default settings. Variant calling, Oncovirus detection, Transcript analysis, and Structural Variant (SV) event interpretation were performed by Oncoanalyser in “WGTS” mode using DRAGEN BAMs and RNA reads. 

Tertiary analysis of WGS data (variant annotation and prioritization/filtering, fusion detection, mutational signatures, biomarker characterization) was carried out with [sash v0.6](https://github.com/scwatts/sash) and oncoanalyser v2.2. 

### Single Nucleotide Variants (SNV) and Insertion-Deletions (Indels)
SNVs and Indels were annotated using the [personalized cancer genome reporter](https://github.com/sigven/pcgr).  SNV/Indels were prioritised according to a four-tiered structure, adopting the joint consensus recommendation by AMP/ACMG ([PMID: 27993330](https://pubmed.ncbi.nlm.nih.gov/27993330/)). SNVs are reported only if they occur at a variant allele frequency ≥10%, with 4 independent IGV reads.
Whole Genome Sequencing at effective 80x sequencing depth, this assay has >95% sensitivity for SNPs and >90% sensitivity for Indels in coding sequences.

### Identification of hotspot variants at low allele frequency
High-confidence variants at allele frequency below the validated 10% may be included in the report if they occur in mutation hotspot regions, however these results should be interpreted with caution. Pathogenic hotspot variants in the coding region of the genome (sourced from CIVIC, OncoKb and CGI) are analysed using Somatic Alterations in Genome (SAGE) [PMID: 31645765](https://pubmed.ncbi.nlm.nih.gov/31645765/).  SAGE searches for each hotspot in the tumour BAM files directly and reports evidence of hotspot SNVs and small in frame indels. Variants in these regions are filtered by variant read base and mapping quality, exact variant match, allele frequency, tumour allelic depth and whether detected in the germline.

### Copy number and Structural variant analysis
Whole genome analysis of Copy number was performed using the PURPLE package. Structural variants were detected using the [eSVee](https://github.com/hartwigmedical/hmftools/tree/master/esvee) SV caller and annotated using SnpEff. 

### Known genomic targets mediating tumour immune checkpoint inhibitor resistance
Based on [PMID: 33819567](https://pubmed.ncbi.nlm.nih.gov/33819567/), the following genes are curated for mutations in tumour samples that display at least one of the following features: Microsatellite instability, high tumour mutational burden (>10 mutations/Mb) and ultramutation due to hotspot mutations in the proofreading domains of POLD1 or POLE: B2M, IFNGR1, IFNGR2, JAK1, JAK2, STK11.

### CNV & SV annotation
CNV/SV changes are annotated with the svprioritize framework, assigning priority to fusion events, whole exon loss or upstream/downstream changes for a list of 806 cancer-associated genes only. The genes reviewed for CNV and SVs: ABCB1, ABI1, ABI2, ABL1, ABL2, ABRAXAS1, ACSL3, ACTA2, ACTB, ACVR1, ACVR1B, ACVR2A, AFDN, AFF1, AFF3, AFF4, AHCYL1, AIP, AJUBA, AKAP9, AKT1, AKT2, AKT3, ALB, ALDH2, ALK, AMER1, ANKRD26, APC, AR, ARAF, ARHGAP26, ARHGAP35, ARHGEF12, ARHGEF7, ARID1A, ARID1B, ARID2, ARID5B, ARMC5, ASPSCR1, ASXL1, ASXL2, ATF1, ATF7IP, ATIC, ATM, ATR, ATRX, AURKA, AURKB, AXIN1, AXIN2, AXL, B2M, BAP1, BARD1, BAX, BCL10, BCL11B, BCL2, BCL2L1, BCL2L11, BCL2L14, BCL6, BCL9, BCOR, BCORL1, BCR, BEND2, BIRC3, BLM, BMPR1A, BORCS5, BRAF, BRCA1, BRCA2, BRD3, BRD4, BRIP1, BTG1, BTK, BUB1B, CAD, CALR, CAMTA1, CANT1, CARD11, CARS1, CASP8, CASR, CBFA2T3, CBFB, CBL, CCDC6, CCNB1IP1, CCNB3, CCND1, CCND2, CCND3, CCNE1, CD19, CD274, CD28, CD74, CD79A, CD79B, CDC73, CDH1, CDK12, CDK4, CDK6, CDK8, CDKN1A, CDKN1B, CDKN1C, CDKN2A, CDKN2B, CDKN2C, CDX2, CEBPA, CEP83, CEP89, CHCHD7, CHD1, CHD4, CHEK1, CHEK2, CHIC2, CIC, CIITA, CLIP1, CLTC, CNBP, CNOT3, COL1A1, COL2A1, COX6C, CPSF6, CREB1, CREB3L1, CREB3L2, CREBBP, CRKL, CRLF2, CRTC1, CRTC3, CSDE1, CSF1R, CSF3R, CTC1, CTCF, CTLA4, CTNNA1, CTNNB1, CTR9, CUL1, CUL3, CUX1, CXCR4, CYLD, CYSLTR2, DAB2IP, DAXX, DAZAP1, DCTN1, DCUN1D1, DDB2, DDIT3, DDR2, DDX10, DDX3X, DDX41, DDX5, DEK, DGCR8, DICER1, DIS3, DIS3L2, DKC1, DNAJB1, DNMT1, DNMT3A, DOT1L, DPYD, DROSHA, DUSP22, DUSP4, DUX4, DUX4L1, E2F3, EBF1, EED, EGFR, EGLN1, EGLN2, EGR1, EHF, EIF1AX, EIF3E, EIF4A1, EIF4A2, ELANE, ELF3, ELK4, ELL, ELOC, EML4, ENG, EP300, EPAS1, EPCAM, EPHA2, EPHA3, EPHA7, EPHB1, EPOR, EPS15, ERBB2, ERBB3, ERBB4, ERC1, ERCC1, ERCC2, ERCC3, ERCC4, ERCC5, ERCC6, ERG, ERRFI1, ESR1, ESRP1, ETNK1, ETS1, ETV1, ETV4, ETV5, ETV6, EWSR1, EXT1, EXT2, EZH1, EZH2, EZR, FAM131B, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, FAS, FAT1, FBXW7, FEV, FGF1, FGF19, FGF3, FGF4, FGFR1, FGFR2, FGFR3, FGFR4, FH, FHIT, FIP1L1, FLCN, FLI1, FLT1, FLT3, FLT4, FOXA1, FOXA2, FOXL2, FOXO1, FOXO3, FOXO4, FOXP1, FUBP1, FUS, FYN, GAB2, GAS7, GATA1, GATA2, GATA3, GATA4, GATA6, GLI1, GLIS2, GMPS, GNA11, GNA13, GNAQ, GNAS, GOLGA5, GOPC, GPC3, GPHN, GREM1, GRIN2A, GSK3B, GTF2I, H1-2, H1-3, H1-4, H1-5, H2BC5, H3-3A, H3-3B, H3-5, H3C2, HERPUD1, HEY1, HFE, HGF, HIF1A, HIP1, HLA-A, HLA-B, HLA-C, HLF, HMGA1, HMGA2, HNF1A, HNRNPA2B1, HNRNPC, HNRNPUL1, HOOK3, HOXA9, HOXB13, HRAS, ID3, IDH1, IDH2, IFNGR1, IGF1R, IGH, IGK, IGL, IKBKE, IKZF1, IKZF2, IKZF3, IL2RA, IL2RB, IL3, IL6R, IL7R, INPP4B, IRF2, IRF2BP2, IRF4, IRF8, IRS2, ITK, JAK1, JAK2, JAK3, JARID2, JAZF1, JMJD1C, JUN, KANK1, KAT6A, KDM5A, KDM5C, KDM6A, KDR, KEAP1, KIAA1549, KIF1B, KIF5B, KIT, KLF4, KLF5, KLF6, KLK2, KMT2A, KMT2B, KMT2C, KMT2D, KRAS, KTN1, LAG3, LASP1, LATS1, LATS2, LCK, LIFR, LIMA1, LMNA, LMO1, LMO2, LMO3, LPP, LRP1B, LSM14A, LYL1, LYN, LZTR1, MAD2L2, MAF, MAFB, MALT1, MAML2, MAP2K1, MAP2K2, MAP2K4, MAP3K1, MAP3K13, MAP3K8, MAPK1, MAPK3, MAX, MBD4, MCL1, MCPH1, MDH2, MDM2, MDM4, MECOM, MED12, MEF2B, MEF2D, MEN1, MET, MGA, MIR143, MITF, MLF1, MLH1, MLLT1, MLLT10, MLLT11, MLLT3, MLLT6, MN1, MNX1, MPL, MRE11, MRTFA, MSH2, MSH3, MSH6, MSN, MST1, MST1R, MTAP, MTOR, MUS81, MUTYH, MYB, MYC, MYCL, MYCN, MYD88, MYH11, MYO5A, MYOD1, NAB2, NBN, NCKIPSD, NCOA1, NCOA2, NCOA3, NCOA4, NCOR1, NCOR2, NDRG1, NF1, NF2, NFATC2, NFE2L2, NFIB, NFIX, NFKB2, NFKBIA, NFKBIE, NID1, NIN, NKX2-1, NKX3-1, NONO, NOTCH1, NOTCH2, NOTCH3, NPM1, NR4A3, NRAS, NRG1, NSD1, NSD2, NSD3, NT5C2, NTHL1, NTN1, NTRK1, NTRK2, NTRK3, NUP214, NUP93, NUP98, NUTM1, NUTM2B, OMD, P2RY8, PAG1, PAK1, PALB2, PARG, PARP1, PATZ1, PAX3, PAX5, PAX7, PAX8, PBRM1, PBX1, PCBP1, PCM1, PDCD1, PDCD1LG2, PDE4DIP, PDGFB, PDGFRA, PDGFRB, PDK1, PDS5B, PGR, PHF1, PHF6, PHOX2B, PICALM, PIGA, PIK3CA, PIK3CB, PIK3CD, PIK3CG, PIK3R1, PIK3R2, PIM1, PLA2R1, PLAG1, PLCG1, PLCG2, PMEL, PML, PMS1, PMS2, POLD1, POLE, POT1, POU5F1, PPARG, PPFIBP1, PPM1D, PPP2R1A, PPP2R2A, PPP6C, PRCC, PRDM1, PRDM16, PRDM9, PREX2, PRKACA, PRKAR1A, PRKCE, PRKCI, PRKN, PSIP1, PTCH1, PTEN, PTK2B, PTPN11, PTPN13, PTPRD, PTPRK, PTPRT, PVT1, PWWP2A, QKI, RABEP1, RAC1, RAD21, RAD50, RAD51, RAD51B, RAD51C, RAD51D, RAD54L, RAF1, RANBP2, RAP1GDS1, RARA, RARG, RASA1, RB1, RBL2, RBM10, RBM14, RBM15, RECQL, RECQL4, REL, RET, RHEB, RHOA, RICTOR, RIF1, RIT1, RNF43, ROS1, RPL11, RPL5, RPN1, RPS15, RPS27, RPS6KB1, RPTOR, RRAS2, RSPO2, RSPO3, RTEL1, RUNX1, RUNX1T1, SAMD9L, SAMHD1, SAV1, SBDS, SCG5, SDC4, SDHA, SDHAF2, SDHB, SDHC, SDHD, SEC31A, SEMA6A, SEPTIN5, SEPTIN6, SEPTIN9, SET, SETBP1, SETD2, SETDB1, SF3B1, SFPQ, SFTPA1, SFTPB, SFTPC, SGK1, SH2B3, SH3GL1, SHLD1, SHLD2, SHLD3, SLC34A2, SLC45A3, SLFN11, SLX4, SMAD2, SMAD3, SMAD4, SMARCA1, SMARCA2, SMARCA4, SMARCA5, SMARCB1, SMARCD1, SMARCE1, SMC1A, SMC3, SMO, SMURF1, SND1, SOCS1, SOS1, SOX2, SOX4, SOX9, SPEN, SPI1, SPOP, SPRED1, SPTBN1, SQSTM1, SRC, SRGAP3, SRSF2, SS18, SS18L1, SSBP2, SSX1, SSX2, SSX4, STAG1, STAG2, STAT1, STAT3, STAT5A, STAT5B, STAT6, STIL, STK11, STRN, SUFU, SUZ12, SYK, TACC1, TAF15, TAL1, TAL2, TAP2, TBL1XR1, TBX3, TCEA1, TCF12, TCF3, TCF7L2, TCL1A, TCL1B, TENT5C, TERC, TERT, TET1, TET2, TFCP2, TFE3, TFEB, TFG, TGFBR1, TGFBR2, THRAP3, TLX1, TLX3, TMEM127, TMPRSS2, TNFAIP3, TNFRSF14, TNKS2, TOP1, TP53, TP53BP1, TP63, TPM3, TPM4, TPR, TRA, TRAF2, TRAF3, TRAF7, TRB, TRIM24, TRIM27, TRIM33, TRIM37, TSC1, TSC2, TSHR, TSLP, TYK2, U2AF1, U2AF2, UBTF, USP2, USP6, USP9X, VEGFA, VHL, VMP1, VTCN1, VTI1A, WDR48, WIF1, WRN, WT1, WWTR1, XPA, XPC, XPO1, XRCC2, YAP1, YWHAE, YY1, ZBTB16, ZBTB7A, ZCCHC7, ZCCHC8, ZEB1, ZFHX3, ZFP36L1, ZFP36L2, ZFYVE19, ZMYM2, ZNF217, ZNF384, ZNF692, ZNF703, ZNF750, ZNRF3, ZRSR2.


### Variant Annotation sources
The annotated variants are based on The Personal Cancer Genome Reporter (PCGR) [PMID: 29272339](https://pubmed.ncbi.nlm.nih.gov/29272339/)

 PCGR draws variant data from the following underlying tools and knowledge resources:
- [VEP](http://asia.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor (v101)
- [GENCODE](https://www.gencodegenes.org/) - reference gene annotation and experimental validation (release 35/19)
- [dbNSFP](https://pubmed.ncbi.nlm.nih.gov/23843252/) - Database of non-synonymous functional predictions (20200616 (v4.1))
- [Pfam](https://www.ebi.ac.uk/interpro/) - Collection of protein families/domains (2020_05 (v33.1))
- [TCGA](https://www.genome.gov/Funded-Programs-Projects/Cancer-Genome-Atlas) - The Cancer Genome Atlas - somatic mutations (20200722 (v25))
- [ICGC-PCAWG](https://docs.icgc.org/pcawg/) - ICGC-Pancancer Analysis of Whole Genomes - somatic mutations (2020_01)
- [TCGA-PCDM](https://www.cancer.gov/ccg/research/genome-sequencing/tcga) - Putative Cancer Driver Mutations on multiple discovery approaches (2019)
- [UniProtKB](https://www.uniprot.org/) - Resource of protein sequence and functional information (2020_04)
- [CORUM](http://mips.helmholtz-muenchen.de/corum/) - Resource of mammalian protein complexes (20180903 (release 3.0))
- [gnomAD](https://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (r2.1 (October 2018))
- [COSMIC](https://cancer.sanger.ac.uk/cosmic) - Catalogue of somatic mutations in cancer (90)
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) - Database of short genetic variants (153)
- [1000Genomes](https://www.internationalgenome.org/) - Germline variant frequencies genome-wide (20130502 (phase 3))
- [DisGenet](https://www.disgenet.org/) - Database of gene-disease associations (2020_05 (v7.0))
- [DoCM](http://www.docm.info/) - Database of curated mutations (release 3.2)
- [CancerHotspots](http://www.cancerhotspots.org/#/home) - A resource for statistically significant mutations in cancer (2017)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - Database of genomic variants of clinical significance (20200803)
- [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - Literature database tumor suppressor/oncogenes (20200908 (release 28))
- [OncoTree](http://oncotree.mskcc.org/) - Ontology for standardization of cancer type diagnosis (2020_04_01)
- [DiseaseOntology](https://disease-ontology.org/) - Standardized ontology for human disease (20200915)
- [EFO](https://github.com/EBISPOT/efo) - Experimental Factor Ontology (v3.22)
- [OncoScore](https://bioconductor.org/packages/release/bioc/html/OncoScore.html) - Literature-based ranking of gene-cancer association (20191227)
- [OpenTargetPlatform](https://platform.opentargets.org/) - Data integration for access to drug targets (2020_06)
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - Manually curated database of bioactive molecules (20200518 (v27))
- [KEGG](https://www.genome.jp/kegg/pathway.html) - Knowledgebase molecular interaction, reaction and relation networks (20200910)
- [CIViC](https://civicdb.org/) - Clinical interpretations of variants in cancer (20200920)
- [CGI](https://cancergenomeinterpreter.org/biomarkers) - Cancer Genome Interpreter Cancer Biomarkers Database (20180117)

The report is generated with PCGR version 1.4.1.

### Tumor mutational burden (TMB) estimate
The size of the targeted (coding) genomic region has been defined as: 34 Mb. For estimation of TMB, we count all coding, somatic base substitutions and indels in the targeted regions, including synonymous alterations, as outlined previously [PMID: 28420421](https://pubmed.ncbi.nlm.nih.gov/28420421/). TMB-low is defined as <10 mutations/Mb; TMB-high is defined as ≥10 mutations/Mb; for cases with TMB >50 mutations/Mb, TMB will be reported as “TMB-High (50 mutations/Megabase)”.

### Microsatellite instability (MSI) detection
MSI was estimated using the hmftools [PURPLE estimator](https://github.com/hartwigmedical/hmftools/tree/master/purple). MSI is classified as either MSI-Unstable or MSI-Stable. 

### Somatic Signature detection
COSMIC V3.2 base substitution signatures are assigned signatures using [MutationalPatterns v3.8.0](https://cancer.sanger.ac.uk/signatures/; [PMID: 32025018](https://pubmed.ncbi.nlm.nih.gov/32025018/).

### Homologous recombination deficiency detection
Homologous recombination deficiency (HRD) was detected by DRAGEN HRD, HRDetect (signature.tools.lib 2.1.2) [PMID: 28288110](https://pubmed.ncbi.nlm.nih.gov/28288110/) and CHORD [PMID: 33149131](https://pubmed.ncbi.nlm.nih.gov/33149131/). DRAGEN HRD is composed of the sum of three components: loss of heterozygosity, telomeric allele imbalance, and large-scale state transition. HRDetect and CHORD consider mutational patterns (SNVs, InDels and SVs) that are characteristic of HR-deficient tumours. DRAGEN HRD scores of ≥ 42 are categorised as HR-deficient. HRDetect grants a score from 0 to 1; tumour samples with a score > 0.7 are categorized as HR-deficient. CHORD classifies tumours in BRCA1-deficient and BRCA2-deficient categories; tumours with a combined probability of these HRD categories > 50%, are categorised as HR-deficient.

### Whole transcriptome sequencing (WTS)
Whole transcriptome sequencing was performed on RNA extracted from the same fresh or FFPE tissue used for WGS (dual DNA/RNA extraction). Approximately 100 million reads were generated per RNA sequencing library. WTS data was aligned and quantified using the Illumina [DRAGEN RNA pipeline (v4.4)](https://help.dragen.illumina.com/). Fusions were confirmed in RNA-seq data using [Arriba v2.5](https://github.com/oicr-gsi/arriba) and the DRAGEN fusion caller. Tertiary analysis of WTS data was performed using [RNAsum v2](https://github.com/umccr/RNAsum), an R package that can post-process, summarize and visualize WTS data outputs. 

## Disclaimers
**The analysis assumes sample identification, family relationships, and clinical diagnoses are as provided. The evidence compiled in the report is derived from third-party data sources and should be interpreted in the context of all other information for this patient. While every effort is made to ensure accuracy, statements regarding clinical trials and therapeutic recommendations made on this report should not be considered exhaustive, or applicable in all clinical situation. For further information, please contact the laboratory.**

Please contact the laboratory on UMCCR-Genomics@unimelb.edu.au if you wish to discuss this report further.
