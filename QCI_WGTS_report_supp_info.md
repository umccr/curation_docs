
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
CNV/SV changes are annotated with the svprioritize framework, assigning priority to fusion events, whole exon loss or upstream/downstream changes for a list of 1246 cancer-associated genes only. The genes reviewed for CNV and SVs:  A1CF,  ABCB1, ABCB11,  ABI1,  ABL1, ABL2, ABRAXAS1, ACD, ACKR3, ACSL3, ACSL6, ACTB,  ACTG1, ACVR1, ACVR1B, ACVR2A, ADGRA2, AFDN, AFF1, AFF3, AFF4, AGO2, AIP, AJUBA, AKAP9, AKT1, AKT2, AKT3, ALDH2, ALK, ALOX12B, ALOX15B, AMER1, ANK1, ANKRD11, ANKRD26, APC, APOB, APOBEC3B, AR, ARAF, ARFRP1, ARHGAP26, ARHGAP35, ARHGAP5, ARHGEF10, ARHGEF10L, ARHGEF12, ARHGEF28, ARID1A, ARID1B, ARID2, ARID3A, ARID3B, ARID3C, ARID4A, ARID4B, ARID5A, ARID5B, ARNT, ARTN, ASPSCR1, ASXL1, ASXL2, ATF1, ATF2, ATF3, ATF7IP, ATIC, ATM, ATP1A1, ATP2B3, ATP6AP1, ATP6V1B2, ATR, ATRX, ATXN2, AURKA, AURKB, AXIN1, AXIN2, AXL, B2M, BABAM1, BACH2, BAP1, BARD1, BAX, BAZ1A, BBC3, BCL10, BCL11A, BCL11B, BCL2, BCL2L1, BCL2L11, BCL2L12, BCL2L2, BCL3, BCL6, BCL6B, BCL7A, BCL9, BCL9L, BCLAF1, BCOR, BCORL1, BCR, BECN1, BIRC3, BIRC6, BLM, BMI1, BMP5, BMPR1A, BRAF, BRCA1, BRCA2, BRD3, BRD4, BRINP1, BRIP1, BTG1, BTG2, BTK, BUB1B, C15orf65, C2orf40, C2orf44, CACNA1D, CAD, CADM1, CALR, CAMTA1, CANT1, CARD11, CARM1, CARS, CASC5, CASP3, CASP8, CASP9, CASR, CBFA2T3, CBFB, CBL, CBLB, CBLC, CCAT1, CCDC6, CCNB1IP1, CCNC, CCND1, CCND2, CCND3, CCNE1, CCR4, CCR7, CD19, CD209, CD22, CD274, CD276, CD28, CD58, CD70, CD74, CD79A, CD79B, CDC42, CDC73, CDH1, CDH10, CDH11, CDH17, CDK12, CDK16, CDK4, CDK6, CDK8, CDKN1A, CDKN1B, CDKN1C, CDKN2A, CDKN2B, CDKN2C, CDX2, CEBPA, CENPA, CEP57, CEP89, CHCHD7, CHD2, CHD4, CHD5, CHEK1, CHEK2, CHIC2, CHST11, CIC, CIITA, CKS1B, CLIP1, CLP1, CLTC, CLTCL1, CNBD1, CNBP, CNOT3, CNTNAP2, CNTRL, COL1A1, COL2A1, COL3A1, COL7A1, COX6C, CPEB3, CRBN, CREB1, CREB3L1, CREB3L2, CREBBP, CRKL, CRLF2, CRNDE, CRNKL1, CRTC1, CRTC3, CSDE1, CSF1R, CSF3R, CSMD3, CTCF, CTLA4, CTNNA1, CTNNA2, CTNNB1, CTNND1, CTNND2, CTR9, CTRC, CUL3, CUX1, CXCR4, CYLD, CYP2C8, CYP2D6, CYR61, CYSLTR2, DAB2, DAPK1, DAXX, DCAF12L2, DCC, DCTN1, DCUN1D1, DDB2, DDIT3, DDR2, DDX10, DDX3X, DDX5, DDX6, DDX41, DEK, DGCR8, DICER1, DIRAS3, DIS3, DIS3L2, DKC1, DLC1, DNAJB1, DNM2, DNMT1, DNMT3A, DNMT3B, DOCK8, DOT1L, DPYD, DROSHA, DTX1, DTX3L, DUSP22, DUSP4, DUX4L1, E2F3, EBF1, ECT2L, EED, EGF, EGFL7, EGFR, EGR1, EIF1AX, EIF3E, EIF4A2, EIF4E, ELANE, ELF3, ELF4, ELK4, ELL, ELN, ELOC, EML4, EMP3, EMSY, ENG, EP300, EP400, EPAS1, EPCAM, EPHA2, EPHA3, EPHA5, EPHA7, EPHB1, EPHB2, EPOR, EPS15, ERBB2, ERBB3, ERBB4, ERC1, ERCC1, ERCC2, ERCC3, ERCC4, ERCC5, ERF, ERG, ERRFI1, ESCO2, ESR1, ETNK1, ETS1, ETV1, ETV4, ETV5, ETV6, EWSR1, EXT1, EXT2, EZH1, EZH2, EZR, FAH, FAM131B, FAM135B, FAM175A, FAM46C, FAM47C, FAM58A, FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, FAS, FASLG, FASN, FAT1, FAT3, FAT4, FBLN2, FBXO11, FBXW7, FCGR2B, FCRL4, FEN1, FES, FEV, FGF1, FGF10, FGF14, FGF19, FGF2, FGF23, FGF3, FGF4, FGF5, FGF6, FGF7, FGF8, FGF9, FGFR1, FGFR1OP, FGFR2, FGFR3, FGFR4, FH, FHIT, FIP1L1, FKBP9, FLCN, FLG, FLI1, FLNA, FLT1, FLT3, FLT4, FNBP1, FOS, FOXA1, FOXA2, FOXL2, FOXM1, FOXO1, FOXO3, FOXO4, FOXP1, FOXP3, FOXQ1, FOXR1, FRS2, FSTL3, FTO, FUBP1, FURIN, FUS, FYN, GAB2, GABRA6, GALNT12, GAS7, GATA1, GATA2, GATA3, GATA4, GATA6, GBA, GEN1, GID4, GJB2, GKN1, GLI1, GMPS, GNA11, GNA12, GNA13, GNAQ, GNAS, GNB1, GNMT, GOLGA5, GOPC, GPC3, GPC5, GPHN, GPRC5A, GPS2, GREM1, GRIN2A, GRM3, GSK3B, GSTP1, GTF2I, GZMA, H3F3A, H3F3B, H3F3C, HABP2, HDAC1, HDAC4, HDAC7, HERPUD1, HEY1, HFE, HGF, HIF1A, HIP1, HIST1H1B, HIST1H1C, HIST1H1D, HIST1H1E, HIST1H2AC, HIST1H2AG, HIST1H2AL, HIST1H2AM, HIST1H2BC, HIST1H2BD, HIST1H2BG, HIST1H2BJ, HIST1H2BK, HIST1H2BO, HIST1H3A, HIST1H3B, HIST1H3C, HIST1H3D, HIST1H3E, HIST1H3F, HIST1H3G, HIST1H3H, HIST1H3I, HIST1H3J, HIST1H4I, HIST2H3A, HIST2H3C, HIST2H3D, HIST3H3, HLA-A, HLA-B, HLA-C, HLF, HMBS, HMGA1, HMGA2, HNF1A, HNF1B, HNRNPA2B1, HOOK3, HOXA11, HOXA13, HOXA9, HOXB13, HOXC11, HOXC13, HOXD11, HOXD13, HPGD, HRAS, HSD3B1, HSP90AA1, HSP90AB1, ICOSLG, ID3, ID4, IDH1, IDH2, IFNGR1, IGF1, IGF1R, IGF2, IGF2BP2, IGFBP2, IGFBP7, IGH, IGK, IGL, IKBKB, IKBKE, IKZF1, IKZF3, IL10, IL2, IL21R, IL2RA, IL6, IL6R, IL6ST, IL7R, INHA, INHBA, INPP4A, INPP4B, INPPL1, INSR, IRF1, IRF2, IRF4, IRF8, IRS1, IRS2, IRS4, ISX, ITGAV, ITK, JAK1, JAK2, JAK3, JARID2, JAZF1, JMJD1C, JUN, KAT6A, KAT6B, KAT7, KCNJ5, KDM5A, KDM5C, KDM6A, KDR, KDSR, KEAP1, KEL, KIAA1549, KIAA1598, KIF14, KIF1B, KIF5B, KISS1, KIT, KL, KLF4, KLF5, KLF6, KLHL6, KLK2, KLLN, KMT2A, KMT2B, KMT2C, KMT2D, KMT5A, KNL1, KNSTRN, KRAS, KSR2, KTN1, LAMP1, LARP4B, LASP1, LATS1, LATS2, LCK, LCP1, LEF1, LEPROTL1, LHFP, LHFPL6, LIFR, LIN28B, LMNA, LMO1, LMO2, LOX, LPP, LRIG3, LRP1B, LRRC4, LRRK2, LSM14A, LTB, LYL1, LYN, LZTR1, MACC1, MAF, MAFB, MAGI2, MALAT1, MALT1, MAML2, MAP2K1, MAP2K2, MAP2K4, MAP3K1, MAP3K13, MAP3K14, MAP3K6, MAP3K7, MAPK1, MAPK3, MAPKAP1, MARK2, MAX, MB21D2, MCC, MCL1, MDC1, MDM2, MDM4, MDS2, MECOM, MED12, MEF2B, MEIS1, MEN1, MET, MGA, MGAM, MGMT, MIB1, MITF, MKI67, MKL1, MLF1, MLH1, MLH3, MLLT1, MLLT10, MLLT11, MLLT3, MLLT4, MLLT6, MN1, MNX1, MOB3B, MPEG1, MPL, MRE11, MRE11A, MS4A1, MSH2, MSH3, MSH6, MSI1, MSI2, MSN, MST1, MST1R, MTAP, MTCP1, MTDH, MTOR, MTSS1, MUC1, MUC16, MUC4, MUTYH, MYB, MYC, MYCL, MYCN, MYD88, MYH11, MYH9, MYO18B, MYO5A, MYOD1, N4BP2, NAB2, NACA, NAV3, NBEA, NBN, NCKIPSD, NCOA1, NCOA2, NCOA3, NCOA4, NCOR1, NCOR2, NCSTN, NDRG1, NDRG2, NEGR1, NEK1, NF1, NF2, NFATC2, NFE2, NFE2L2, NFE2L3, NFIB, NFKB1, NFKB2, NFKBIA, NFKBIE, NHP2, NIN, NKX2-1, NKX3-1, NONO, NOP10, NOTCH1, NOTCH2, NOTCH3, NOTCH4, NPM1, NR0B2, NR4A3, NRAS, NRG1, NSD1, NSD2, NSD3, NT5C2, NTHL1, NTRK1, NTRK2, NTRK3, NUF2, NUMA1, NUP214, NUP93, NUP98, NUTM1, NUTM2A, NUTM2B, OGG1, OLIG2, OMD, OTX2, P2RY8, PABPC1, PAFAH1B2, PAK1, PAK3, PAK5, PAK7, PALB2, PARK2, PARK7, PARP1, PATZ1, PAWR, PAX3, PAX5, PAX6, PAX7, PAX8, PBRM1, PBX1, PCBP1, PCLO, PCM1, PDCD1, PDCD1LG2, PDCD4, PDE4DIP, PDGFB, PDGFRA, PDGFRB, PDK1, PDL1, PDL2, PDLIM7, PDPK1, PDS5B, PDX1, PER1, PGR, PHF1, PHF6, PHOX2B, PI3, PICALM, PIGA, PIK3C2B, PIK3C2G, PIK3C3, PIK3CA, PIK3CB, PIK3CD, PIK3CG, PIK3R1, PIK3R2, PIK3R3, PIM1, PINK1, PLAG1, PLAGL2, PLCG1, PLCG2, PLK2, PLXNB1, PMAIP1, PML, PMS1, PMS2, PNRC1, POLD1, POLE, POLG, POLH, POLQ, POT1, POU2AF1, POU2F2, POU3F2, POU5F1, PPARG, PPFIBP1, PPM1D, PPP2R1A, PPP2R2A, PPP4R2, PPP6C, PRCC, PRDM1, PRDM14, PRDM16, PRDM2, PRDM9, PREX2, PRF1, PRKACA, PRKAR1A, PRKCB, PRKCI, PRKD1, PRKDC, PRKN, PRPF40B, PRR11, PRRX1, PRSS1, PRSS8, PSIP1, PTCH1, PTEN, PTGDR, PTGS2, PTK6, PTP4A1, PTPN1, PTPN11, PTPN13, PTPN2, PTPN6, PTPRB, PTPRC, PTPRD, PTPRK, PTPRO, PTPRS, PTPRT, PTTG1, PWWP2A, QKI, RAB25, RAB35, RABEP1, RAC1, RAC2, RAD17, RAD21, RAD50, RAD51, RAD51B, RAD51C, RAD51D, RAD52, RAD54L, RAF1, RALGDS, RANBP2, RAP1GAP, RAP1GDS1, RARA, RARRES1, RASA1, RASAL1, RASSF1, RB1, RBM10, RBM15, RCC2, RECQL, RECQL4, REL, RELB, RELN, RET, RFWD2, RFWD3, RGPD3, RGS7, RHBDF2, RHEB, RHNO1, RHOA, RHOH, RICTOR, RING1, RINT1, RIT1, RMI2, RMRP, RNF213, RNF43, ROBO1, ROBO2, ROS1, RPL10, RPL22, RPL5, RPN1, RPS15, RPS6KA4, RPS6KB1, RPS6KB2, RPTOR, RRAGC, RRAS, RRAS2, RSPO2, RSPO3, RTEL1, RUNX1, RUNX1T1, RUNX3, RXRA, RYBP, S100A7, SALL4, SAMHD1, SASH1, SBDS, SCG5, SDC4, SDHA, SDHAF2, SDHB, SDHC, SDHD, SEC31A, SEMA3B, SEPT5, SEPT6, SEPT9, SERPINA1, SERPINB5, SESN1, SESN2, SESN3, SET, SETBP1, SETD1A, SETD1B, SETD2, SETD3, SETD4, SETD5, SETD6, SETD7, SETDB1, SETDB2, SF3B1, SFPQ, SFRP1, SFRP4, SGK1, SH2B3, SH2D1A, SH3GL1, SHH, SHOC2, SHQ1, SHTN1, SIN3A, SIRPA, SIRT1, SIRT2, SIRT3, SIRT6, SIX1, SIX2, SKI, SKP2, SLC25A13, SLC34A2, SLC45A3, SLIT2, SLX4, SMAD2, SMAD3, SMAD4, SMARCA1, SMARCA4, SMARCB1, SMARCD1, SMARCE1, SMC1A, SMC3, SMG1, SMO, SMURF1, SMYD3, SNCAIP, SND1, SNX29, SOCS1, SOS1, SOX1, SOX10, SOX17, SOX2, SOX21, SOX4, SOX9, SP140, SPARC, SPDEF, SPECC1, SPEN, SPI1, SPINK1, SPOP, SPRED1, SPRTN, SPTA1, SRC, SRF, SRGAP3, SRSF2, SRSF3, SRY, SS18, SS18L1, SSX1, SSX2, SSX4, STAG1, STAG2, STAT1, STAT3, STAT4, STAT5A, STAT5B, STAT6, STIL, STK11, STK19, STK40, STRN, SUFU, SUZ12, SYK, TACC1, TAF1, TAF15, TAL1, TAL2, TAP1, TAP2, TAZ, TBL1XR1, TBX3, TCEA1, TCEB1, TCF12, TCF3, TCF7L2, TCL1A, TEC, TEK, TENT5C, TERC, TERF2IP, TERT, TET1, TET2, TET3, TFCP2, TFE3, TFEB, TFF1, TFG, TFPT, TFRC, TGFBR1, TGFBR2, THRAP3, THY1, TLR4, TLX1, TLX3, TMEM127, TMPRSS2, TNC, TNFAIP3, TNFRSF14, TNFRSF17, TNFRSF6, TOP1, TOP2A, TP53, TP53BP1, TP63, TPD52, TPM3, TPM4, TPMT, TPR, TRA, TRAF2, TRAF3, TRAF5, TRAF7, TRB, TRD, TRIM24, TRIM27, TRIM3, TRIM33, TRIM37, TRIP11, TRRAP, TSC1, TSC2, TSHR, TTF1, TUSC3, TWIST1, TWSG1, TYK2, TYMS, U2AF1, U2AF2, UBR5, UCHL1, UGT1A1, UHRF1, UPF1, UROD, USP44, USP6, USP8, USP9X, VAV1, VAV2, VEGFA, VHL, VTCN1, VTI1A, WAS, WDCP, WEE1, WHSC1, WHSC1L1, WIF1, WISP1, WISP3, WNK2, WNT1, WNT5A, WRN, WT1, WWOX, WWTR1, XBP1, XIAP, XPA, XPC, XPO1, XRCC2, YAP1, YES1, YWHAE, ZBTB16, ZBTB2, ZBTB7A, ZCCHC8, ZEB1, ZFHX3, ZMYM2, ZMYM3, ZNF198, ZNF217, ZNF278, ZNF331, ZNF384, ZNF429, ZNF479, ZNF521, ZNF703, ZNF750, ZNRF3, ZRSR2.


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
