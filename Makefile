.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

R_SCRIPT=Rscript
SCRIPTS_DIR=scripts/


# 1) From raw data
DATA/$(species)/$(tissue)/$(tissue).txt:
	$(R_SCRIPT) $(SCRIPTS_DIR)GeneID_Transform.R $(species) $(tissue)
##1)GeneIDTransform: From raw_tissue.txt file, transform original_IDs to GeneIDs (not unique) into tissue.txt file.
.PHONY : GeneIDTransform
GeneIDTransform : DATA/$(species)/$(tissue)/$(tissue).txt


##2)Running algorithms:
##RAIN: run RAIN
DATA/$(species)/$(tissue)/RAIN.txt: DATA/$(species)/$(tissue)/$(tissue).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)RAIN_script.R $(species) $(tissue)
.PHONY : RAIN
RAIN : DATA/$(species)/$(tissue)/RAIN.txt

##empJTK: run empJTK
DATA/$(species)/$(tissue)/empJTK.txt: DATA/$(species)/$(tissue)/$(tissue).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)empJTK_script.R $(species) $(tissue)
.PHONY : empJTK
empJTK : DATA/$(species)/$(tissue)/empJTK.txt

##MetaCycle: run MetaCycle
DATA/$(species)/$(tissue)/JTK.txt: DATA/$(species)/$(tissue)/$(tissue).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)MetaCycle_script.R $(species) $(tissue)
.PHONY : MetaCycle
MetaCycle : DATA/$(species)/$(tissue)/JTK.txt

##GeneCycle: run GeneCycle
DATA/$(species)/$(tissue)/GeneCycle.txt: DATA/$(species)/$(tissue)/$(tissue).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)GeneCycle_script.R $(species) $(tissue)
.PHONY : GeneCycle
GeneCycle : DATA/$(species)/$(tissue)/GeneCycle.txt


# 3) Normalization of p-values per gene (ProbIDs/TranscriptsIDs -> uniqueGeneID) : Normalize p-values produced by rhythm detection algorithms
##3)pvaluePerGeneNormalization: Normalization of p-values per gene (ProbIDs/TranscriptsIDs -> uniqueGeneID) : For each gene with several data (ProbIDs or TranscriptsIDs), normalize p-values produced by rhythm detection algorithms. Normalization by mean/fisher/kost/brown
.PHONY : pvaluePerGeneNormalization
pvaluePerGeneNormalization : DATA/$(species)/$(tissue)/normalized_BH.Q/normalized_pvalue_per_gene.txt
	$(R_SCRIPT) $(SCRIPTS_DIR)normalization_per_gene.R $(species) $(tissue)


# 4) Generate the plots of density distribution of p-values
##4)pvaluesDistributionAnalysis: Plot the density distribution of p-values
.PHONY : pvaluesDistributionAnalysis
pvaluesDistributionAnalysis : DATA/$(species)/$(tissue)/normalized_BH.Q/normalized_pvalue_per_gene.txt
	$(R_SCRIPT) $(SCRIPTS_DIR)pvalues_distrib_analysis.R $(species) $(tissue)


# 5) UpSet plot of rhythmic genes called in common:
##5)UpsetDiagramRhythmicGenes: Upset diagram of genes detected rhythmic by the different algorithms. p-values can be adjusted by changing p_adj_method argument. By default p_adj_method=none. Other choices are: holm, hochberg, hommel, bonferroni, BH, BY, fdr
threshold_1=0.05
threshold_2=0.01
p_adj_method=none
.PHONY : UpsetDiagramRhythmicGenes
UpsetDiagramRhythmicGenes : DATA/$(species)/$(tissue)/normalized_default.pvalue/normalized_pvalue_per_gene.txt
	$(R_SCRIPT) $(SCRIPTS_DIR)upset_rhythmicgenes.R $(species) $(tissue) $(threshold_1) $(threshold_2) $(p_adj_method)


# 6) Generate files with orthology informations:
OMA_DIR=DATA/ORTHOLOGY/oma/
ORTHO_DIR=DATA/ORTHOLOGY/result/

##6a)UpdateOMA: update recover OMA orthologs pairs between 2 species
.PHONY : UpdateOMA
UpdateOMA :
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=MOUSE&p2=PAPAN&p3=EnsemblGene" >> $(OMA_DIR)mouse_baboon_tmp.txt
	{ echo "mouse_ID\tbaboon_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)mouse_baboon_tmp.txt; } > $(OMA_DIR)mouse_baboon_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=MOUSE&p2=RATNO&p3=EnsemblGene" >> $(OMA_DIR)mouse_rat_tmp.txt ;
	{ echo "mouse_ID\trat_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)mouse_rat_tmp.txt; } > $(OMA_DIR)mouse_rat_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=RATNO&p2=PAPAN&p3=EnsemblGene" >> $(OMA_DIR)rat_baboon_tmp.txt ;
	{ echo "rat_ID\tbaboon_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)rat_baboon_tmp.txt; } > $(OMA_DIR)rat_baboon_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=MOUSE&p2=DANRE&p3=EnsemblGene" >> $(OMA_DIR)mouse_zebrafish_tmp.txt ;
	{ echo "mouse_ID\tzebrafish_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)mouse_zebrafish_tmp.txt; } > $(OMA_DIR)mouse_zebrafish_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=DANRE&p2=PAPAN&p3=EnsemblGene" >> $(OMA_DIR)zebrafish_baboon_tmp.txt ;
	{ echo "zebrafish_ID\tbaboon_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)zebrafish_baboon_tmp.txt; } > $(OMA_DIR)zebrafish_baboon_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=DROME&p2=AEDAE&p3=Source" >> $(OMA_DIR)drosophila_aedes_tmp.txt ;
	{ echo "drosophila_ID\taedes_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)drosophila_aedes_tmp.txt; } > $(OMA_DIR)drosophila_aedes_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=DROME&p2=ANOGA&p3=Source" >> $(OMA_DIR)drosophila_anopheles_tmp.txt ;
	{ echo "drosophila_ID\tanopheles_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)drosophila_anopheles_tmp.txt; } > $(OMA_DIR)drosophila_anopheles_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=AEDAE&p2=ANOGA&p3=Source" >> $(OMA_DIR)aedes_anopheles_tmp.txt ;
	{ echo "aedes_ID\tanopheles_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)aedes_anopheles_tmp.txt; } > $(OMA_DIR)aedes_anopheles_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=MOUSE&p2=DROME&p3=Source" >> $(OMA_DIR)mouse_drosophila_tmp.txt ;
	{ echo "mouse_ID\tdrosophila_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)mouse_drosophila_tmp.txt; } > $(OMA_DIR)mouse_drosophila_oma.txt ;
	curl "https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=PAPAN&p2=HUMAN&p3=Source" >> $(OMA_DIR)baboon_human_tmp.txt ;
	{ echo "baboon_ID\thuman_ID\torthotype\tOMAgroup"; cat $(OMA_DIR)baboon_human_tmp.txt; } > $(OMA_DIR)baboon_human_oma.txt ;
	rm $(OMA_DIR)*tmp.txt


# by default, same homologous tissue is compared :
tissue1=$(tissue)
tissue2=$(tissue)

##6b)OrthologyPerGene: Retrieve all data of species1 with orthology data of species2. When one gene of species 1 has several orthologs in species 2, P-values of orthologs of species 2 are normalized by fisher or sidak method. See orthology_per_gene.csv file.

DATA/$(species1)/$(tissue1)/normalized_default.pvalue/normalized_pvalue_per_gene.txt:
	$(R_SCRIPT) $(SCRIPTS_DIR)normalization_per_gene.R $(species1) $(tissue1)

DATA/$(species2)/$(tissue2)/normalized_default.pvalue/normalized_pvalue_per_gene.txt:
	$(R_SCRIPT) $(SCRIPTS_DIR)normalization_per_gene.R $(species2) $(tissue2)

$(ORTHO_DIR)default.pvalue_$(species1)_$(tissue1)_$(species2)_$(tissue2).txt: DATA/$(species1)/$(tissue1)/normalized_default.pvalue/normalized_pvalue_per_gene.txt DATA/$(species2)/$(tissue2)/normalized_default.pvalue/normalized_pvalue_per_gene.txt
	$(R_SCRIPT) $(SCRIPTS_DIR)orthology_per_gene.R $(species1) $(tissue1) $(species2) $(tissue2)

$(ORTHO_DIR)default.pvalue_$(species2)_$(tissue2)_$(species1)_$(tissue1).txt: DATA/$(species1)/$(tissue1)/normalized_default.pvalue/normalized_pvalue_per_gene.txt DATA/$(species2)/$(tissue2)/normalized_default.pvalue/normalized_pvalue_per_gene.txt
	$(R_SCRIPT) $(SCRIPTS_DIR)orthology_per_gene.R $(species2) $(tissue2) $(species1) $(tissue1)

.PHONY : OrthologyPerGene
	OrthologyPerGene : $(ORTHO_DIR)default.pvalue_$(species1)_$(tissue1)_$(species2)_$(tissue2).txt $(ORTHO_DIR)default.pvalue_$(species2)_$(tissue2)_$(species1)_$(tissue1).txt


# 7) p-values density distribution plots of rhythmic orthologs vs non-rhythmic orthologs
##7)pvaluesDistribOrthologs: Plot the p-values density distribution of rhythmic orthologs vs non-rhythmic orthologs (See method). argument P_VAL should be among: raw.pvalue, default.pvalue, or BH.Q
.PHONY : pvaluesDistribOrthologs
pvaluesDistribOrthologs : $(ORTHO_DIR)default.pvalue_$(species1)_$(tissue1)_$(species2)_$(tissue2).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)pvalues_distrib_ORTHO.R $(species1) $(tissue1) $(species2) $(tissue2) $(species2.method) $(species2.threshold) $(ONLY_WITHIN_CONSERVED_GENES) $(ONLY_ORTHO_ONE_TO_ONE) $(P_VAL)


# 8) Variation of the proportion of rhythmic orthologs among all species1-species2 orthologs as a function of the number of orthologs detected rhythmic in species1.
##8)ProportionRhythmicOrthoPlots: Plot the variation of the proportion of rhythmic orthologs among all species1-species2 orthologs as a function of the number of orthologs detected rhythmic in species1.
.PHONY : ProportionRhythmicOrthoPlots
ProportionRhythmicOrthoPlots : #$(ORTHO_DIR)default.pvalue_$(species1)_$(tissue1)_$(species2)_$(tissue2).txt
	$(R_SCRIPT) $(SCRIPTS_DIR)rhythmic_orthologs_proportions.R $(species1) $(tissue1) $(species2) $(tissue2) $(species2.method) $(species2.threshold) $(species2.pvalue) $(point.species1.threshold) $(ONLY_WITHIN_CONSERVED_GENES) $(ONLY_ORTHO_ONE_TO_ONE)





##9)Other functions:

#.PHONY : orderedgenesUpsetDiagramRhythmicGenes
#orderedgenesUpsetDiagramRhythmicGenes : DATA/$(species)/$(tissue)/normalized_default.pvalue/normalized_pvalue_per_gene.txt
#	$(R_SCRIPT) $(SCRIPTS_DIR)upset_rhythmicgenes_ordered.R $(species) $(tissue) $(threshold)

##°circadGenesRunning: Run algorithms and normalization on known circadian genes list
.PHONY : circadGenesRunning
circadGenesRunning :
	$(R_SCRIPT) $(SCRIPTS_DIR)algo_and_normalization_through_known_circad_genes.R $(species) $(tissue)

##°restrictedTimepointsRunning: Run algorithms and normalization on mouse microarray dataset with timepoints restricted to timepoints of the RNAseq dataset
.PHONY : restrictedTimepointsRunning
restrictedTimepointsRunning :
	$(R_SCRIPT) $(SCRIPTS_DIR)algo_and_normalization_through_limited_timepoints_mouse_microarray.R $(species) $(tissue)

##°circadGenesAnalysisANDrestrictedTimepointsAnalysis: Run algorithms on known circadian genes microarray dataset with timepoints restricted to timepoints of the RNAseq dataset
.PHONY : circadGenesANDrestrictedTimepointsRunning
circadGenesANDrestrictedTimepointsRunning :
	$(R_SCRIPT) $(SCRIPTS_DIR)algo_and_normalization_through_known_circad_genes_AND_limited_timepoints_mouse_microarray.R $(species) $(tissue)

##°RNAseqDatasetWithLowExpressedGenesRemovedRunning: Run algorithms and normalization on mouse RNAseq dataset with low expressed genes removed
.PHONY : RNAseqDatasetWithLowExpressedGenesRemovedRunning
RNAseqDatasetWithLowExpressedGenesRemovedRunning :
	$(R_SCRIPT) $(SCRIPTS_DIR)algo_and_normalization_through_RNAseq_dataset_with_low_expressed_genes_removed.R $(species) $(tissue)

# pvaluesDistributionQuartiles : Plot the distribution of initial pvalues of 1st and 4th quartiles of gene expression levels
#.PHONY : pvaluesDistributionRNAseqMicroarray
#pvaluesDistributionRNAseqMicroarray :
#	$(R_SCRIPT) $(SCRIPTS_DIR)pvalues_distrib_RNAseq_vs_microarray_mouse.R $(tissue)

# pvaluesDistributionQuartiles : Plot the distribution of initial pvalues of 1st and 4th quartiles of gene expression levels
#.PHONY : pvaluesDistributionRNAseqMicroarrayCircadGenes
#pvaluesDistributionRNAseqMicroarrayCircadGenes :
#	$(R_SCRIPT) $(SCRIPTS_DIR)pvalues_distrib_RNAseq_vs_microarray_mouse_CircadGenes.R $(tissue)

#.PHONY : ProportionRhythmicOrthoWithinQuartilesPlots
#ProportionRhythmicOrthoWithinQuartilesPlots : $(ORTHO_DIR)raw.pvalue_$(species1)_$(tissue1)_$(species2)_$(tissue2).txt
#	$(R_SCRIPT) $(SCRIPTS_DIR)rhythmic_orthologs_proportions_within_quartiles.R $(species1) $(tissue1) $(species2) $(tissue2) $(species2.method) $(species2.threshold) $(species2.pvalue) $(ONLY_WITHIN_CONSERVED_GENES) $(ONLY_ORTHO_ONE_TO_ONE)
