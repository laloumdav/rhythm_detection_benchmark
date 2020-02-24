make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=aorta species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=brain_stem species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=heart species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=cerebellum species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=kidney species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=hypothalamus species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=lung species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=muscle species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=scn species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait

make pvaluesDistribOrthologs species1=baboon species2=zebrafish tissue=liver species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=baboon species2=rat tissue=lung species2.method=meta2d species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait

make pvaluesDistribOrthologs species1=mouse_microarray species2=zebrafish tissue=liver species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=rat tissue=lung species2.method=meta2d species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait

make pvaluesDistribOrthologs species1=mouse_microarray species2=zebrafish tissue=liver species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=raw.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=rat tissue=lung species2.method=meta2d species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=raw.pvalue &
wait

make pvaluesDistribOrthologs species1=mouse_RNAseq species2=zebrafish tissue=liver species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_RNAseq species2=rat tissue2=lung species2.method=RAIN species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait

make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=JTK species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=JTK species2.threshold=0.05 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=ARS species2.threshold=0.05 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=mouse_microarray species2=baboon tissue=liver species2.method=RAIN species2.threshold=0.05 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait


make pvaluesDistribOrthologs species1=drosophila species2=anopheles tissue=head species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=drosophila species2=anopheles tissue=body species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=drosophila species2=aedes tissue=head species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=anopheles species2=drosophila tissue=body species2.method=JTK species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=aedes species2=drosophila tissue=head species2.method=JTK species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue &
wait
make pvaluesDistribOrthologs species1=aedes species2=anopheles tissue=head species2.method=ARS species2.threshold=0.01 ONLY_WITHIN_CONSERVED_GENES=TRUE ONLY_ORTHO_ONE_TO_ONE=TRUE P_VAL=default.pvalue





