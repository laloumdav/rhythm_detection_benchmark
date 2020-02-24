# VERTEBRATES
make OrthologyPerGene species1=baboon species2=mouse_RNAseq tissue=liver &
make OrthologyPerGene species1=baboon species2=zebrafish tissue=liver &
make OrthologyPerGene species1=baboon species2=mouse_RNAseq tissue=lung &
wait

make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=adrenal_gland &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=aorta &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=brain_stem &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=heart &
wait
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=cerebellum &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=kidney &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=hypothalamus &
wait
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=lung &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=liver &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=muscle &
make OrthologyPerGene species1=mouse_RNAseq species2=baboon tissue=white_adipose &
wait

make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=adrenal_gland &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=aorta &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=brain_stem &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=heart &
wait
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=cerebellum &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=kidney &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=hypothalamus &
wait
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=lung &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=liver &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=muscle &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=white_adipose &
wait
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue=scn &

make OrthologyPerGene species1=baboon species2=zebrafish tissue=liver &
make OrthologyPerGene species1=baboon species2=rat tissue=lung &

make OrthologyPerGene species1=mouse_microarray species2=zebrafish tissue=liver &
wait
make OrthologyPerGene species1=mouse_microarray species2=rat tissue=lung &

make OrthologyPerGene species1=mouse_RNAseq species2=zebrafish tissue=liver &
make OrthologyPerGene species1=mouse_RNAseq species2=rat tissue2=lung &


make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=liver tissue2=aorta &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=liver tissue2=lung &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=liver tissue2=kidney &
wait
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=liver tissue2=adrenal_gland &

make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=lung tissue2=aorta &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=lung tissue2=liver &
wait
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=lung tissue2=kidney &
make OrthologyPerGene species1=mouse_microarray species2=baboon tissue1=lung tissue2=adrenal_gland &

wait
# INSECTS
make OrthologyPerGene species1=drosophila species2=anopheles tissue=head &
make OrthologyPerGene species1=drosophila species2=anopheles tissue=body &
wait
make OrthologyPerGene species1=drosophila species2=aedes tissue=head &
make OrthologyPerGene species1=aedes species2=anopheles tissue=head