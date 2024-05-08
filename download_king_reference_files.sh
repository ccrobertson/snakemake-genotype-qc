#1KG reference files for KING Ancestry Inference
#https://www.kingrelatedness.com/ancestry/

mkdir -p resources
cd resources
wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz
xz -d -v KGref.bed.xz
xz -d -v KGref.bim.xz
xz -d -v KGref.fam.xz
