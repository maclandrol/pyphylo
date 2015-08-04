
for file in $(ls $1/*fasta); do
	out="$file"-distmat
	sed -i 's/>.*#/>/g' $file
	fastprot -I fasta -o $out -O phylip $file
	ntax=$(grep '>' $file | wc -l)
	echo "#NEXUS\n\nBEGIN taxa;
   DIMENSIONS ntax=$ntax;
TAXLABELS
$(grep '>' $file | tr -d '>'|awk '{printf("[%d]   \x27%s\x27 \n", NR,$0)}')
;
END;

BEGIN distances;
   DIMENSIONS ntax=$ntax;
   FORMAT
   triangle=BOTH
   diagonal
   labels
   ;

   MATRIX\n$(cat $out | tail -n +2|  sed -e "s/[^ ][^ ]*/'&'/g" |sed -e 's/^/   /' )\n   ;\n   END;" > $out
done

myfastafile="$1/myfastafile"
mydistfile="$1/mydistfile"
ls $1/*fasta > $myfastafile &&
ls $1/*fasta-distmat > $mydistfile &&
distR -d $mydistfile -f $myfastafile -o "$1"
rm -f $mydistfile
rm -f $myfastafile
rm -f "$1/*fasta-distmat"