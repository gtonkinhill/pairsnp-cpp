EXE="../pairsnp"

echo "check version"
$EXE -v

echo "check good"
$EXE ./good.aln > temp.txt
cmp temp.txt good.out

echo "check sparse"
$EXE -s ./good.aln > temp.txt
cmp temp.txt sparse.out

echo "check filter"
$EXE -s -d 2 ./good.aln > temp.txt
cmp temp.txt filter.out

echo "check singleton"
$EXE ./singleton.aln > temp.txt
cmp temp.txt singleton.out

echo "check lowercase"
$EXE ./lowercase.aln > temp.txt
cmp temp.txt lowercase.out

echo "check ambig"
$EXE ./ambig.aln > temp.txt
cmp temp.txt ambig.out

echo "check ambig2"
$EXE ./ambig_cons.aln > temp.txt
cmp temp.txt ambig_cons.out

echo "check bad"
$EXE ./bad.aln > temp.txt || true
cmp temp.txt bad.out

rm temp.txt
