EXE="../src/pairsnp"

echo "check version"
$EXE -v

echo "check good"
$EXE ./good.aln > temp.txt
cmp temp.txt test_good.txt

echo "check sparse"
$EXE -s ./good.aln > temp.txt
cmp temp.txt test_sparse.txt

echo "check filter"
$EXE -s -d 2 ./good.aln > temp.txt
cmp temp.txt test_filter.txt

echo "check singleton"
$EXE ./singleton.aln > temp.txt
cmp temp.txt test_singleton.txt

echo "check lowercase"
$EXE ./lowercase.aln > temp.txt
cmp temp.txt test_lowercase.txt

echo "check ambig"
$EXE ./ambig.aln > temp.txt
cmp temp.txt test_ambig.txt

$EXE -n ./ambig.aln > temp.txt
cmp temp.txt test_ambig_with_n.txt

$EXE ./ambig_cons.aln > temp.txt
cmp temp.txt test_ambig_cons.txt

echo "check bad"
$EXE ./bad.aln > temp.txt || true
cmp temp.txt test_bad.txt

rm temp.txt
