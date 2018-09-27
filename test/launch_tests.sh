EXE="../src/pairsnp"

echo "check version"
$EXE -v

echo "check singleton"
$EXE ./singleton.aln
$EXE -n ./singleton.aln

echo "check lowercase"
$EXE ./lowercase.aln
$EXE -n ./lowercase.aln

echo "check ambig"
$EXE ./ambig.aln
$EXE -n ./ambig.aln

echo "check good"
$EXE ./good.aln
$EXE -n ./good.aln
$EXE -c ./good.aln
$EXE -s ./good.aln

echo "check bad"
$EXE ./bad.aln
$EXE -n ./bad.aln

