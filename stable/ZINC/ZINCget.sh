while read p; do wget -nc $p;done<ZINC-downloader-2D-txt.uri 
find . -size 0 -delete
mkdir input
mv *.txt input/
cd input
find . -type f -size +1G | for i in *; do split -l 10000000 --additional-suffix=".txt" -d "$i" "${i%.txt}_"; rm $i;done
