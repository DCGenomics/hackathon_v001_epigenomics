BASEDIR=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/

FILES_ARRAY=$(cat gmxxx_rbbs.txt)

for FILE in ${FILES_ARRAY}
do
    echo $${FILE}
    wget ${BASEDIR}${FILE}
done