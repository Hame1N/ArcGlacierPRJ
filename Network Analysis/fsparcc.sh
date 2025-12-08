###Author:zhangzh
###DATE:2023.02.16

helpDoc(){
cat <<EOF
Usage: bash $0 <otutab>  [output file prefix]

    otutab file:
        relative abundance
    output file:
        default: <otutab >".txt"
EOF
}

if [ $# -lt 2 ]
then
    helpDoc
    exit 1
fi

FF=$1
OUT=$2
echo "first step" 
fastspar  -c ${FF} -i 100 -x 20 -t 10 -e 0.1 -r ${OUT}_cor.txt -a ${OUT}_cov.txt
echo "second step" 
mkdir bootstrap_counts
mkdir bootstrap_correlation
fastspar_bootstrap --otu_table ${FF} --number 1000 --prefix bootstrap_counts/${OUT}_data
echo "third step" 
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*
echo "final step" 
fastspar_pvalues --otu_table ${FF}  --correlation ${OUT}_cor.txt  --prefix bootstrap_correlation/cor_${OUT}_data_ --permutations 1000 --outfile pvalues.tsv
echo "All done"