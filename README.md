# snakme-runhic

adapt from xiaotao's repo [runHiC](http://xiaotaowang.github.io/HiC_pipeline/).

for gcb cluster submit the jobs using sbatch 

modify the resource for each rule as needed.


snakemake -j 10  --cluster-config cluster.json \
 --cluster "sbatch -J {cluster.job} --mem={cluster.mem} -N 1 \
 -n {threads}  -o {cluster.out} -e  {cluster.err} "  &> log & 


