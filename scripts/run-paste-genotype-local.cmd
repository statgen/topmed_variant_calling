# out : log/paste-geno
# list : INTERVAL : index/intervals/b38.intervals.X.10Mb.1Mb.txt
# var : ROOT : ..
# var : PREFIX : $INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$
# name : example-paste-genotype
# target : out/genotypes/merged/$PREFIX$.genotypes.bcf out/genotypes/merged/$PREFIX$.genotypes.bcf.csi out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf.csi out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf.csi out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf.csi
mkdir -p out/genotypes/merged/$INTERVAL$1$/
mkdir -p out/genotypes/minDP0/$INTERVAL$1$/
mkdir -p out/genotypes/minDP10/$INTERVAL$1$/
mkdir -p out/genotypes/hgdp/$INTERVAL$1$/
cut -f 1,20 out/index/list.107.local.crams.vb_xy.index | tail -n +2 > out/genotypes/merged/$PREFIX$.sex_map.txt
cat index/seq.batches.by.20.txt | xargs -I {} echo 'out/genotypes/batches/{}/b{}.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf' > out/genotypes/merged/$PREFIX$.bcflist.txt
$ROOT$/cramore/cramore vcf-paste-calls --vcf-list out/genotypes/merged/$PREFIX$.bcflist.txt --num-pc 0 --sex-map out/genotypes/merged/$PREFIX$.sex_map.txt --xLabel chrX --yLabel chrY --mtLabel chrM --xStart 2781479 --xStop 155701383 --skip-tmp-info --region $INTERVAL$1$:$INTERVAL$4$-$INTERVAL$5$ --out out/genotypes/merged/$PREFIX$.genotypes.bcf > out/genotypes/merged/$PREFIX$.genotypes.bcf.out 2> out/genotypes/merged/$PREFIX$.genotypes.bcf.err
$ROOT$/bcftools/bcftools index -f out/genotypes/merged/$PREFIX$.genotypes.bcf
$ROOT$/cramore/cramore vcf-squeeze --in out/genotypes/merged/$PREFIX$.genotypes.bcf --minDP 0 --out out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf > out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf.out 2> out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf.err
$ROOT$/bcftools/bcftools index -f out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf
$ROOT$/cramore/cramore vcf-squeeze --in out/genotypes/merged/$PREFIX$.genotypes.bcf --minDP 10 --out out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf > out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf.out 2> out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf.err
$ROOT$/bcftools/bcftools index -f out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf
$ROOT$/cramore/cramore vcf-extract --vcf out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf --site resources/ref/HGDP_938.hg38.sites.vcf.gz --out out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf > out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf.out 2>  out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf.err
$ROOT$/bcftools/bcftools index -f out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf
