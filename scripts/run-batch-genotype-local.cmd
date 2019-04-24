# out : log/batch-geno
# list : BATCH : index/seq.batches.by.20.txt
# list : INTERVAL : index/intervals/b38.intervals.X.10Mb.10Mb.txt
# var : ROOT : ..
# var : PREFIX : out/genotypes/batches/$BATCH$1$/b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$
# name : example-batch-genotype
# target : $PREFIX$.genotypes.bcf $PREFIX$.genotypes.bcf.csi
mkdir --p out/genotypes/batches/$BATCH$1$/
cut -f 1,20 out/index/list.107.local.crams.vb_xy.index | tail -n +2 | tail -n +$BATCH$1$ | head -n 20 > $PREFIX$.sex_map.txt
cut -f 1,2,5 out/index/list.107.local.crams.vb_xy.index | tail -n +2 | tail -n +$BATCH$1$ | head -n 20 > $PREFIX$.cram_index.txt
bash -c 'set -o pipefail; REF_PATH=resources/ref/md5/%2s/%2s/%s $ROOT$/cramore/cramore dense-genotype --in-cram-list $PREFIX$.cram_index.txt --in-vcf out/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf --unit 6000000 --region $INTERVAL$1$:$INTERVAL$2$-$INTERVAL$3$ --sex-map $PREFIX$.sex_map.txt --xLabel chrX --yLabel chrY --xStart 2781479 --xStop 155701383 --print-tmp-info --out $PREFIX$.genotypes.bcf --min-mq 1 > $PREFIX$.genotypes.bcf.out 2> $PREFIX$.genotypes.bcf.err'
$ROOT$/bcftools/bcftools index -f $PREFIX$.genotypes.bcf
