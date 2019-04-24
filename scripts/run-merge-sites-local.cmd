# out : log/merge
# list : BATCH : index/seq.batches.by.20.txt
# list : INTERVAL : index/intervals/b38.intervals.X.10Mb.10Mb.txt
# var : ROOT : ..
# var : PREFIX : out/union/$BATCH$1$/b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$
# target : $PREFIX$.merged.sites.bcf $PREFIX$.merged.sites.bcf.csi
# name: example-merge
mkdir --p out/union/$BATCH$1$/
cut -f 3 out/index/list.107.local.crams.vb_xy.index | tail -n +2 | tail -n +$BATCH$1$ | head -n 20 > $PREFIX$.bcflist.txt
$ROOT$/cramore/cramore vcf-merge-candidate-variants --in-vcf-list $PREFIX$.bcflist.txt --region $INTERVAL$1$:$INTERVAL$2$-$INTERVAL$3$ --out-vcf $PREFIX$.merged.sites.bcf > $PREFIX$.merged.sites.bcf.out 2> $PREFIX$.merged.sites.bcf.err
$ROOT$/bcftools/bcftools index $PREFIX$.merged.sites.bcf
