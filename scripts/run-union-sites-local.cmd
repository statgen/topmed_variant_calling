# out : log/merge
# list : INTERVAL : index/intervals/b38.intervals.X.10Mb.10Mb.txt
# var : ROOT : /net/fantasia/home/hmkang/code/working/topmed_variant_calling/
# var : PREFIX : out/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$
# target : $PREFIX$.sites.bcf $PREFIX$.sites.bcf.csi
# name: example-union
bash -c 'cat index/seq.batches.by.20.txt | xargs -I {} echo out/union/{}/b{}.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf > $PREFIX$.bcflist.txt'
$ROOT$/cramore/cramore vcf-merge-candidate-variants --in-vcf-list $PREFIX$.bcflist.txt --region $INTERVAL$1$:$INTERVAL$2$-$INTERVAL$3$ --out-vcf $PREFIX$.merged.sites.bcf > $PREFIX$.merged.sites.bcf.out 2> $PREFIX$.merged.sites.bcf.err
$ROOT$/bcftools/bcftools index -f $PREFIX$.merged.sites.bcf
bash -c 'set -o pipefail; $ROOT$/vt-topmed/vt annotate_indels -r resources/ref/hs38DH.fa $PREFIX$.merged.sites.bcf -o + 2> $PREFIX$.annotated.sites.bcf.err | $ROOT$/vt-topmed/vt consolidate_variants + -o $PREFIX$.sites.bcf > $PREFIX$.sites.bcf.out 2> $PREFIX$.bcf.sites.err'
$ROOT$/bcftools/bcftools index -f $PREFIX$.sites.bcf
