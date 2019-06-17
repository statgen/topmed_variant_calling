# out : log/discover
# list : SMIDX : index/list.107.local.crams.index
# var : ROOT : ..
# target : out/sm/$SMIDX$1$/$SMIDX$1$.vb2 out/sm/$SMIDX$1$/$SMIDX$1$.norm.xy out/sm/$SMIDX$1$/$SMIDX$1$.bcf out/sm/$SMIDX$1$/$SMIDX$1$.bcf.csi
# name: example-discovery
mkdir -p out/sm/$SMIDX$1$/
bash -c 'set -o pipefail; REF_PATH=resources/ref/md5/%2s/%2s/%s $ROOT$/samtools/samtools view -uh -T resources/ref/hs38DH.fa $SMIDX$2$ 2> out/sm/$SMIDX$1$/$SMIDX$1$.bcf.samtools_err | $ROOT$/bamUtil/bin/bam clipoverlap --poolSize 100000000 --in -.ubam --out -.ubam 2> out/sm/$SMIDX$1$/$SMIDX$1$.bcf.bamUtil_err | $ROOT$/vt-topmed/vt discover2 -z -q 20 -b + -r resources/ref/hs38DH.fa -s $SMIDX$1$ -o out/sm/$SMIDX$1$/$SMIDX$1$.bcf 2> out/sm/$SMIDX$1$/$SMIDX$1$.bcf.vt_err'
$ROOT$/bcftools/bcftools index -f out/sm/$SMIDX$1$/$SMIDX$1$.bcf
REF_PATH=resources/ref/md5/%2s/%2s/%s $ROOT$/cramore/cramore cram-verify-bam --svd resources/ref/HGDP_938.b38.genotypes.svd --sam $SMIDX$2$ --cap-DP 100 --out out/sm/$SMIDX$1$/$SMIDX$1$.vb2 --num-PC 4 > out/sm/$SMIDX$1$/$SMIDX$1$.vb2.stdout 2> out/sm/$SMIDX$1$/$SMIDX$1$.vb2.stderr
$ROOT$/cramore/cramore vcf-normalize-depth --xy --vcf out/sm/$SMIDX$1$/$SMIDX$1$.bcf --known resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz --gc resources/ref/hs38DH.gc.w150.s5.gz --xLabel chrX --yLabel chrY --xStart 2781479 --xStop 15570138 --out out/sm/$SMIDX$1$/$SMIDX$1$.norm > out/sm/$SMIDX$1$/$SMIDX$1$.norm.out 2> out/sm/$SMIDX$1$/$SMIDX$1$.norm.err
