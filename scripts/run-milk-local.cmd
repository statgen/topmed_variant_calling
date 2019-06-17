# out : log/milk
# list : INTERVAL : index/intervals/b38.intervals.X.10Mb.1Mb.txt
# var : ROOT : ..
# var : IN_PREFIX : out/genotypes/merged/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$
# var : OUT_PREFIX : out/milk/$INTERVAL$1$/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$
# name : example-milk
# target : $OUT_PREFIX$.full.vcf.gz $OUT_PREFIX$.sites.vcf.gz $OUT_PREFIX$.sites.vcf.gz.tbi
mkdir -p out/milk/$INTERVAL$1$/
$ROOT$/vt-topmed/vt milk_filter -f out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.inferred.ped -b $IN_PREFIX$.genotypes.bcf -o $OUT_PREFIX$.full.vcf.gz -g $IN_PREFIX$.sex_map.txt --xLabel chrX --yLabel chrY --mtLabel chrM --xStart 2781479 --xStop 155701383 --af-field AF
zcat $OUT_PREFIX$.full.vcf.gz | cut -f 1-8 | $ROOT$/htslib/bgzip -c > $OUT_PREFIX$.sites.vcf.gz
$ROOT$/htslib/tabix -pvcf $OUT_PREFIX$.sites.vcf.gz
