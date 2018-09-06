# out : /net/topmed8/working/call_sets/joint.131k/out/geno
# list : BATCH : /net/topmed8/working/call_sets/joint.131k/index/batches/seq.batches.east.txt
# list : INTERVAL : /net/topmed8/working/call_sets/joint.131k/index/intervals/b38.intervals.X.10Mb.10Mb.txt
# target : gs://topmed-vt-us/joint.131k/genotypes/batches/east.b$BATCH$1$/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf gs://topmed-vt-us/joint.131k/genotypes/batches/east.b$BATCH$1$/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf.csi
# name: joint-east
# image: ubuntu-20180312-vt
# disk: 150
# zone: us-east1-b us-east1-c us-east1-d
# gs: gs://topmed-vt-us/joint.131k/log/geno
# opts: --custom-cpu 1 --custom-memory 6656MiB
# preemptible: false
# keep: false
# checkpoint-secs: 300
# max-hours: 22
# sleep-secs: 10
sudo chown -R hmkang:hmkang /home/hmkang
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-ccdg-us
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-ccdg-east1
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang gcsfuse --implicit-dirs topmed-ccdg-us /mnt/gcsfuse/topmed-ccdg-us
sudo -u hmkang gcsfuse --implicit-dirs topmed-ccdg-east1 /mnt/gcsfuse/topmed-ccdg-east1
sudo -u hmkang gcsfuse --implicit-dirs topmed-vt-us /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang bash -c 'cut -f 1,3,9 /mnt/gcsfuse/topmed-vt-us/joint.131k/index/batches/joint.131k.2018_02b.pass.east.b$BATCH$1$.index > /home/hmkang/cram_index.txt'
sudo -u hmkang bash -c 'cut -f 1,20 /mnt/gcsfuse/topmed-vt-us/joint.131k/index/batches/joint.131k.2018_02b.pass.east.b$BATCH$1$.index > /home/hmkang/sex_map.txt'
sudo -u hmkang cp /mnt/gcsfuse/topmed-vt-us/joint.131k/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf /mnt/gcsfuse/topmed-vt-us/joint.131k/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf.csi /home/hmkang/
sudo -u hmkang bash -c 'set -o pipefail; REF_PATH=/home/alignment/ref/md5/%2s/%2s/%s cramore dense-genotype --in-cram-list /home/hmkang/cram_index.txt --in-vcf /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf --unit 100000 --region $INTERVAL$1$:$INTERVAL$2$-$INTERVAL$3$ --sex-map /home/hmkang/sex_map.txt --xLabel chrX --yLabel chrY --xStart 2781479 --xStop 155701383 --print-tmp-info --min-mq 1 --out /home/hmkang/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf > /home/hmkang/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf.out'
sudo -u hmkang bcftools index -f /home/hmkang/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf
sudo -u hmkang gsutil cp /home/hmkang/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf /home/hmkang/east.b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf.csi gs://topmed-vt-us/joint.131k/genotypes/batches/east.b$BATCH$1$/
