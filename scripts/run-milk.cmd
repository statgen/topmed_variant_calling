# out : /net/topmed8/working/call_sets/joint.131k/out/milk
# list : INTERVAL : /net/topmed8/working/call_sets/joint.131k/index/intervals/b38.intervals.X.10Mb.100kb.txt
# target : gs://topmed-vt-us/joint.131k/milk/$INTERVAL$1$/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.full.vcf.gz gs://topmed-vt-us/joint.131k/milk/$INTERVAL$1$/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.sites.vcf.gz gs://topmed-vt-us/joint.131k/milk/$INTERVAL$1$/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.sites.vcf.gz.tbi
# name: milk-131k
# image: ubuntu-20180417-vt
# disk: 200
# zone: us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d
# gs: gs://topmed-vt-us/joint.131k/log/milk
# opts: --custom-cpu 1 --custom-memory 6656MiB
# preemptible: true
# keep: false
# checkpoint-secs: 300
# max-hours: 22
# sleep-secs: 10
sudo chown -R hmkang:hmkang /home/hmkang
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang gcsfuse --implicit-dirs topmed-vt-us /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang bash -c 'cut -f 1,20 /mnt/gcsfuse/topmed-vt-us/joint.131k/index/joint.131k.2018_02b.pass.index > /home/hmkang/sex_map.txt'
sudo -u hmkang vt milk_filter -f /mnt/gcsfuse/topmed-vt-us/joint.131k/index/hgdp.autosomes.gtonly.minDP0.king.inferred.ped -b /mnt/gcsfuse/topmed-vt-us/joint.131k/genotypes/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf -o /home/hmkang/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.full.vcf.gz -g /home/hmkang/sex_map.txt --xLabel chrX --yLabel chrY --mtLabel chrM --xStart 2781479 --xStop 155701383 --af-field AF
sudo -u hmkang bash -c 'zcat /home/hmkang/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.full.vcf.gz | cut -f 1-8 | bgzip -c > /home/hmkang/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.sites.vcf.gz'
sudo -u hmkang /home/hmkang/tools/htslib.orig/tabix -pvcf /home/hmkang/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.sites.vcf.gz
sudo -u hmkang gsutil cp /home/hmkang/milk.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.* gs://topmed-vt-us/joint.131k/milk/$INTERVAL$1$/
