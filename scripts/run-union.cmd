# out : /net/topmed8/working/call_sets/joint.131k/out/union
# list : INTERVAL : /net/topmed8/working/call_sets/freeze6/index/intervals/b38.intervals.XYM.10Mb.10Mb.txt
# target : gs://topmed-vt-us/joint.131k/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf gs://topmed-vt-us/joint.131k/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.csi gs://topmed-vt-us/joint.131k/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf gs://topmed-vt-us/ccdg.freeze1/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf.csi
# name: union-131k
# image: ubuntu-20180127-vt
# disk: 150
# zone: us-east1-c us-east1-d us-east1-b us-central1-a us-central1-b us-central1-c us-central1-f
# gs: gs://topmed-vt-us/joint.131k/log/union
# opts: --custom-cpu 2 --custom-memory 13312MiB
# preemptible: true
# keep: false
# checkpoint-secs: 300
# max-hours: 4
# sleep-secs: 60
sudo chown -R hmkang:hmkang /home/hmkang
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang gcsfuse --implicit-dirs topmed-vt-us /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang cp /mnt/gcsfuse/topmed-vt-us/topmed.freeze6/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf /home/hmkang/
sudo -u hmkang bcftools index -f /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf
sudo -u hmkang gsutil cp /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.csi gs://topmed-vt-us/topmed.freeze6/union/
sudo -u hmkang cp /mnt/gcsfuse/topmed-vt-us/ccdg.freeze1/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf /home/hmkang/
sudo -u hmkang bcftools index -f /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf
sudo -u hmkang gsutil cp /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.csi gs://topmed-vt-us/ccdg.freeze1/union/
sudo -u hmkang cp /mnt/gcsfuse/topmed-vt-us/inpsyght.wave5/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf /home/hmkang/
sudo -u hmkang bcftools index -f /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf
sudo -u hmkang gsutil cp /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.csi gs://topmed-vt-us/inpsyght.wave5/union/
sudo -u hmkang cramore vcf-merge-candidate-variants --in-vcf /mnt/gcsfuse/topmed-vt-us/topmed.freeze6/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf --in-vcf /mnt/gcsfuse/topmed-vt-us/ccdg.freeze1/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf --in-vcf /mnt/gcsfuse/topmed-vt-us/inpsyght.wave5/union/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf --region $INTERVAL$1$:$INTERVAL$2$-$INTERVAL$3$ --out-vcf /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf > /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.out
sudo -u hmkang bcftools index -f /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf
sudo -u hmkang bash -c 'set -o pipefail; vt annotate_indels -r /home/alignment/ref/hs38DH.fa /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf -o + 2> /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.annotated.sites.bcf.err | vt consolidate_variants + -o /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf 2> /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf.err'
sudo -u hmkang bcftools index -f /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf
sudo -u hmkang gsutil cp /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.merged.sites.bcf.csi /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf /home/hmkang/union.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.sites.bcf.csi gs://topmed-vt-us/joint.131k/union/
