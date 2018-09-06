# out : /net/topmed8/working/call_sets/joint.131k/out/paste
# list : INTERVAL : /net/topmed8/working/call_sets/joint.131k/index/intervals/b38.intervals.X.10Mb.100kb.txt
# target : gs://topmed-vt-us/joint.131k/genotypes/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf gs://topmed-vt-us/joint.131k/genotypes/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf.csi gs://topmed-vt-us/joint.131k/gtonly/minDP0/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf gs://topmed-vt-us/joint.131k/gtonly/minDP0/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.csi gs://topmed-vt-us/joint.131k/hgdp/minDP0/$INTERVAL$1$/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf gs://topmed-vt-us/joint.131k/hgdp/minDP0/$INTERVAL$1$/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.csi gs://topmed-vt-us/joint.131k/gtonly/minDP10/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf gs://topmed-vt-us/joint.131k/gtonly/minDP10/$INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf.csi
# name: paste-131k
# image: ubuntu-20180417-vt
# disk: 200
# zone: us-central1-a us-central1-b us-central1-c us-central1-f
# gs: gs://topmed-vt-us/joint.131k/log/paste
# opts: --custom-cpu 1 --custom-memory 6656MiB
# preemptible: false
# keep: false
# checkpoint-secs: 300
# max-hours: 22
# sleep-secs: 10
sudo chown -R hmkang:hmkang /home/hmkang
sudo -u hmkang mkdir --p /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang gcsfuse --implicit-dirs topmed-vt-us /mnt/gcsfuse/topmed-vt-us
sudo -u hmkang bash -c 'cut -f 1,10,11 /mnt/gcsfuse/topmed-vt-us/joint.131k/index/joint.131k.2018_02b.pass.index > /home/hmkang/sample_evecs.txt'
sudo -u hmkang bash -c 'cut -f 1,20 /mnt/gcsfuse/topmed-vt-us/joint.131k/index/joint.131k.2018_02b.pass.index > /home/hmkang/sex_map.txt'
sudo -u hmkang cat /mnt/gcsfuse/topmed-vt-us/joint.131k/index/batches/seq.batches.joint.txt | xargs -I {} echo '/mnt/gcsfuse/topmed-vt-us/joint.131k/genotypes/batches/{}/{}.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$.genotypes.bcf' > /home/hmkang/bcf_list.txt
sudo -u hmkang bash -c 'set -o pipefail; cramore vcf-paste-calls --vcf-list /home/hmkang/bcf_list.txt --evec /home/hmkang/sample_evecs.txt --num-pc 2 --sex-map /home/hmkang/sex_map.txt --xLabel chrX --yLabel chrY --mtLabel chrM --xStart 2781479 --xStop 155701383 --skip-tmp-info --out /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf --region $INTERVAL$1$:$INTERVAL$4$-$INTERVAL$5$ > /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf.out'
sudo -u hmkang bcftools index -f /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf
sudo -u hmkang cramore vcf-squeeze --in /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf --minDP 0 --out /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf >  /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.out
sudo -u hmkang bcftools index -f /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf
sudo -u hmkang cramore vcf-squeeze --in /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf --minDP 10 --out /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf >  /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf.out
sudo -u hmkang bcftools index -f /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf
sudo -u hmkang cramore vcf-extract --vcf /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf --site /home/alignment/ref/HGDP_938.hg38.sites.vcf.gz --out /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf > /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.out 2> /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.err
sudo -u hmkang bcftools index -f /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf
sudo -u hmkang gsutil cp /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.genotypes.bcf.csi gs://topmed-vt-us/joint.131k/genotypes/$INTERVAL$1$/
sudo -u hmkang gsutil cp /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.csi gs://topmed-vt-us/joint.131k/gtonly/minDP0/$INTERVAL$1$/
sudo -u hmkang gsutil cp /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf /home/hmkang/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP10.bcf.csi gs://topmed-vt-us/joint.131k/gtonly/minDP10/$INTERVAL$1$/
sudo -u hmkang gsutil cp /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf /home/hmkang/hgdp.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$.gtonly.minDP0.bcf.csi gs://topmed-vt-us/joint.131k/hgdp/minDP0/$INTERVAL$1$/
