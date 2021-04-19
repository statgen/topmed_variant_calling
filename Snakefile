import os
import math
#import pdb

#pdb.set_trace()

configfile: "config.yml"
singularity: "topmed-variant-calling.sif"

mem_step_size = 6400

rule discover_bcf:
    input:
        config["input_expression"]
    output:
        "out/sm/{sample_id}/{sample_id}.bcf"
        #"out/sm/{sample_id}/{sample_id}.bcf.csi"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt:  mem_step_size + mem_step_size * attempt
    shell:
        """
        set +e

        tmp_dir=/tmp/snk_`head -c 8 <(pwd | md5sum)`_discover_{wildcards.sample_id}
        mkdir -p $tmp_dir

        set -o pipefail 
        export REF_PATH=resources/ref/md5/%2s/%2s/%s 
        {config[exe_root]}/samtools/samtools view -uh -T resources/ref/hs38DH.fa {input} 2> $tmp_dir/{wildcards.sample_id}.bcf.samtools_err \
          | {config[exe_root]}/bamUtil/bin/bam clipoverlap --poolSize 100000000 --in -.ubam --out -.ubam 2> $tmp_dir/{wildcards.sample_id}.bcf.bamUtil_err \
          | {config[exe_root]}/vt-topmed/vt discover2 -z -q 20 -b + -r resources/ref/hs38DH.fa -s {wildcards.sample_id} -o $tmp_dir/{wildcards.sample_id}.bcf 2> $tmp_dir/{wildcards.sample_id}.bcf.vt_err
        rc=$?
        
        if [[ $rc == 0 ]]; then 
          {config[exe_root]}/bcftools/bcftools index -f $tmp_dir/{wildcards.sample_id}.bcf
          rc=$?
        fi
        
        mv $tmp_dir/*_err out/sm/{wildcards.sample_id}/        
        if [[ $rc == 0 ]]; then 
          mv $tmp_dir/{wildcards.sample_id}.bcf $tmp_dir/{wildcards.sample_id}.bcf.csi out/sm/{wildcards.sample_id}/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """
        
rule sample_qc:
    input:
        config["input_expression"]
    output:
        "out/sm/{sample_id}/{sample_id}.vb2"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e

        tmp_dir=/tmp/snk_`head -c 8 <(pwd | md5sum)`_vb2_{wildcards.sample_id}
        mkdir -p $tmp_dir

        export REF_PATH=resources/ref/md5/%2s/%2s/%s
        
        {config[exe_root]}/cramore/bin/cramore cram-verify-bam \
          --svd resources/ref/HGDP_938.b38.genotypes.svd \
          --sam {input} \
          --cap-DP 100 \
          --out $tmp_dir/{wildcards.sample_id}.vb2 \
          --num-PC 4 \
          > $tmp_dir/{wildcards.sample_id}.vb2.stdout \
          2> $tmp_dir/{wildcards.sample_id}.vb2.stderr
        rc=$?

        mv $tmp_dir/*.stderr $tmp_dir/*.stdout out/sm/{wildcards.sample_id}/
        if [[ $rc == 0 ]]; then
          mv $tmp_dir/{wildcards.sample_id}.vb2 out/sm/{wildcards.sample_id}/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule discover_bcf_xy_depths:
    input:
        "out/sm/{sample_id}/{sample_id}.bcf" #rules.discover_bcf.output #"out/sm/{sample_id}/{sample_id}.bcf"
    output:
        "out/sm/{sample_id}/{sample_id}.norm.xy"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e
        
        tmp_dir=/tmp/snk_`head -c 8 <(pwd | md5sum)`_norm_{wildcards.sample_id}
        mkdir -p $tmp_dir

        export REF_PATH=resources/ref/md5/%2s/%2s/%s

        {config[exe_root]}/cramore/bin/cramore vcf-normalize-depth --xy \
          --vcf {input} \
          --known resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz \
          --gc resources/ref/hs38DH.gc.w150.s5.gz \
          --xLabel chrX \
          --yLabel chrY \
          --xStart 2781479 \
          --xStop 15570138 \
          --out $tmp_dir/{wildcards.sample_id}.norm \
          > $tmp_dir/{wildcards.sample_id}.norm.out \
          2> $tmp_dir/{wildcards.sample_id}.norm.err
        rc=$?

        mv $tmp_dir/*.out $tmp_dir/*.err out/sm/{wildcards.sample_id}/
        if [[ $rc == 0 ]]; then
          mv $tmp_dir/{wildcards.sample_id}.norm.xy out/sm/{wildcards.sample_id}/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

def get_batch_sample_ids(wc):
    batch_size = config["batch_size"]
    beg = batch_size * int(wc.batch)
    end = beg + batch_size
    return config["sample_ids"][beg:end]

def get_batch_cram_files(wc):
    batch_size = config["batch_size"]
    beg = batch_size * int(wc.batch)
    end = beg + batch_size
    return [config["input_expression"].format(sample_id=sid) for sid in config["sample_ids"][beg:end]]


def get_batch_norm_files(wc):
    return ["out/sm/" + sid + "/" + sid + ".norm.xy" for sid in get_batch_sample_ids(wc)]

def get_batch_vb_files(wc):
    return ["out/sm/" + sid + "/" + sid + ".vb2" for sid in get_batch_sample_ids(wc)]

rule group_vb_xy_table:
    input:
        get_batch_norm_files, get_batch_vb_files
    params:
        cram_files = get_batch_cram_files,
        sample_ids = get_batch_sample_ids
    output:
        cram_list = "out/index/cram_list_b{batch}.tsv",
        table = "out/index/vb_xy_{batch}.tsv"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        #set +e
        tmp_dir=`mktemp -d`

        echo {params.sample_ids} | xargs -n1 > $tmp_dir/ids.txt
        echo $tmp_dir/ids.txt
        echo {params.cram_files} | xargs -n1 > $tmp_dir/paths.txt
        echo $tmp_dir/paths.txt
        paste $tmp_dir/ids.txt $tmp_dir/paths.txt > {output.cram_list}
        {config[exe_root]}/apigenome/bin/cram-vb-xy-index --index {output.cram_list} --dir out/sm --out {output.table}
        rc=$?
        rm -r $tmp_dir
        exit $rc
        """

rule batch_site_list:
    input:
       rules.group_vb_xy_table.output.table #"out/index/vb_xy_{batch}.tsv"
    output:
        "out/union/b{batch}/b{batch}_{chrom}_{beg}_{end}.merged.sites.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e

        tmp_dir=/tmp/snk_`head -c 8 <(pwd | md5sum)`_$(basename {output} .bcf)
        mkdir -p $tmp_dir

        {config[exe_root]}/cramore/bin/cramore vcf-merge-candidate-variants \
          --in-vcf-list <(cut -f 3 {input} | tail -n+2) \
          --region "{wildcards.chrom}:{wildcards.beg}-{wildcards.end}" \
          --out-vcf $tmp_dir/out.bcf \
          > $tmp_dir/out.bcf.out \
          2> $tmp_dir/out.bcf.err \
        && {config[exe_root]}/bcftools/bcftools index $tmp_dir/out.bcf
        
        rc=$?
        
        mv $tmp_dir/out.bcf.out {output}.out
        mv $tmp_dir/out.bcf.err {output}.err
        if [[ $rc == 0 ]]; then
          mv $tmp_dir/out.bcf {output} && mv $tmp_dir/out.bcf.csi {output}.csi
          rc=$?
        fi


        rm -r $tmp_dir
        exit $rc
        """

def get_batch_site_files(wc):
    num_batches = math.ceil(len(config["sample_ids"]) / config["batch_size"])
    ret = []
    for b in range(0, num_batches):
        ret.append(("out/union/b{}/b{}_" + wc.chrom + "_" + wc.beg + "_" + wc.end + ".merged.sites.bcf").format(b,b))
    return ret
    

rule union_site_list:
    input:
        get_batch_site_files
    output:
        "out/union/union.{chrom}_{beg}_{end}.merged.sites.bcf",
        "out/union/union.{chrom}_{beg}_{end}.sites.bcf"
    threads: 2
    resources:
        mem_mb = lambda wc, attempt: mem_step_size + mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/cramore/bin/cramore vcf-merge-candidate-variants \
          --in-vcf-list <(echo {input} | tr ' ' '\n') \
          --region {wildcards.chrom}:{wildcards.beg}-{wildcards.end} \
          --out-vcf {output[0]} \
          > {output[0]}.out 2> {output[0]}.err

        {config[exe_root]}/bcftools/bcftools index -f {output[0]}
        {config[exe_root]}/vt-topmed/vt annotate_indels -r resources/ref/hs38DH.fa {output[0]} -o + 2> {output[0]}.annotate.err \
          | {config[exe_root]}/vt-topmed/vt consolidate_variants + -o {output[1]} > {output[1]}.out 2> {output[1]}.err
        {config[exe_root]}/bcftools/bcftools index -f {output[1]}
        """

def get_regions_for_contig(c_length, region_size):
    ret = []
    num_regions = int(math.ceil(c_length / region_size))
    for i in range(0, num_regions):
        start = i * region_size
        end = start + region_size
        ret.append(str(start + 1) + "_" + str(min(end, c_length)))
    return ret

def get_region_strings(contigs):
    ret = []
    region_size = config["region_size"]
    for c, c_length in contigs.items():
        for r in get_regions_for_contig(c_length, region_size):
            ret.append(c + "_" + r)
    return ret

def get_merge_region_strings(contigs):
    ret = []
    region_size = config["merge_region_size"]
    for c, c_length in contigs.items():
        for r in get_regions_for_contig(c_length, region_size):
            ret.append(c + "_" + r)
    return ret

rule all_union_site_lists:
    input:
        [ancient("out/union/union." + r + ".sites.bcf") for r in get_region_strings(config["contigs"])]

rule genotyped_batch:
    input:
        cram_list = rules.group_vb_xy_table.output.cram_list, #"out/index/cram_list_b{batch}.tsv",
        sex_map = rules.group_vb_xy_table.output.table, #"out/index/vb_xy_{batch}.tsv",
        sites_bcf = rules.union_site_list.output[1] #"out/union/union.{chrom}_{region}.sites.bcf"
    output:
        "out/genotypes/b{batch}/b{batch}.{chrom}_{beg}_{end}.genotypes.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        # out : log/batch-geno
        # list : BATCH : index/seq.batches.by.20.txt
        # list : INTERVAL : index/intervals/b38.intervals.X.10Mb.10Mb.txt
        # var : ROOT : ..
        # var : PREFIX : out/genotypes/batches/$BATCH$1$/b$BATCH$1$.$INTERVAL$1$_$INTERVAL$2$_$INTERVAL$3$
        # name : example-batch-genotype
        # target : $PREFIX$.genotypes.bcf $PREFIX$.genotypes.bcf.csi
        #mkdir -p out/genotypes/batches/$BATCH$1$/
        #cut -f 1,20 out/index/list.107.local.crams.vb_xy.index | tail -n +2 | tail -n +$BATCH$1$ | head -n 20 > $PREFIX$.sex_map.txt
        #cut -f 1,2,5 out/index/list.107.local.crams.vb_xy.index | tail -n +2 | tail -n +$BATCH$1$ | head -n 20 > $PREFIX$.cram_index.txt
        tmp_dir=`mktemp -d`
        set +e
        
        VB_DEPTH_IDX=4; FREEMIX_IDX=5; FRAC_DP10_IDX=16
        tail -n+2 {input.sex_map} | awk "\$$VB_DEPTH_IDX > {config[thresholds][vb_depth]} && \$$FREEMIX_IDX < {config[thresholds][freemix]} && \$$FRAC_DP10_IDX > {config[thresholds][frac_dp10]} {{print ;}}" | cut -f 1,20 > ${{tmp_dir}}/sex_map.txt
        tail -n+2 {input.sex_map} | awk "\$$VB_DEPTH_IDX > {config[thresholds][vb_depth]} && \$$FREEMIX_IDX < {config[thresholds][freemix]} && \$$FRAC_DP10_IDX > {config[thresholds][frac_dp10]} {{print ;}}" | cut -f 1,2,5 > ${{tmp_dir}}/cram_index.txt

        #cut -f 1,20 {input.sex_map} | grep -v "^SAMPLE_ID" > ${{tmp_dir}}/sex_map.txt
        #cut -f 1,2,5 {input.sex_map} | grep -v "^SAMPLE_ID" > ${{tmp_dir}}/cram_index.txt
       
        tmp_out=${{tmp_dir}}/$(basename {output})
 
        export REF_PATH=resources/ref/md5/%2s/%2s/%s 
        {config[exe_root]}/cramore/bin/cramore dense-genotype \
          --in-cram-list ${{tmp_dir}}/cram_index.txt \
          --in-vcf {input.sites_bcf} \
          --unit 1000000 \
          --region {wildcards.chrom}:{wildcards.beg}-{wildcards.end} \
          --sex-map ${{tmp_dir}}/sex_map.txt \
          --xLabel chrX \
          --yLabel chrY \
          --xStart 2781479 \
          --xStop 155701383 \
          --print-tmp-info \
          --out ${{tmp_out}} \
          --min-mq 1 \
          > ${{tmp_out}}.out 2> ${{tmp_out}}.err &&
        {config[exe_root]}/bcftools/bcftools index -f ${{tmp_out}}
        rc=$?
        mv ${{tmp_out}}.out ${{tmp_out}}.err $(dirname {output})/
        if [[ $rc == 0 ]]; then
          mv ${{tmp_out}} ${{tmp_out}}.csi $(dirname {output})/
          rc=$?
        fi
        rm -r $tmp_dir
        exit $rc
        """

rule entire_batch:
    input:
       ["out/genotypes/b{batch}/b{batch}." + r + ".genotypes.bcf" for r in get_region_strings(config["contigs"])]
    output:
        "out/b{batch}.done"
    threads: 1
    shell:
        """
        touch {output}
        """

rule all_batches:
    input:
        ["out/b" + str(b) + ".done" for b in range(0, math.ceil(len(config["sample_ids"]) / config["batch_size"]))]


rule full_sex_map:
    input:
        ["out/index/vb_xy_{batch}.tsv".format(batch=b) for b in range(0, math.ceil(len(config["sample_ids"]) / config["batch_size"]))]
    output:
        "out/index/sex_map.txt",
        "out/index/merged_vb_xy.tsv"
    threads: 1
    resources:
        mem_mb = mem_step_size
    shell:
        """
        for f in {input}; do
          if [[ $f == "out/index/vb_xy_0.tsv" ]]; then
            cat $f
          else
            tail -n+2 $f
          fi
        done > {output[1]}

        tail -n+2 {output[1]} | cut -f1,20 > {output[0]}
        """


def get_paste_input(wc):
    ret = []
    region_sz=int(config["region_size"])
    beg = int((int(wc.beg) - 1) / region_sz) * region_sz + 1
    end = min(beg + region_sz - 1, config["contigs"][wc.chrom])

    for b in range(0, math.ceil(len(config["sample_ids"]) / config["batch_size"])):
        ret.append("out/genotypes/b{batch}/b{batch}.{chrom}_{beg}_{end}.genotypes.bcf".format(batch=b, chrom=wc.chrom, beg=beg, end=end))
    return ret


rule pasted_genotype_region:
    input:
        sex_map = rules.full_sex_map.output,
        bcfs = get_paste_input
    output:
        "out/genotypes/merged/{chrom}/merged.{chrom}_{beg}_{end}.genotypes.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        # out : log/paste-geno
        # list : INTERVAL : index/intervals/b38.intervals.X.10Mb.1Mb.txt
        # var : ROOT : ..
        # var : PREFIX : $INTERVAL$1$/merged.$INTERVAL$1$_$INTERVAL$4$_$INTERVAL$5$
        # name : example-paste-genotype
        # target : out/genotypes/merged/$PREFIX$.genotypes.bcf out/genotypes/merged/$PREFIX$.genotypes.bcf.csi out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf out/genotypes/minDP0/$PREFIX$.gtonly.minDP0.bcf.csi out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf out/genotypes/minDP10/$PREFIX$.gtonly.minDP10.bcf.csi out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf out/genotypes/hgdp/$PREFIX$.gtonly.minDP0.hgdp.bcf.csi
        set +e
        tmp_dir=`mktemp -d`

        #mkdir -p out/genotypes/merged/{wildcards.chrom}/
        #mkdir -p out/genotypes/minDP0/{wildcards.chrom}/
        #mkdir -p out/genotypes/minDP10/{wildcards.chrom}/
        #mkdir -p out/genotypes/hgdp/{wildcards.chrom}/
        
        #cut -f1,20 `ls -v out/index/vb_xy_*.tsv` | grep -v "^SAMPLE_ID" > ${{tmp_dir}}/sex_map.txt
        echo {input.bcfs} | tr ' ' '\n' > ${{tmp_dir}}/merged.bcflist.txt
 
        merged_file=${{tmp_dir}}/$(basename {output})
        
        tmp_info_flag="--skip-tmp-info"
        #tmp_info_flag=""

        {config[exe_root]}/cramore/bin/cramore vcf-paste-calls \
          --vcf-list ${{tmp_dir}}/merged.bcflist.txt \
          --num-pc 0 \
          --sex-map {input.sex_map} \
          --xLabel chrX --yLabel chrY \
          --mtLabel chrM \
          --xStart 2781479 --xStop 155701383 \
          $tmp_info_flag \
          --region {wildcards.chrom}:{wildcards.beg}-{wildcards.end} \
          --out ${{merged_file}} \
          > ${{merged_file}}.out 2> ${{merged_file}}.err &&
        {config[exe_root]}/bcftools/bcftools index -f ${{merged_file}}
        rc=$?
        mv ${{merged_file}}.out ${{merged_file}}.err $(dirname {output})/
        
        if [[ $rc == 0 ]]; then
          mv ${{merged_file}} ${{merged_file}}.csi $(dirname {output})/
          rc=$?
        fi
        rm -r ${{tmp_dir}}
        exit $rc
        """


rule all_merged_genotypes:
    input:
        ["out/genotypes/merged/" + r.split("_")[0] + "/merged." + r + ".genotypes.bcf" for r in get_merge_region_strings(config["contigs"])]



def get_merge_regions_for_region(chrom, mbeg, mend):
    ret = []
    region_sz=int(config["region_size"])
    merge_region_sz=int(config["merge_region_size"])

    for i in range(0, region_sz // merge_region_sz):
        beg = int(mbeg) + i * merge_region_sz
        end = min(beg + merge_region_sz - 1, config["contigs"][chrom])
        ret.append(str(beg) + "_" + str(end))
        if end == config["contigs"][chrom]:
            break
    return ret


rule mindp_region:
    input:
        sex_map = rules.full_sex_map.output,
        bcfs = rules.pasted_genotype_region.output lambda wc: [rules.pasted_genotype_region.output.format(chrom=wc.chrom, region=r) for r in get_merge_regions_for_region(wc.chrom, wc.beg, wc.end)]
    output:
        "out/genotypes/minDP{min_dp}/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP{min_dp}.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell: 
        """
        set +e
        tmp_dir=`mktemp -d`
        mkdir ${{tmp_dir}}/chunks        
        mindp_file=${{tmp_dir}}/$(basename {output})
        
        #cut -f1,20 `ls -v out/index/vb_xy_*.tsv` | grep -v "^SAMPLE_ID" > ${{tmp_dir}}/sex_map.txt

        for input_chunk in {input.bcfs}; do
          chunk_mindp_file=${{tmp_dir}}/chunks/$(basename ${{input_chunk}} .bcf).mindp_chunk.bcf
          {config[exe_root]}/cramore/bin/cramore vcf-squeeze \
            --in {input} \
            --sex-map {input.sex_map} \
            --x-label chrX --y-label chrY --mt-label chrM \
            --x-start 2781479 --x-stop 155701383 \
            --minDP-male-X $(( {wildcards.min_dp} / 2 )) \
            --minDP {wildcards.min_dp} \
            --out ${{chunk_mindp_file}} \
            > ${{mindp_file}}.out 2> ${{mindp_file}}.err
          rc=$?
          [[ $rc != 0 ]] && break;
        done

        if [[ $rc == 0 ]]; then
          {config[exe_root]}/bcftools/bcftools concat --naive $(ls -v ${{tmp_dir}}/chunks/*.mindp_chunk.bcf) -O -o ${{mindp_file}} &&
          {config[exe_root]}/bcftools/bcftools index -f ${{mindp_file}}
        fi
 
        mv ${{mindp_file}}.out ${{mindp_file}}.err $(dirname {output})/

        if [[ $rc == 0 ]]; then
          mv ${{mindp_file}} ${{mindp_file}}.csi $(dirname {output})/
          rc=$?
        fi
        rm -r ${{tmp_dir}}
        exit $rc
        """

rule hgdp_region:
    input:
        "out/genotypes/minDP0/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP0.bcf"
    output:
        "out/genotypes/hgdp/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP0.hgdp.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e
        tmp_dir=`mktemp -d`
        hgdp_file=${{tmp_dir}}/$(basename {output})

        {config[exe_root]}/cramore/bin/cramore vcf-extract \
          --vcf {input} \
          --site resources/ref/HGDP_938.hg38.sites.vcf.gz \
          --out ${{hgdp_file}} \
          > ${{hgdp_file}}.out 2> ${{hgdp_file}}.err &&
        {config[exe_root]}/bcftools/bcftools index -f ${{hgdp_file}}
        rc=$?
        mv ${{hgdp_file}}.out ${{hgdp_file}}.err $(dirname {output})/

        if [[ $rc == 0 ]]; then
          mv ${{hgdp_file}} ${{hgdp_file}}.csi $(dirname {output})/
          rc=$?
        fi
        rm -r ${{tmp_dir}}
        exit $rc
        """


def get_autosome_region_strings(contigs, autosome_contigs, region_size):
    ret = []
    for c in autosome_contigs:
        for r in get_regions_for_contig(contigs[c], region_size):
            ret.append(c + "_" + r)
    return ret


rule all_gtonly_genotypes:
    input:
        ["out/genotypes/minDP0/" + r.split("_")[0] + "/merged." + r + ".gtonly.minDP0.bcf" for r in get_merge_region_strings(config["contigs"])] + ["out/genotypes/minDP10/" + r.split("_")[0] + "/merged." + r + ".gtonly.minDP10.bcf" for r in get_merge_region_strings(config["contigs"])]

rule all_hgdp_genotypes:
    input:
        ["out/genotypes/hgdp/" + r.split("_")[0] + "/merged." + r + ".gtonly.minDP0.hgdp.bcf" for r in get_autosome_region_strings(config["contigs"], config["autosome_contigs"], config["region_size"])]


rule concatenated_hgdp_autosomes:
    input:
        ["out/genotypes/hgdp/" + r.split("_")[0] + "/merged." + r + ".gtonly.minDP0.hgdp.bcf" for r in get_autosome_region_strings(config["contigs"], config["autosome_contigs"], config["region_size"])]
    output:
        "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/bcftools/bcftools concat -n -Ob -o {output} {input}
        """

rule autosomes_bed_file:
    input: rules.concatenated_hgdp_autosomes.output
    output: "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.plink.bed"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell: 
        """
        output_file={output}
        plink-1.9 --bcf {input} --make-bed --out ${{output_file%.bed}} --allow-extra-chr --double-id
        """

rule kinship:
    input: rules.autosomes_bed_file.output
    output: "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.kin0"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        output_file={output}
        {config[exe_root]}/king/king -b {input} --cpus {threads} --related --degree 1 --kinship --prefix ${{output_file%.kin0}}
        """

rule pedigree:
    input:
        sex_map = rules.full_sex_map.output,
        kin = rules.kinship.output
    output: 
        "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.inferred.ped"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/apigenome/bin/vcf-infer-ped --kin0 {input.kin} --sex {input.sex_map} --out {output}
        """

rule filtered_pedigree:
    input:
        ped = rules.pedigree.output,
        bcf = rules.concatenated_hgdp_autosomes.output
    params:
        mendel_prefix = "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.mendel"
    output: "out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.king.inferred.filtered.ped" #"out/genotypes/hgdp/merged.autosomes.gtonly.minDP0.hgdp.mendel.ind.dup.conc"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        cut -f 1,16 out/index/merged_vb_xy.tsv | tail -n+2 | perl -lane 'print $F[0] if ( $F[1] < 0.95)' > out/genotypes/hgdp/exclude.fdp10_lt_95_percent.ids

        {config[exe_root]}/cramore/bin/cramore mendel-dup-conc --vcf {input.bcf} --ped {input.ped} --out {params.mendel_prefix}

        {config[exe_root]}/apigenome/bin/vcf-filter-ped-by-mendel --ped {input.ped} --mendel {params.mendel_prefix} --exclude out/genotypes/hgdp/exclude.fdp10_lt_95_percent.ids --out {output}
        
        #scripts_dir=/net/wonderland/home/lefaivej/topmed_variant_calling_master
        #perl $scripts_dir/scripts/c02-filter-ped.pl {input.ped} {params.mendel_prefix} out/index/merged_vb_xy.tsv > {output}
        """

rule milk_filtered_region:
    input: 
        pedigree = rules.filtered_pedigree.output,
        bcf = lambda wc: ["out/genotypes/merged/{chrom}/merged.{chrom}_{region}.genotypes.bcf".format(chrom=wc.chrom, region=r) for r in get_merge_regions_for_region(wc.chrom, wc.beg, wc.end)]
    output:
        "out/milk/{chrom}/milk.{chrom}_{beg}_{end}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e
        tmp_dir=`mktemp -d`

        for input_file in {input.bcf}; do
          tmp_chunk_out=${{tmp_dir}}/chunks/$(basename $input_file .bcf).milk.vcf.gz
          {config[exe_root]}/vt-topmed/vt milk_filter \
            -f {input.pedigree} \
            -b ${{input_file}} \
            -o ${{tmp_chunk_out}} \
            -g out/index/sex_map.txt \
            --xLabel chrX --yLabel chrY --mtLabel chrM \
            --xStart 2781479 --xStop 155701383 \
            --af-field AF
          rc=$?
        done
       
        tmp_out=${{tmp_dir}}/$(basename {output} .sites.vcf.gz).full.vcf.gz

        if [[ $rc == 0 ]]; then
          {config[exe_root]}/bcftools/bcftools concat --naive $(ls -v ${{tmp_dir}}/chunks/*.milk.vcf.gz) -Oz -o ${{tmp_out}} &&
          rm -r ${{tmp_dir}}/chunks/ &&
          zcat ${{tmp_out}} | cut -f 1-8 | {config[exe_root]}/htslib/bgzip -c > ${{tmp_dir}}/$(basename {output})
          rc=$?
        fi
 
        if [[ $rc == 0 ]]; then
          # TODO: only move sites file
          mv ${{tmp_dir}}/* $(dirname {output})/
          rc=$?
        fi
        
        rm -r ${{tmp_dir}}
        exit $rc
        """


def get_concatenated_milk_chromome_input(wc):
    return ["out/milk/" + wc.chrom + "/milk." + wc.chrom + "_" + r + ".sites.vcf.gz" for r in get_regions_for_contig(config["contigs"][wc.chrom], config["region_size"])]

rule concatenated_milk_chromosome:
    input:
       get_concatenated_milk_chromome_input
    output:
       "out/milk/milk.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        {config[exe_root]}/bcftools/bcftools concat {input} -Oz -o $tmp_out &&
        {config[exe_root]}/htslib/tabix -f -pvcf $tmp_out
        rc=$?
        
        if [[ $rc == 0 ]]; then
          mv $tmp_dir/* $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule all_concatenated_milk_chromosomes:
    input:
        ["out/milk/milk." + c + ".sites.vcf.gz" for c, l in config["contigs"].items()]


rule missingness_statistics:
    input:
        rules.mindp_region.output
    output:
       "out/missingness_update/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP{min_dp}.missigness_update.sites.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt 
    shell:
        """
        set +e
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        {config[exe_root]}/cramore/bin/cramore vcf-update-info \
          --in {input} \
          --out $tmp_out \
          --update-AC --update-GC --update-FMIS --site-only 
        rc=$?
        
        if [[ $rc == 0 ]]; then
          mv $tmp_out $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule concatenated_missignenss_statistics:
    input:
       lambda wc: ["out/missingness_update/" + wc.chrom + "/merged." + wc.chrom + "_" + r + ".gtonly.minDP10.missigness_update.sites.bcf" for r in get_regions_for_contig(config["contigs"][wc.chrom], config["region_size"])]
    output:
        "out/missingness_update/merged.{chrom}.gtonly.minDP10.missigness_update.sites.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/bcftools/bcftools concat --naive {input} > {output}
        """


rule milk_with_fmis_sites:
    input:
        rules.concatenated_milk_chromosome.output, rules.concatenated_missignenss_statistics.output
    output:
        "out/svm/input/milk_with_fmis.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        set +e
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})

        export EXE_PREFIX={config[exe_root]}
        perl {config[exe_root]}/scripts/d13-add-fmis-to-frz9.pl {input} > $tmp_out &&
        {config[exe_root]}/htslib/tabix -f -pvcf $tmp_out;
        rc=$?
         
        if [[ $rc == 0 ]]; then
          mv $tmp_out $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """

rule svm_model:
    input:
        "out/svm/input/milk_with_fmis.chr2.sites.vcf.gz" #rules.milk_with_fmis_sites.output
    params:
       output_prefix = "out/svm/milk_svm.chr2"
    output:
        "out/svm/milk_svm.chr2.svm.model"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/apigenome/bin/vcf-svm-milk-filter \
          --ignore AVG_IF --ignore LLK0 --ignore BETA_IF \
          --in-vcf {input} \
          --out {params.output_prefix} \
          --ref resources/ref/hs38DH.fa \
          --dbsnp resources/ref/dbsnp_142.b38.vcf.gz \
          --posvcf resources/ref/hapmap_3.3.b38.sites.vcf.gz \
          --posvcf resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz \
          --train --cutoff -0.5 --min-disc-count 5 --min-disc-frac 0.1 \
          --centromere resources/ref/hg38.centromere.bed.gz \
          --bgzip {config[exe_root]}/htslib/bgzip \
          --tabix {config[exe_root]}/htslib/tabix \
          --invNorm {config[exe_root]}/invNorm/bin/invNorm \
          --svm-train {config[exe_root]}/libsvm/svm-train \
          --svm-predict {config[exe_root]}/libsvm/svm-predict
        """


rule svm_filtered_vcf:
    input:
        vcf = rules.milk_with_fmis_sites.output,
        model = rules.svm_model.output
    params:
        output_prefix = "out/svm/milk_svm.{chrom}"
    output:
        "out/svm/milk_svm.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        {config[exe_root]}/apigenome/bin/vcf-svm-milk-filter \
          --ignore AVG_IF --ignore LLK0 --ignore BETA_IF \
          --in-vcf {input.vcf} \
          --out {params.output_prefix} \
          --ref resources/ref/hs38DH.fa \
          --dbsnp resources/ref/dbsnp_142.b38.vcf.gz \
          --posvcf resources/ref/hapmap_3.3.b38.sites.vcf.gz \
          --posvcf resources/ref/1000G_omni2.5.b38.sites.PASS.vcf.gz \
          --cutoff -0.5 --min-disc-count 5 --min-disc-frac 0.1 \
          --model {input.model} \
          --centromere resources/ref/hg38.centromere.bed.gz \
          --bgzip {config[exe_root]}/htslib/bgzip \
          --tabix {config[exe_root]}/htslib/tabix \
          --invNorm {config[exe_root]}/invNorm/bin/invNorm \
          --svm-train {config[exe_root]}/libsvm/svm-train \
          --svm-predict {config[exe_root]}/libsvm/svm-predict
#          --svm-train {config[exe_root]}/libsvm/svm-train \
#          --svm-predict {config[exe_root]}/libsvm/svm-predict
        """

rule hard_filtered_sites:
    input:
        "out/svm/milk_svm.{chrom}.sites.vcf.gz", "out/milk/milk.{chrom}.sites.vcf.gz" #rules.svm_filtere_vcf.output
    params:
       output_prefix = "out/hard_filter/milk_svm_hard.{chrom}.sites"
    output:
        "out/hard_filter/milk_svm_hard.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        export EXE_PREFIX={config[exe_root]}
        perl {config[exe_root]}/scripts/e04-filter-vars-fixedX.pl {input} {wildcards.chrom} {params.output_prefix}
        """

rule all_hard_called_sites:
    input:
        [rules.hard_filtered_sites.output[0].format(chrom="chr" + str(c)) for c in range(1,23)]


rule gwas_catalog_override:
    input:
        rules.hard_filtered_sites.output
    params:
        output_prefix = "out/catalog/milk_svm_hard_catalog.{chrom}.sites"
    output:
        "out/catalog/milk_svm_hard_catalog.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        export EXE_PREFIX={config[exe_root]}
        perl {config[exe_root]}/scripts/e05-whitelist-gwas-variants.pl {wildcards.chrom} {input} {params.output_prefix}
        """

rule anno_database:
    output: directory("out/anno/GRCh38.86")
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        # -dataDir needs to be absolute path. i
        database_dir=$(dirname {output})
        if [[ "${{database_dir:0:1}}" != / ]]; then
          database_dir=$(pwd)/$database_dir
        fi

        java -jar {config[exe_root]}/snpEff/snpEff.jar download -v -dataDir $database_dir GRCh38.86
        """

#https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download
rule annotated_sites:
    input:
        vcf = rules.gwas_catalog_override.output,
        database = rules.anno_database.output
    #params:
    #    output_prefix = "out/catalog/milk_svm_hard_catalog.{chrom}.sites"
    output:
        "out/anno/milk_svm_hard_catalog_anno.{chrom}.sites.vcf.gz"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        ls resources/dbSNP/b151/b38/All_20180418_withchr.vcf.gz

        database_dir=$(dirname {input.database})
        if [[ "${{database_dir:0:1}}" != / ]]; then
          database_dir=$(pwd)/$database_dir
        fi
 
        export EXE_PREFIX={config[exe_root]}
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
        set -euo pipefail
       
        zcat {input.vcf} | java -Xmx4g -jar {config[exe_root]}/snpEff/snpEff.jar -v -stats $(dirname {output})/$(basename {output} .vcf.gz).summary.html GRCh38.86 -dataDir $database_dir -nodownload | {config[exe_root]}/htslib/bgzip -c > $tmp_out

        {config[exe_root]}/apigenome/bin/vcf-add-rsid --vcf $tmp_out --db resources/dbSNP/b151/b38/All_20180418_withchr.vcf.gz | {config[exe_root]}/htslib/bgzip -c > {output}
        
        rm -r $tmp_dir
        """


rule filtered_genotypes:
    input:
        sites = rules.gwas_catalog_override.output,
        genotypes = "out/genotypes/minDP{min_dp}/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP{min_dp}.bcf" #rules.mindp_region.output
    output: "out/filtered_genotypes/minDP{min_dp}/{chrom}/merged.{chrom}_{beg}_{end}.gtonly.minDP{min_dp}.filtered.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        #sudo -u hmkang bash -c 'tabix -h /mnt/gcsfuse/topmed-vt/freeze9/release/sites/freeze9.merged.$INTERVAL$1$.filtered.anno.gwas.sites.vcf.gz $INTERVAL$1$:$INTERVAL$4$-$INTERVAL$5$ | bgzip -c > /home/hmkang/sites.vcf.gz'
        
        #{config[exe_root]}/cramore/bin/cramore vcf-update-sites --md-info GWAS \
        #  --rm-info AF --rm-info GC --rm-info GN --rm-info HWEAF_P --rm-info AVG_IF \
        #  --md-info SVM --md-info ANN --replace-filter --replace-id \
        #  --md-info FMIS10 --md-info MAX_HM3_R2 --md-info LOG_HM3_BUDDY_03 --md-info LOG_HM3_BUDDY_05 \
        #  --md-info LOG_HM3_BUDDY_08 --md-info LOG_HM3_DIST_03 --md-info LOG_HM3_DIST_05 \
        #  --md-info LOG_HM3_DIST_08 --md-info DUP_NH_ALL --md-info DUP_NH_DIS --md-info TRI_NH_ALL \
        #  --md-info TRI_NH_DIS \
        #  --vcf {input.genotypes} \
        #  --site {input.sites} \
        #  --region {wildcards.chrom}:{wildcards.region} \
        #  --out {output}

        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
        set -uo pipefail
        
        {config[exe_root]}/cramore/bin/cramore vcf-update-sites \
          --rm-info AF \
          --rm-info GC \
          --rm-info GN \
          --rm-info HWEAF_P \
          --rm-info AVG_IF \
          --rm-info FLT20 \
          --rm-info NS_NREF \
          --replace-filter --replace-id \
          --md-info SVM \
          --md-info ANN \
          --md-info FMIS10 \
          --md-info DUP_NH_ALL \
          --md-info DUP_NH_DIS \
          --md-info TRI_NH_ALL \
          --md-info TRI_NH_DIS \
          --md-info GWAS \
          --md-info CLINVAR \
          --md-info CLINVARB \
          --vcf {input.genotypes} \
          --site {input.sites} \
          --region {wildcards.chrom}:{wildcards.beg}-{wildcards.end} \
          --out $tmp_out
        rc=$?
        
        if [[ $rc == 0 ]]; then
          mv ${{tmp_out}} $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        
        #{config[exe_root]}/bcftools/bcftools index -m 1 -f {output}
        """




rule concatenated_filtered_genotypes:
    input:
        lambda wc: [rules.filtered_genotypes.output[0].format(chrom=wc.chrom, min_dp=wc.min_dp, region=r) for r in get_regions_for_contig(config["contigs"][wc.chrom], config["region_size"])]
    output: "out/filtered_genotypes/minDP{min_dp}/merged.{chrom}.gtonly.minDP{min_dp}.filtered.bcf"
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: mem_step_size * attempt
    shell:
        """
        tmp_dir=`mktemp -d`
        tmp_out=$tmp_dir/$(basename {output})
        set -uo pipefail

        {config[exe_root]}/bcftools/bcftools concat --naive {input} -o $tmp_out &&
        {config[exe_root]}/bcftools/bcftools index -m 1 $tmp_out
        rc=$?

        if [[ $rc == 0 ]]; then
          mv ${{tmp_out}} ${{tmp_out}}.csi $(dirname {output})/
          rc=$?
        fi

        rm -r $tmp_dir
        exit $rc
        """
