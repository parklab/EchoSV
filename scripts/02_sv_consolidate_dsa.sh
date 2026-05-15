for s in colo829 h1437 h2009 hcc1937 hcc1954 hg008; do 
    for r in hap1 hap2; do
        # Genotype SVs from the $r call set
        echosv genotype --longread -i ${r}_${s}_merge_15callset_highconfi.vcf.gz -b \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}_hifi_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}_ont_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}bl_hifi_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}bl_ont_haplotagged.bam \
        -o dsa_merging/${r}_${s}_merge_15callset_highconfi_lr.vcf.gz;

        echosv genotype --shortread -i ${r}_${s}_merge_15callset_highconfi.vcf.gz -b \
        /home/yuz006/disk/srs_alignment/cancer_public/${r}_${s}_ill.bam \
        /home/yuz006/disk/srs_alignment/cancer_public/${r}_${s}bl_ill.bam \
        -o dsa_merging/${r}_${s}_merge_15callset_highconfi_sr.vcf.gz;
    done;
    cat > "dsa_merging/dsa_merge_${s}.json" <<EOF
{
  "refs": {
    "1": "hap1",
    "2": "hap2"
  },
  "vcfs": {
    "1": "dsa_merging//hap1_${s}_merge_15callset_highconfi_lr.vcf.gz",
    "2": "dsa_merging//hap2_${s}_merge_15callset_highconfi_lr.vcf.gz"
  },
  "chains": {
    "2_to_1": "/home/yuz006/disk/t2t_assembly/cancer_public/${s}bl_hap2_hap1.chain.gz"
  },
  "output": "dsa_merging/dsa_merge_${s}_comparison.txt"
}
EOF
    echosv match -i dsa_merging/dsa_merge_${s}.json --multiplat --merge --filter ;
    mv dsa_merging/merge_${s}_merge_15callset_highconfi_lr.vcf.gz ./dsa_${s}_merge_15callset_highconfi_lr.vcf.gz ;
    mv dsa_merging/merge_${s}_merge_15callset_highconfi_sr.vcf.gz ./dsa_${s}_merge_15callset_highconfi_sr.vcf.gz ;
done