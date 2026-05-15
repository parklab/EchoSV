for s in colo829 h1437 h2009 hcc1937 hcc1954 hg008; do 
    for r in grch38 chm13; do
        # Genotype SVs from the $r call set
        echosv genotype --longread -i ./${r}_${s}_merge_15callset_highconfi.vcf.gz -b \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}_hifi_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}_ont_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}bl_hifi_haplotagged.bam \
        /home/yuz006/disk/alignment/cancer_public/${r}_${s}bl_ont_haplotagged.bam \
        -o ./${r}_${s}_merge_15callset_highconfi_lr.vcf.gz;

        echosv genotype --shortread -i ./${r}_${s}_merge_15callset_highconfi.vcf.gz -b \
        /home/yuz006/disk/srs_alignment/cancer_public/${r}_${s}_ill.bam \
        /home/yuz006/disk/srs_alignment/cancer_public/${r}_${s}bl_ill.bam \
        -o ./${r}_${s}_merge_15callset_highconfi_sr.vcf.gz;
    done;
    cat > "../final_truthset/refs_cmp_${s}.json" <<EOF
{
  "refs": {
    "1": "grch38",
    "2": "chm13",
    "3": "dsa"
  },
  "vcfs": {
    "1": "./grch38_${s}_merge_15callset_highconfi_lr.vcf.gz",
    "2": "./chm13_${s}_merge_15callset_highconfi_lr.vcf.gz",
    "3": "./dsa_${s}_merge_15callset_highconfi_lr.vcf.gz"
  },
  "chains": {
    "2_to_1": "/home/yuz006/disk/t2t_assembly/cancer_public/chm13_to_grch38.chain.gz",
    "3_to_1": "/home/yuz006/disk/t2t_assembly/cancer_public/${s}bl_hap*_grch38.chain.gz"
  },
  "output": "../final_truthset/refs_${s}_corrMerge.txt"
}
EOF
    echosv match -i ../final_truthset/refs_cmp_${s}.json --multiplat ;
done