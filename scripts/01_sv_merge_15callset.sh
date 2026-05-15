for r in grch38 chm13 hap1 hap2; do 
    for s in colo829 h1437 h2009 hcc1937 hcc1954 hg008; do 
        j=`for i in $(ls ../raw_calls/*/${r}_${s}_*vcf.gz ); do echo $i; done`; 
        echosv merge -i $j -o ${r}_${s}_merge_15callset.vcf.gz --merge --new ; 
        echosv merge -i ${r}_${s}_merge_15callset.vcf.gz -o ${r}_${s}_merge_15callset_highconfi.vcf.gz --extract ;
    done ; 
done