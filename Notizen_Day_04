
# Day 3

## Vorlesung

### From bins to species and abundance estimation

## Tutorial

```
(base) kurs@Kurs006:~$ ssh -X sunam232@caucluster-old.rz.uni-kiel.de
[sunam232@caucluster2 ~]$ conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1
(anvio-7.1) [sunam232@caucluster2 ~]$ cd $WORK
(anvio-7.1) [sunam232@caucluster2 sunam232]$ cd ./day_04
(anvio-7.1) [sunam232@caucluster2 day_04]$ sbatch anviscript_contigs-db 
```
### anviscript_contigs-db

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=anvio
#SBATCH --output=anvio_out
#SBATCH --error=anvio.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


#anvio


conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

anvi-gen-contigs-database -f /work_beegfs/sunam232/day_03/4_mapping/contigs.anvio.fa -o contigs.db -n 'biol217'
```


```
(anvio-7.1) [sunam232@caucluster2 day_04]$ sbatch anviscript_HMM 
```

### anviscript_HMM
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=hmm
#SBATCH --output=hmm_out
#SBATCH --error=hmm.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


#hmm


anvi-run-hmms -c /work_beegfs/sunam232/day_04/contigs.db
```
```
(anvio-7.1) [sunam232@caucluster2 day_04]$ srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=all /bin/bash
[sunam232@node010 day_04]$ conda activate /home/sunam225/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

(anvio-7.1) [sunam232@node010 day_04]$ anvi-display-contigs-stats contigs.db
```
Open new terminal


```
ssh -L 8060:localhost:8080 sunam323@caucluster-old.rz.uni-kiel.de
ssh -L 8080:localhost:8080 node'010

```

Firefox


```
http://127.0.0.1:8060
```

zweites Terminal schließen
erstes Terminal strg C und dann exit

```
(anvio-7.1) [sunam232@caucluster2 day_04]$ 
```
sbatch anviscript_samtools
```

cd /work_beegfs/sunam232/day_03/4_mapping/4_mapping/

for i in *.bam; do anvi-init-bam $i -o /work_beegfs/sunam232/day_03/5_anvio-profiles/"$i".sorted.bam; done
```

```
cp /home/sunam226/Day3/contigs.db .
(anvio-7.1) [sunam232@caucluster2 5_anvio-profiles]$ 
```
sbatch anviscript_anvios_profile

```
cd /work_beegfs/sunam232/day_03/5_anvio-profiles/
mkdir /work_beegfs/sunam232/day_03/5_anvio-profiles/profiling/
for i in `ls *.sorted.bam | cut -d "." -f 1`; do anvi-profile -i "$i".bam.sorted.bam -c contigs.db -o /work_beegfs/sunam232/day_03/5_anvio-profiles/profiling/”$i”; done
```


sbatch anviscript_anvi-merge 
```

cd /work_beegfs/sunam232/day_03/5_anvio-profiles/

anvi-merge /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/BGR_130305/PROFILE.db /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/BGR_130527/PROFILE.db /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/BGR_130708/PROFILE.db -o /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/merged_profiles_02 -c /work_beegfs/sunam232/day_03/5_anvio-profiles/contigs.db --enforce-hierarchical-clustering

```


sbatch anviscript_metabat2 

```
cd /work_beegfs/sunam232/day_03/5_anvio-profiles/

anvi-cluster-contigs -p /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam232/day_03/5_anvio-profiles/contigs.db -C METABAT --driver metabat2 --just-do-it --log-file log-metabat2

anvi-summarize -p /work_beegfs/sunam232/day_03/5_anvio-profiles/5_anvio_profiles/merged_profiles/PROFILE.db -c /work_beegfs/sunam232/day_03/5_anvio-profiles/contigs.db -o SUMMARY_METABAT -C METABAT
```

```
(anvio-7.1) [sunam232@caucluster2 5_anvio_profiles]$  anvi-estimate-genome-completeness -p merged_profiles/PROFILE.db -c /work_beegfs/sunam232/day_03/5_anvio-profiles/contigs.db --list-collections
c(anvio-7.1) [sunam232@caucluster2 5_anvio-profiles]$ p /home/sunam226/Day3/5_anvio_profiles/SUMMARY*/*.html /work_beegfs/sunam232/day_03/5_anvio-profiles/download_summary/
```


Questions

Number of bins you got from MetaBAT2? 3
Number of bins you got from CONCOCT? 2
Number of bins you got after consolidating the bins? 2