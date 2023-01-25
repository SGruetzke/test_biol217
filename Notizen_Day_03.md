# Day 2
## From raw reades to contigs
### **Vorlesung**

Binning
- reference-based
- reference-free
- Tetranucleotide (o k-mer composition)


- fragmentation
- completeness
- contamination
- heterogenetic

### **Folders**

2_fastq
3_coassembly
3_metaquast
4_mapping
5_anvio_profiles

### **Script**

kann sein, dass manchmal hin und hernavigieren fehlt

```

ssh -X sunam232@caucluster.rz.uni-kiel.de
cd $WORK
mkdir day_03
cd ./day_03
mkdir 3_assembly
mkdir 3_metaquast
mkdir 4_mapping
mkdir 5_anvio-profiles
cp -r /home/sunam226/Day3/* /work_beegfs/sunam232/day_03
conda activate /home/sunam226/.conda/envs/anvio
cd ./3_coassembly/
grep -c ">" final.contigs.fa
megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg
cd ../
sbatch anviscript_metaquast

```
didn´t work, results copy worked:

```
cp -r /home/sunam226/Day3/3_metaquast_out/* /work_beegfs/sunam232/day_03/3_metaquast_2

```

new Terminal not on CAU Cluster

```
(base) kurs@Kurs006:~$ cd Downloads/^C
(base) kurs@Kurs006:~$ ^C
(base) kurs@Kurs006:~$ cd Desktop/Bandage/
(base) kurs@Kurs006:~/Desktop/Bandage$ ls
Bandage  dependencies  sample_LastGraph
(base) kurs@Kurs006:~/Desktop/Bandage$ ./Bandage
```

dann wieder im CAU Cluster

```
(anvio) [sunam232@caucluster2 day_03]$ cp -r /home/sunam226/Day3/contigs.(anvio) [sunam232@caucluster2 day_03]$ anvio.fa /work_beegfs/sunam232/day_03/4_mapping/
(anvio) [sunam232@caucluster2 day_03]$ sbatch anviscript_binning 
(anvio) [sunam232@caucluster2 day_03]$ sbatch anviscript_bowtie2 
```

abgebrochen, stattdessen heruntergeladen

```
(anvio) [sunam232@caucluster2 day_03]$ cp -r /home/sunam226/Day3/4_mapping /work_beegfs/sunam232/day_03/4_mapping
```


##**File formats**

**fastq**
-> (clean up) -> **fastq** -> (assembly) -> **fasta** -> (mapping) -> **sam** -> (converte) -> **bam**