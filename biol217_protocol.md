# Day 1 - Introduction to linux

- become familiar with the Linux terminal
- try some commands

current directory
```
pwd
```
list files in current directory
```
ls -l
```
create new directory
```
mkdir day_01
```
change directory
```
cd day_01
cd ../ #upwards
```
move file 
```
mv test.txt dir/day_01 
```
remove directory
```
rmdir day_01
```


# **Day 2 - From raw reades to contigs**

## **Data**
https://ami-journals.onlinelibrary.wiley.com/doi/full/10.1111/1751-7915.13313

## **Miniconda**
 ```
module load miniconda3/4.7.12.1
conda activate anvio
cd /work_beegfs/sunam232/day_02
for i in *.gz; do fastqc $i; done#
 ```
### **Script (fastqc)**

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH -D ./
#SBATCH --output=anviscript.out
#SBATCH --error=anviscript.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


#load your anvio environment (path needs to be adjusted)
source activate /home/sunamXXX/miniconda3/miniconda4.9.2/usr/etc/profile.d/conda.sh/envs/anvio-7.1

#set environmental variables here (anything you want to run multiple times)
outdir=/work_beegfs/sunamXXX/out_dir
meta=/work_beegfs/sunamXXX/.txt
runname=Run2023-01-20
trunclength=240
`fehlt noch was`
```
### **Script (fastp)**

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=all
#SBATCH --reservation=biol217



#fastp
cd /work_beegfs/sunam232/day_02
mkdir ../clean_reads

for i in `ls *_R1.fastq.gz`;
do
    second=`echo ${i} | sed 's/_R1/_R2/g'`
    fastp -i ${i} -I ${second} -R _report -o ../clean_reads/"${i}" -O ../clean_reads/"${second}" -t 6 -q 20

done
```
### **Script (megahit)**

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=all
#SBATCH --reservation=biol217



#megahit

cd /work_beegfs/sunam232/day_02/clean

                                       
megahit -1 BGR_130305_R1.fastq.gz -1 BGR_130527_R1.fastq.gz -1 BGR_130708_R1.fastq.gz -2 BGR_130305_R2.fastq.gz -2 BGR_130527_R2.fastq.gz -2 BGR_130708_R2.fastq.gz --min-contig-len 1000 --presets meta-large -m 0.85 -o /work_beegfs/sunam232/day_02/assembly -t 20                      

```

# **Day 3 - From raw reads to contigs**

### **folder structure**

2_fastq 
3_coassembly 
3_metaquast 
4_mapping 
5_anvio_profiles


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

in first Terminal (CAU Cluster)

```
(anvio) [sunam232@caucluster2 day_03]$ cp -r /home/sunam226/Day3/contigs.(anvio) [sunam232@caucluster2 day_03]$ anvio.fa /work_beegfs/sunam232/day_03/4_mapping/
(anvio) [sunam232@caucluster2 day_03]$ sbatch anviscript_binning 
(anvio) [sunam232@caucluster2 day_03]$ sbatch anviscript_bowtie2 
```
abgebrochen, stattdessen heruntergeladen

```
(anvio) [sunam232@caucluster2 day_03]$ cp -r /home/sunam226/Day3/4_mapping /work_beegfs/sunam232/day_03/4_mapping
```

# **Day 4 - From bins to species and abundance estimation**











