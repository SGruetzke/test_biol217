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
copy file
```
cp test02.txt day_01
```
remove directory
```
rmdir day_01
```
'*' is a wildcard. For example, if you want to copy all the text files to a directory, the command could be like this:
```
cp *.txt day_01
```
editing files
```
cat >> test.txt #command
This is a test file. #write the text
#press ctrl + C 
```
merge files
```
cat test.txt test02.txt
```
`slide 39 weiterbearbeiten`

# **Day 2 - From raw reades to contigs**

## **Data**
* samples from mesophilic agricultural biogas plant (near Cologne, Germany)
* over a time period of 587 days, monthly intervals
* by high-throughput amplicon sequencing of the archaeal and bacterial 16S rRNA genes
* https://ami-journals.onlinelibrary.wiley.com/doi/full/10.1111/1751-7915.13313
* we focused on 3 samples in the course

## **Miniconda**
needed packeges and programs were already installed into one conda environment which we activated:
 ```
module load miniconda3/4.7.12.1
conda activate anvio
 ```
 
## Quality control of raw reads
To assess the quality of the sequences, we used FastQC.
This process is very computationally intensive, which is why we did not use the Front End. In such cases, you can use the compute nodes (which have higher computing power) by using a batch script.
Every job script starts with the directive #!/bin/bash on the first line. In the next lines of the script, the job parameters are defined. These always begin with #SBATCH. Job parameters in the batch script for the quality control:
* number of nodes
* numbers of cores per task
* Real memory required per node; G for gigabytes
* walltime
* jobname
* standard output file
* standard error file
* Slurm partition (~batch class)
* use ressources reserved for course 

The job parameters are followed by the command lines for the job.
Our comand line for quality control looked like this: for i in *.gz; do fastqc $i; done
Since the same command was to be executed for several files, a loop was used in which fastqc was applied to all files with the .gz extension. -o always specifies where the newly acquired data should be stored (output folder).

### **Script (fastqc)**

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=anviscript.out
#SBATCH --error=anviscript.err
#SBATCH --partition=all
#SBATCH --reservation=biol217


for i in *.gz; do fastqc $i; done

```
To start the script, enter the following command into the terminal:
```
sbatch fastqc
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











