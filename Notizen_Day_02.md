# Day 2
## From raw reades to contigs
### **Vorlesung**


Metagenome: hole Genes -> wrong to say amplicon Metagenomic

### **Link**

https://github.com/AammarTufail/Bioinformatics_Master_Module2023/blob/main/Day-2/Tutorial_Day2.md

### **Data**

https://ami-journals.onlinelibrary.wiley.com/doi/full/10.1111/1751-7915.13313

### **Miniconda**

module load miniconda3/4.7.12.1

conda activate anvio


### **fastqc**
cd /work_beegfs/sunam232/day_02
for i in *.gz; do fastqc $i; done#

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

#commands are printed into the log
set -xu

#navigate to working directory
cd $HOME

# write your command here, if you want to use some environemtal variables, these can be accessed via $nameofvariable

# HERE YOU WRITE YOUR COMMANDS
#
#
rm -rf /home/sunam226/.conda/envs/viromics

# Finish each script with this (prints Done and a Unicorn into your logile), then you know everything has run through
echo "Done removing"

printf '\U1F984\n'

#this prints the required resources into your logfile
jobinfo

```

### **fastp**

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

### **megahit**

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

