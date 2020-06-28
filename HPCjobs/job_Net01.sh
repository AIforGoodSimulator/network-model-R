#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=40:mem=480gb
#PBS -N Net01
module load anaconda3/personal

cd ~/network-model-R/

Rscript -e 'rmarkdown::render("HPCjobs/Netv01_cluster.Rmd")'
