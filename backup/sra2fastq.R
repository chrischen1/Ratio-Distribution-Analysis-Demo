system('module load sratoolkit/sratoolkit-2.8.1-2')


for(i in list.files('/storage/htc/bdm/Collaboration/minRen_rnaseq/sra',full.names = T)){
  system(paste('sbatch ~/configs/bash_sub.sbatch fastq-dump -I --split-files',i))
}

rawdata_path = '/storage/htc/bdm/Collaboration/minRen_rnaseq/fastq/'
trim_data_path = rawdata_path
script_dir = '/storage/htc/bdm/ccm3x/rnaseq_pipeline/'
ref_genome = '/storage/htc/bdm/ref_genomes/tair10/star_ref'
alignment_result = '/storage/htc/bdm/Collaboration/minRen_rnaseq/results_tair10/'

# for(i in list.files(rawdata_path,pattern = '.+.fastq$')){
#   infile <- paste(rawdata_path,i,sep = '')
#   outfile <- paste(trim_data_path,i,sep = '') 
#   system(paste('sbatch ',script_dir,'P2_trim_data.sbatch ',infile,' ',outfile,sep = ''))
# }

system('module load star/star-2.5.2b')
fastq_files <- list.files(trim_data_path,pattern = '.+.fastq.*')
for(i in unique(gsub('_\\d.fastq','',fastq_files))){
  out_path_i <- paste(alignment_result,gsub('.fastq.*','',i),'/',sep = '')
  dir.create(out_path_i,showWarnings = F)
  infile <- paste(trim_data_path,i,c('_1','_2'),'.fastq',sep = '')
  system(paste('sbatch ',script_dir,'P3_alignment_with_star.sbatch ',ref_genome,' "',infile[1],' ',infile[2],'" ',out_path_i,sep = ''))
}

ercc_percentage <- c()
for(i in list.dirs('/storage/htc/bdm/Collaboration/minRen_rnaseq/results_0.66/',recursive=F)){
  ercc_percentage <- c(ercc_percentage,read.delim(paste0(i,'/Log.final.out'),sep='\t',as.is=T)[8,2])
}
names(ercc_percentage) <- gsub('.+//','',list.dirs('/storage/htc/bdm/Collaboration/minRen_rnaseq/results_0.66/',recursive=F))
meta <- read.csv('/storage/htc/bdm/Collaboration/minRen_rnaseq/meta.csv',as.is = T,header = F,row.names = 1)
names(ercc_percentage) <- grep('RNA',value = T,meta$V2)