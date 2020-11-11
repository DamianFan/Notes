# script for terminal


import os


# get file names
# set pair file names
# set paths - index, sample_1, sample_2, dir_path
# loop and run terminal script

# set dir
index_dir = '/home/disk/caijingtao/0-index/hg19/STAR-ENSEMBL-HG19-PE150'
sample_dir = '/home/disk/caijingtao/5-jichenbo/RNA/trim_out/'
dir_root_path = '/home/disk/caijingtao/5-jichenbo/RNA/star_out/'
gzfile = 'gznames.txt'

# get gz file names and return a list
# the input attr is '.txt'
def get_gz_names(files):
    gz_filenames = []
    files = open(files,'r')
    for line in files.readlines():
        if line.find('fq.gz') != -1 and line.find('.txt') == -1:
            gz_filenames.append(line.replace('\n',''))
    files.close
    gz_filenames.sort()
    return gz_filenames

gz_filenames = get_gz_names(gzfile)
# print('check get_gzfilenames: ', gz_filenames)


for i in range(0,len(gz_filenames),2):
    # print('current i: ',i)
    # print(sample_dir+gz_filenames[i])
    sample_1 = sample_dir+gz_filenames[i]
    # print(sample_dir+gz_filenames[i+1])
    sample_2 = sample_dir+gz_filenames[i+1]
    tar_dir = dir_root_path + gz_filenames[i][:len(gz_filenames)-5]
    terminal_script = '/home/software/bin/STAR --runThreadN 20 --genomeDir ' + index_dir + ' --readFilesCommand gunzip -c --readFilesIn ' + sample_1 + ' ' + sample_2+ ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + tar_dir
    
    os.system(terminal_script)
    
    print('execuated target sample pair: ', sample_1,sample_2)
    # i += 1

print('all steps have been finished')

# run terminal
# def run_terminal(index_dir,sample_dir,target_dir,gz_filenames):
    
#     sample_1 = '/home/disk/caijingtao/5-jichenbo/RNA/trim_out/lncRNA-mRNA-hfBAT-3-m-2_R2_val_2.fq.gz'
#     sample_2 = '/home/disk/caijingtao/5-jichenbo/RNA/trim_out/lncRNA-mRNA-hfBAT-3-m-2_R1_val_1.fq.gz'
#     dir_path = '/home/disk/caijingtao/5-jichenbo/RNA/star_out/lncRNA-mRNA-hfBAT-3-m-2-'
#     terminal_script = '/home/software/bin/STAR --runThreadN 20 --genomeDir ' + index_dir + ' --readFilesCommand gunzip -c --readFilesIn ' + sample_1 + ' ' + sample_2+ ' --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ' + dir_path + ' &'
#     # print('the script is: ',terminal_script)
#     os.system(terminal_script)
#     print('successfully loaded script')
#     return ''












