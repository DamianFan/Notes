# Author: Damain.F
# open files
# set sig status - up/down/non
# save csv files for next steps - gene names, heatmaps, volcano


import math

# open files

files = open('gene_exp.csv','r')


# gene_exp file
# 0 test_id, 1 gene_id, 2 gene[gene name], 3 locus, 4 sample_1, 5 sample_2, 6 status, 7 value_1, 8 value_2, 9 log2(fold_change), 10 test_stat, 11 p_value, 12 q_value, 13 sig[add by Damian]

geneset = [] # fliter sig gene names for next step
diffgeneset = [] # filter only diff genes in geneset
upgeneset = [] # up gene names
downgeneset = [] # down gene names
outputset = [] # output for csv
heatmapset = [] # output for heatmap
volset = [] # output for volcano plot

# start read files
for line in files.readlines():
    row = line.split(',')
    # filter head row
    row[12]=row[12].replace('\n','')
    if row[0] == 'test_id':
        continue
    if row[9] == 'NOTEST' or row[9] == '#NAME?' or  row[9] == '#VALUE!' or row[9] == 'inf':
        continue
    # filter sig gene names and set sig status, relogfc>0 and p<0.05 - up, locfc<0 and p<0.05 - down
    # print('*******************************************',row[9])
    log2reverse = -1 * float(row[9]) # reverse logfc value
    pvalue = float(row[11])
    
    if pvalue <= 0.05 and abs(log2reverse) >1:
        if row[2] not in geneset:
            geneset.append(row[2])
        if row[2] not in upgeneset and log2reverse > 0:
            upgeneset.append(row[2])
            row.append('up')
            heatmapset.append(','.join([row[2],str(math.log(float(row[7]),10)),str(math.log(float(row[8]),10))]))
            diffgeneset.append(','.join(row))
        if row[2] not in downgeneset and log2reverse < 0:
            downgeneset.append(row[2])
            row.append('down')
            heatmapset.append(','.join([row[2],str(math.log(float(row[7]),10)),str(math.log(float(row[8]),10))]))
            diffgeneset.append(','.join(row))
    else:
        row.append('non')
    if float(row[11]) < 1:
            volset.append(','.join([row[9],row[11],row[13]]))
        

    outputset.append(','.join(row))

files.close


# set up csv files of results
# sum gene result

genesigfile = 'gene_exp_all.csv'
with open(genesigfile,'w') as file_object:
    for i in outputset:
        file_object.write(i)
        file_object.write('\n')

diffgenefile = 'diff_gene_updown.csv'
with open(diffgenefile,'w') as diffupdown_object:
    for i in diffgeneset:
        diffupdown_object.write(i)
        diffupdown_object.write('\n')

# enrichgo and kegg - gene name 
genenamefile = 'diff_gene.txt'
with open(genenamefile,'w') as genename_object:
    for i in geneset:
        genename_object.write(i)
        genename_object.write('\n')

updowngenefile = 'updown_gene.csv'
with open(updowngenefile,'w') as updown_object:
    for i in upgeneset:
        updown_object.write(i)
        updown_object.write(',up')
        updown_object.write('\n')
    for i in downgeneset:
        updown_object.write(i)
        updown_object.write(',down')
        updown_object.write('\n')

upfile = 'up_gene.txt'
with open(upfile,'w') as up_object:
    for i in upgeneset:
        up_object.write(i)
        up_object.write('\n')

downfile = 'down_gene.txt'
with open(downfile,'w') as down_object:
    for i in downgeneset:
        down_object.write(i)
        down_object.write('\n')


# heatmap file
heatmapfile = 'hm.csv'
with open(heatmapfile,'w') as hm_object:
    hm_object.write('gene,value1,value2\n')
    for i in heatmapset:
        hm_object.write(i)
        hm_object.write('\n')


# vol plot file
volfile = 'vol.csv'
with open(volfile,'w') as vol_object:
    vol_object.write('Log2FC,p_value,significant\n')
    for i in volset:
        vol_object.write(i)
        vol_object.write('\n')



print('finish')




