import os
from Bio import SeqIO
from Bio import Entrez
import sys


'''fi = open('/mnt/disk2_workspace/lichen/li/IBD/ibd_cluster.supp.geneID.txt', 'r')
processed_list = []
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    bkp_info_part = buf.split(',')
    bkp = bkp_info_part[0]+','+bkp_info_part[1]+','+bkp_info_part[2]+','+bkp_info_part[3]
    processed_list.append(bkp)
    #processed_dict.update({bkp:1})
fi.close()'''

gene_list = []
Entrez.email = 'whatever@mail.com'
#hgt_cluster_file_list = ['/mnt/disk2_workspace/lichen/li/infant/infant_cluster.hgt.bkp.txt']

'''ibd_ref_list = []
fi = open('/mnt/disk2_workspace/lichen/li/IBD/ibd_cluster.hgt.bkp.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    bkp_info_part = buf.split(',')
    ref_1 = bkp_info_part[0]
    ref_2 = bkp_info_part[2]
    if ref_1 not in ibd_ref_list:
        ibd_ref_list.append(ref_1)
    if ref_2 not in ibd_ref_list:
        ibd_ref_list.append(ref_2)   

found = False
ref_list = []
for hgt_cluster_file in hgt_cluster_file_list:
    fi = open(hgt_cluster_file, 'r')
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        bkp_info_part = buf.split(',')
        ref_1 = bkp_info_part[0]
        ref_2 = bkp_info_part[2]
        if ref_1 not in ref_list and ref_1 not in ibd_ref_list:
            ref_list.append(ref_1)
            gbk_1 = '/mnt/disk2_workspace/lichen/li/genebank/'+ref_1+'.gbk'
            if os.path.exists(gbk_1) is False:
                handle = Entrez.efetch(db='nucleotide', id=ref_1, rettype='gb')
                local_file = open(gbk_1, 'w')
                local_file.write(handle.read())
                handle.close()
                local_file.close()
            for seq_record in SeqIO.parse(gbk_1, "genbank"):
                for seq_feature in seq_record.features:
                    if 'db_xref' in seq_feature.qualifiers:
                        for val in seq_feature.qualifiers['db_xref']:
                            info_part = val.split(':')
                            if info_part[0] == 'GeneID' or info_part[0] == 'GI':
                                feature_start = str(seq_feature.location.start)
                                feature_end = str(seq_feature.location.end)
                                fo = open('/mnt/disk2_workspace/lichen/li/ref_gene_region.txt','a')
                                fo.write(ref_1+','+val+','+feature_start+','+feature_end+'\n')
                                fo.close()
        if ref_2 not in ref_list and ref_2 not in ibd_ref_list:
            ref_list.append(ref_2)
            gbk_2 = '/mnt/disk2_workspace/lichen/li/genebank/'+ref_2+'.gbk'
            if os.path.exists(gbk_2) is False:
                handle = Entrez.efetch(db='nucleotide', id=ref_2, rettype='gb')
                local_file = open(gbk_2, 'w')
                local_file.write(handle.read())
                handle.close()
                local_file.close()
            for seq_record in SeqIO.parse(gbk_2, "genbank"):
                for seq_feature in seq_record.features:
                    if 'db_xref' in seq_feature.qualifiers:
                        for val in seq_feature.qualifiers['db_xref']:
                            info_part = val.split(':')
                            if info_part[0] == 'GeneID' or info_part[0] == 'GI':
                                feature_start = str(seq_feature.location.start)
                                feature_end = str(seq_feature.location.end)
                                fo = open('/mnt/disk2_workspace/lichen/li/ref_gene_region.txt','a')
                                fo.write(ref_2+','+val+','+feature_start+','+feature_end+'\n')
                                fo.close()'''

acc_species_dict = {}
fi = open('/mnt/disk2_workspace/lichen/li/hgt_network/acc2tax.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    acc = buf.split(';')[0]
    species = buf.split(';')[-1]
    acc_species_dict.update({acc:species})

ref_gene_dict = {}
fi = open('/mnt/disk2_workspace/lichen/li/ref_gene_region.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    gene_info_part = buf.split(',')
    ref = gene_info_part[0]
    gene_name = gene_info_part[1]
    start = gene_info_part[2]
    end = gene_info_part[3]
    if ref not in ref_gene_dict:
        ref_gene_dict.update({ref:[[gene_name, start, end]]})
    else:
        if [gene_name, start, end] not in ref_gene_dict[ref]:
            ref_gene_dict[ref].append([gene_name, start, end])


hgt_cluster_file_list = ['/mnt/disk2_workspace/lichen/li/IBD/ibd_cluster.hgt.bkp.txt']
found = False
bkp_gene_dict = {}
for hgt_cluster_file in hgt_cluster_file_list:
    fi = open(hgt_cluster_file, 'r')
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        bkp_info_part = buf.split(',')
        ref_1 = bkp_info_part[0]
        pos_1 = int(bkp_info_part[1])
        ref_2 = bkp_info_part[2] 
        pos_2 = int(bkp_info_part[3])
        sample_list = buf[len(ref_1)+len(ref_2)+len(bkp_info_part[1])+len(bkp_info_part[3])+4:]
        if ref_1 not in ref_gene_dict or ref_2 not in ref_gene_dict:
            continue
        if ref_1 not in acc_species_dict or ref_2 not in acc_species_dict:
            continue
        if acc_species_dict[ref_1] == acc_species_dict[ref_2]:
            continue
        flag_1 = False
        flag_2 = False
        gene1_list = []
        if ref_1 in bkp_gene_dict:
            if pos_1 in bkp_gene_dict[ref_1]:
                gene1_list = bkp_gene_dict[ref_1][pos_1]
                flag_1 = True
        if len(gene1_list) == 0:
            for gene_info in ref_gene_dict[ref_1]:
                if '>' in gene_info[1] or '<' in gene_info[1]:
                    gene_info[1] = gene_info[1][1:]
                if '>' in gene_info[2] or '<' in gene_info[2]:
                    gene_info[2] = gene_info[2][1:]
                feature_start = int(gene_info[1])
                feature_end = int(gene_info[2])
                if pos_1 < feature_start:
                    break
                if pos_1 >= feature_start and pos_1 <= feature_end:
                    if gene_info[0]+','+str(feature_start)+','+str(feature_end) not in gene1_list:
                        gene1_list.append(gene_info[0]+','+str(feature_start)+','+str(feature_end))
                        flag_1 = True
            if flag_1 == True:
                if ref_1 not in bkp_gene_dict:
                    bkp_gene_dict.update({ref_1:{pos_1:gene1_list}})
                else:
                    if pos_1 not in bkp_gene_dict[ref_1]:
                        bkp_gene_dict[ref_1].update({pos_1:gene1_list})
        if flag_1 == False:
            continue
        gene2_list = []
        if ref_2 in bkp_gene_dict:
            if pos_2 in bkp_gene_dict[ref_2]:
                gene2_list = bkp_gene_dict[ref_2][pos_2]
                flag_2 = True
        if len(gene2_list) == 0:
            for gene_info in ref_gene_dict[ref_2]:
                if '>' in gene_info[1] or '<' in gene_info[1]:
                    gene_info[1] = gene_info[1][1:]
                if '>' in gene_info[2] or '<' in gene_info[2]:
                    gene_info[2] = gene_info[2][1:]             
                feature_start = int(gene_info[1])
                feature_end = int(gene_info[2])
                if pos_2 < feature_start:
                    break
                if pos_2 >= feature_start and pos_2 <= feature_end:
                    if gene_info[0]+','+str(feature_start)+','+str(feature_end) not in gene2_list:
                        gene2_list.append(gene_info[0]+','+str(feature_start)+','+str(feature_end))
                        flag_2 = True
            if len(gene2_list) != 0:
                if ref_2 not in bkp_gene_dict:
                    bkp_gene_dict.update({ref_2:{pos_2:gene2_list}})
                else:
                    if pos_2 not in bkp_gene_dict[ref_2]:
                        bkp_gene_dict[ref_2].update({pos_2:gene2_list})
        if flag_1 is True and flag_2 is True:
            outfile = hgt_cluster_file.split('.')[0]+'.fusion.geneID.fast1.txt'
            geneID = ''
            if len(gene1_list) > 0:
                geneID += 'from;'
                for ele in gene1_list:
                    geneID += ele
                    geneID += ';'
            if len(gene2_list) > 0:
                geneID += 'to;'
                for ele in gene2_list:
                    geneID += ele
                    geneID += ';'
            if len(geneID) > 0:
                fo = open(outfile, 'a')
                fo.write(bkp_info_part[0]+','+bkp_info_part[1]+','+bkp_info_part[2]+','+bkp_info_part[3]+';'+geneID[:-1]+';'+sample_list+'\n')
                fo.close()
