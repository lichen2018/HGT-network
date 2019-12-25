import pandas as pd
import os
import fnmatch
import numpy as np
from collections import Counter



mother_sample_list = []
child_sample_list = []
fi = open('/disk2/workspace/lichen/li/infant/infant_cluster.hgt.bkp.txt', 'r')
M2C_bkp_label_dict = {}
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    bkp_info = buf.split(',')
    bkp_name = bkp_info[0]+','+bkp_info[1]+','+bkp_info[2]+','+bkp_info[3]
    label_list = []
    for i in range(4, len(bkp_info)):
        if bkp_info[i][5] == 'M':
            label_list.append('Mother')
            if bkp_info[i] not in mother_sample_list:
                mother_sample_list.append(bkp_info[i])
        else:
            label_list.append('Child')
            if bkp_info[i] not in child_sample_list:
                child_sample_list.append(bkp_info[i])
        #if label not in label_list:
    M2C_bkp_label_dict.update({bkp_name:label_list})

df = pd.read_csv('/disk2/workspace/lichen/li/infant/infant_hgt_hits_clubkp_cluster_0.6_10_clusterPartition_dynamic_tree_cut.csv')
M2C_cluster_label_dict = {}
for index, row in df.iterrows():
    bkp_info = row[0].split(';')
    bkp_name = bkp_info[0]+','+bkp_info[1]+','+bkp_info[2]+','+bkp_info[3]
    cluster_label = row['colorh1']
    if cluster_label == 0:
        continue
    if cluster_label not in M2C_cluster_label_dict:
        M2C_cluster_label_dict.update({cluster_label:M2C_bkp_label_dict[bkp_name]})
    else:
        M2C_cluster_label_dict[cluster_label].extend(M2C_bkp_label_dict[bkp_name])

df = pd.read_csv('/disk2/workspace/lichen/li/infant/infant_hgt_hits_clubkp_cluster_0.6_10_clusterPartition_dynamic_tree_cut.csv')
M2C_bkp_cluster_dict = {}
for index, row in df.iterrows():
    bkp_info = row[0].split(';')
    bkp_name = bkp_info[0]+','+bkp_info[1]+','+bkp_info[2]+','+bkp_info[3]
    cluster_label = row['colorh1']
    if cluster_label == 0:
        continue
    if cluster_label not in M2C_bkp_cluster_dict:
        M2C_bkp_cluster_dict.update({cluster_label:[bkp_name]})
    else:
        M2C_bkp_cluster_dict[cluster_label].append(bkp_name)


count_M2C_cluster_label_dict = {}
for key in M2C_cluster_label_dict:
    count_M2C_cluster_label_dict.update({key:Counter(M2C_cluster_label_dict[key])})

label_sn_dict = {}
for file_name in os.listdir("/disk2/workspace/lichen/li/IBD/new_bkp_result"):
    if fnmatch.fnmatch(file_name, "*.filter.acc.bkp.txt"):
        sn_id = file_name.split('.')[0]
        sn = file_name[:5]
        lb = sample_label_dict[sn]
        if lb not in label_sn_dict:
            label_sn_dict.update({lb:[sn_id]})
        elif sn_id not in label_sn_dict[lb]:
            label_sn_dict[lb].append(sn_id)

acc_phylum_dict = {}
fi = open('/disk2/workspace/lichen/li/hgt_network/acc2tax.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    acc = buf.split(';')[0]
    phylum = buf.split(';')[2].strip('.')
    acc_phylum_dict.update({acc:phylum})

acc_genus_dict = {}
fi = open('/disk2/workspace/lichen/li/hgt_network/acc2tax.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    acc = buf.split(';')[0]
    genus = buf.split(';')[-2]
    if genus == '':
        continue
    genus = genus.strip()
    acc_genus_dict.update({acc:genus})

acc_species_dict = {}
fi = open('/disk2/workspace/lichen/li/hgt_network/acc2tax.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    acc = buf.split(';')[0]
    species = buf.split(';')[-1]
    if species == '':
        continue
    acc_species_dict.update({acc:species})


M2C_phylum_cluster_dict = {}
M2C_genus_cluster_dict = {}
M2C_species_cluster_dict = {}
for key in M2C_bkp_cluster_dict:
    for bkp in M2C_bkp_cluster_dict[key]:
        bkp_info = bkp.split(',')
        ref1 = bkp_info[0]
        ref2 = bkp_info[2]
        if ref1 in acc_phylum_dict:
            ph1 = acc_phylum_dict[ref1]
            if key not in M2C_phylum_cluster_dict:
                M2C_phylum_cluster_dict.update({key:[ph1]})
            else:
                M2C_phylum_cluster_dict[key].append(ph1)
        if ref2 in acc_phylum_dict:
            ph2 = acc_phylum_dict[ref2]
            if key not in M2C_phylum_cluster_dict:
                M2C_phylum_cluster_dict.update({key:[ph2]})
            else:
                M2C_phylum_cluster_dict[key].append(ph2)
        if ref1 in acc_genus_dict:
            g1 = acc_genus_dict[ref1]
            if key not in M2C_genus_cluster_dict:
                M2C_genus_cluster_dict.update({key:[g1]})
            else:
                M2C_genus_cluster_dict[key].append(g1)
        if ref2 in acc_genus_dict:
            g2 = acc_genus_dict[ref2]
            if key not in M2C_genus_cluster_dict:
                M2C_genus_cluster_dict.update({key:[g2]})
            else:
                M2C_genus_cluster_dict[key].append(g2)
        if ref1 in acc_species_dict:
            s1 = acc_species_dict[ref1]
            if key not in M2C_species_cluster_dict:
                M2C_species_cluster_dict.update({key:[s1]})
            else:
                M2C_species_cluster_dict[key].append(s1)
        if ref2 in acc_species_dict:
            s2 = acc_species_dict[ref2]
            if key not in M2C_species_cluster_dict:
                M2C_species_cluster_dict.update({key:[s2]})
            else:
                M2C_species_cluster_dict[key].append(s2)

label_sn_dict={'non':38, 'M2C':109}
threshold = 1/2
fo = open('infant_bkp_clt_per.txt','w')
fo.write('label,'+'cluster,'+'genus,'+'percentage'+'\n')
c1 = 0
c2 = 0
ls = []
for cluster_key in count_M2C_cluster_label_dict:
    clt_label = ''
    comm_cluster_label = []
    rlt = count_M2C_cluster_label_dict[cluster_key]
    label_count = 0
    for key in rlt:
        if key == 'non':
            label_count += rlt[key]/label_sn_dict[key]
        else:
            label_count += rlt[key]/label_sn_dict['M2C']
    if 'non' in rlt and (rlt['non']/label_sn_dict['non'])/ label_count > threshold:
        clt_label = 'non'
        c1 += 1
    else:
        if 'non' not in rlt:
            clt_label = 'M2C'
            c2 += 1
        elif (label_count-(rlt['non']/label_sn_dict['non']))/ label_count > threshold: 
            clt_label = 'M2C'
            c2 += 1
    genus_list = M2C_genus_cluster_dict[cluster_key]
    for genus in Counter(genus_list):
        if genus not in ls and Counter(genus_list)[genus]/len(genus_list) > 0.1:
            ls.append(genus)
    if clt_label == 'non':
        for key in Counter(genus_list):
            fo.write('Non-IBD,'+str(cluster_key)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
    if clt_label == 'M2C':
        for key in Counter(genus_list):
            fo.write('IBD,'+str(cluster_key)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
fo.close()



fo = open('M2C_bkp_clt_species_1.txt','w')
fo.write('label,'+'cluster,'+'species,'+'percentage'+'\n')
ls = []
species_dict = {}
for cluster_key in count_M2C_cluster_label_dict:
    clt_label = ''
    comm_cluster_label = []
    rlt = count_M2C_cluster_label_dict[cluster_key]
    label_count = 0
    for key in rlt:
        if key == 'non':
            label_count += rlt[key]/label_sn_dict[key]
        else:
            label_count += rlt[key]/label_sn_dict['M2C']
    if 'non' in rlt and (rlt['non']/label_sn_dict['non'])/ label_count > threshold:
        clt_label = 'non'
    else:
        if 'non' not in rlt:
            clt_label = 'M2C'
        elif (label_count-(rlt['non']/label_sn_dict['non']))/ label_count > threshold: 
            clt_label = 'M2C'
    species_list = M2C_species_cluster_dict[cluster_key]
    for key in Counter(species_list):
        if Counter(species_list)[key]/len(species_list) < 0.05:
            continue
        if key not in ls:
            ls.append(key)
        if clt_label == 'non':
            fo.write('Non-IBD,'+str(cluster_key)+','+key+','+str(Counter(species_list)[key]/len(species_list))+'\n')
        if clt_label == 'M2C':
            fo.write('IBD,'+str(cluster_key)+','+key+','+str(Counter(species_list)[key]/len(species_list))+'\n')
fo.close()

from statistics import mean
ls = []
for key in species_dict:
    ls.append([key, mean(species_dict[key])])


fo = open('M2C_bkp_clt_genus.txt','w')
fo.write('label,'+'cluster,'+'genus,'+'percentage'+'\n')
ls = []
m_count = 0
c_count = 0
for cluster_key in count_M2C_cluster_label_dict:
    clt_label = ''
    comm_cluster_label = []
    rlt = count_M2C_cluster_label_dict[cluster_key]
    label_count = 0
    for key in rlt:
        if key == 'Child':
            label_count += rlt[key]/len(child_sample_list)
        else:
            label_count += rlt[key]/len(mother_sample_list)
    if 'Child' in rlt and (rlt['Child']/len(child_sample_list))/ label_count > threshold:
        clt_label = 'Child'
    else:
        if 'Child' not in rlt:
            clt_label = 'Mother'
        elif (label_count-(rlt['Child']/len(child_sample_list)))/ label_count > threshold: 
            clt_label = 'Mother'
    genus_list = M2C_genus_cluster_dict[cluster_key]
    other = 0
    for key in Counter(genus_list):
        if Counter(genus_list)[key]/len(genus_list) < 0.05:
            other += Counter(genus_list)[key]/len(genus_list)
            continue
        if key not in ls:
            ls.append(key)
        if clt_label == 'Child':
            fo.write('Child,'+str(c_count)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
        if clt_label == 'Mother':
            fo.write('Mother,'+str(m_count)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
    if other != 0:
        if clt_label == 'Child':
            fo.write('Child,'+str(c_count)+',Other,'+str(other)+'\n')
        if clt_label == 'Mother':
            fo.write('Mother,'+str(m_count)+',Other,'+str(other)+'\n')
    if clt_label == 'Child':
        c_count += 1
    else:
        m_count += 1
fo.close()


fi = open('/disk2/workspace/lichen/li/hgt_network/M2C_bkp_clt_genus.txt', 'r')
c_genus_dict = {}
m_genus_dict = {}
c_count = 0
m_count = 0
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    clt_info = buf.split(',')
    if clt_info[0] == 'label':
        continue
    if clt_info[0] == 'Child':
        c_count += float(clt_info[3])
        if clt_info[2] not in c_genus_dict:
            c_genus_dict.update({clt_info[2]:float(clt_info[3])})
        else:
            c_genus_dict[clt_info[2]] += float(clt_info[3])
    else:
        m_count += float(clt_info[3])
        if clt_info[2] not in m_genus_dict:
            m_genus_dict.update({clt_info[2]:float(clt_info[3])})
        else:
            m_genus_dict[clt_info[2]] += float(clt_info[3])

for key in c_genus_dict:
    c_genus_dict[key] = c_genus_dict[key]/c_count
for key in m_genus_dict:
    m_genus_dict[key] = m_genus_dict[key]/m_count

sorted_c_genus_dict = sorted(c_genus_dict.items(), key=lambda x: x[1], reverse=True)
sorted_m_genus_dict = sorted(m_genus_dict.items(), key=lambda x: x[1], reverse=True)


fo = open('M2C_bkp_clt_species.txt','w')
fo.write('label,'+'cluster,'+'species,'+'percentage'+'\n')
ls = []
m_count = 0
c_count = 0
for cluster_key in count_M2C_cluster_label_dict:
    clt_label = ''
    comm_cluster_label = []
    rlt = count_M2C_cluster_label_dict[cluster_key]
    label_count = 0
    for key in rlt:
        if key == 'Child':
            label_count += rlt[key]/len(child_sample_list)
        else:
            label_count += rlt[key]/len(mother_sample_list)
    if 'Child' in rlt and (rlt['Child']/len(child_sample_list))/ label_count > threshold:
        clt_label = 'Child'
    else:
        if 'Child' not in rlt:
            clt_label = 'Mother'
        elif (label_count-(rlt['Child']/len(child_sample_list)))/ label_count > threshold: 
            clt_label = 'Mother'
    genus_list = M2C_species_cluster_dict[cluster_key]
    other = 0
    for key in Counter(genus_list):
        if Counter(genus_list)[key]/len(genus_list) < 0.05:
            other += Counter(genus_list)[key]/len(genus_list)
            continue
        if key not in ls:
            ls.append(key)
        if clt_label == 'Child':
            fo.write('Child,'+str(c_count)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
        if clt_label == 'Mother':
            fo.write('Mother,'+str(m_count)+','+key+','+str(Counter(genus_list)[key]/len(genus_list))+'\n')
    if other != 0:
        if clt_label == 'Child':
            fo.write('Child,'+str(c_count)+',Other,'+str(other)+'\n')
        if clt_label == 'Mother':
            fo.write('Mother,'+str(m_count)+',Other,'+str(other)+'\n')
    if clt_label == 'Child':
        c_count += 1
    else:
        m_count += 1
fo.close()