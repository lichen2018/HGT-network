import os

import numpy as np
from collections import Counter
import pandas as pd


def getCommunityInfo(file_name):
    c_sample_list = []
    m_sample_list = []
    label_sample_dict = {}
    community_info_dict = {}
    df = pd.read_csv(file_name)
    for index, row in df.iterrows():
        comm_id = row['comm_id']
        community_info_dict.update({comm_id:{'network_id':row['network_id'],'individual_id':row['individual_id'],'timepoint':row['timepoint'],'label':row['diseaseState'],'node_list':row['node_in_comm']}})
    label_sample_dict.update({'non':38, 'ibd':109})
    return community_info_dict, label_sample_dict


acc_phylum_dict = {}
fi = open('C:\\Users\\cli335\\Desktop\\hgt_network\\acc2tax.txt', 'r')
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    acc = buf.split(';')[0]
    phylum = buf.split(';')[2].strip('.')
    acc_phylum_dict.update({acc:phylum})


community_info_dict, label_sample_dict = getCommunityInfo('C:\\Users\\cli335\\Desktop\\hgt_network\\comm_info.csv')
df = pd.read_csv('C:\\Users\\cli335\\Desktop\\hgt_network\\comm_cluster_power3_0.6_10_clusterPartition_dynamic_tree_cut.csv')
comm_cluster_dict = {}
for index, row in df.iterrows():
    comm_name = row[0]
    cluster_label = row['colorh1']
    if cluster_label == 0:
        continue
    if cluster_label not in comm_cluster_dict:
        comm_cluster_dict.update({cluster_label:[[community_info_dict[comm_name]['label'],community_info_dict[comm_name]['network_id'],community_info_dict[comm_name]['node_list']]]})
    else:
        comm_cluster_dict[cluster_label].append([community_info_dict[comm_name]['label'],community_info_dict[comm_name]['network_id'],community_info_dict[comm_name]['node_list']])

comm_ls = {}
leader_label_dict = {}
threshold = 1/2
comm_phylum_cluster_dict = {}
fo = open('group_ibd_comm_clt_per1.txt','w')
fo.write('label,'+'cluster,'+'phylum,'+'percentage'+'\n')
ibd_count = 0
non_count = 0
for cluster_key in comm_cluster_dict:
    comm_cluster_label = []
    clt_label = ''
    for comm in comm_cluster_dict[cluster_key]:
        comm_cluster_label.append(comm[0])
    rlt = Counter(comm_cluster_label)
    label_count = 0
    for key in rlt:
        label_count += rlt[key]/label_sample_dict[key]
    if 'non' in rlt and (rlt['non']/len(label_sample_dict['non']))/ label_count > threshold:
        clt_label = 'non'
    else:
        if 'non' not in rlt:
            clt_label = 'ibd'
        elif (label_count-(rlt['non']/len(label_sample_dict['non'])))/ label_count > threshold: 
            clt_label = 'ibd'
    comm_ls.update({cluster_key:rlt})
    phylum_list = []
    for comm in comm_cluster_dict[cluster_key]:
        comm_label = comm[0]
        ref_name_list = comm[-1].split(';')[:-1]
        for ref_name in ref_name_list:
            if ref_name not in acc_phylum_dict:
                continue
            if acc_phylum_dict[ref_name] == '':
                continue
            phylum_name = acc_phylum_dict[ref_name]
            if acc_phylum_dict[ref_name] == 'ORGANISM  [Clostridium] saccharolyticu' or acc_phylum_dict[ref_name] == 'ORGANISM Clostridium saccharolyticu':
                phylum_name = 'Clostridium saccharolyticu'
            phylum_list.append(phylum_name)
    if clt_label == 'non':
        for key in Counter(phylum_list):
            fo.write('Non-IBD,'+str(non_count)+','+key+','+str(Counter(phylum_list)[key]/len(phylum_list))+'\n')
        non_count += 1
    if clt_label == 'ibd':
        for key in Counter(phylum_list):
            fo.write('IBD,'+str(ibd_count)+','+key+','+str(Counter(phylum_list)[key]/len(phylum_list))+'\n')
        ibd_count += 1
fo.close()
