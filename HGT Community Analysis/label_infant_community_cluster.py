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
        community_info_dict.update({comm_id:{'network_id':row['network_id'],'timepoint':row['timepoint'],'label':row['m_or_c'],'node_list':row['node_in_comm']}})
        if row['m_or_c'] == 'C' and row['network_id'] not in c_sample_list:
            c_sample_list.append(row['network_id'])
        if row['m_or_c'] == 'M' and row['network_id'] not in m_sample_list:
            m_sample_list.append(row['network_id'])
    label_sample_dict.update({'C':len(c_sample_list), 'M':len(m_sample_list)})
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


community_info_dict, label_sample_dict = getCommunityInfo('/disk2/workspace/lichen/li/infant/infant_comm_info.csv')
df = pd.read_csv('/disk2/workspace/lichen/li/IBD/infant_comm_cluster_power3_0.6_10_clusterPartition_dynamic_tree_cut.csv')
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
for cluster_key in comm_cluster_dict:
    comm_cluster_label = []
    for comm in comm_cluster_dict[cluster_key]:
        comm_cluster_label.append(comm[0])
    rlt = Counter(comm_cluster_label)
    label_count = 0
    for key in rlt:
        label_count += rlt[key]/label_sample_dict[key]
    if 'C' in rlt and (rlt['C']/label_sample_dict['C'])/ label_count > threshold:
        if 'C' not in leader_label_dict:
            leader_label_dict.update({'C':[cluster_key]})
        else:
            leader_label_dict['C'].append(cluster_key)
    else:
        if 'C' not in rlt:
            if 'M' not in leader_label_dict:
                leader_label_dict.update({'M':[cluster_key]})
            else:
                leader_label_dict['M'].append(cluster_key)
        elif (label_count-(rlt['C']/label_sample_dict['C']))/ label_count > threshold: 
            if 'M' not in leader_label_dict:
                leader_label_dict.update({'M':[cluster_key]})
            else:
                leader_label_dict['M'].append(cluster_key)
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
    count_dict = {}
    for key in Counter(phylum_list):
        count_dict.update({key:Counter(phylum_list)[key]/len(phylum_list)})
    comm_phylum_cluster_dict.update({cluster_key:count_dict})