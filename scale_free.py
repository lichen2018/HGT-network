from statistics import mean
from statistics import variance

import pandas as pd
import fnmatch
import os
import igraph as ig
import powerlaw
import math


IBD_graph_list = []
IBD_file_list = []
for file_name in os.listdir("/public/lichen/hgt_network"):
    if fnmatch.fnmatch(file_name, "*csv"):
        df = pd.read_csv("/public/lichen/hgt_network/"+file_name)
        ref_name_list = df.columns.values.tolist()
        if len(ref_name_list) < 100:
            continue
        IBD_file_list.append(file_name)
        event_value = []
        for i in range(1, len(ref_name_list)):
            ref_name = ref_name_list[i]
            row_list=df[ref_name].tolist()
            event_value.append(row_list)
        g = ig.Graph.Adjacency(event_value)
        IBD_graph_list.append(g)
IBD_alpha_list = []
compare_list = []
for g in IBD_graph_list:
    d= g.degree()
    result = powerlaw.Fit(d)
    p=result.distribution_compare('power_law', 'lognormal_positive')
    compare_list.append(p)
    IBD_alpha_list.append(result.power_law.alpha)

lognormal_compare_list = []
for g in IBD_graph_list:
    d= g.degree()
    result = powerlaw.Fit(d)
    p=result.distribution_compare('power_law', 'lognormal')
    lognormal_compare_list.append(p)


exp_compare_list = []
for g in IBD_graph_list:
    d= g.degree()
    result = powerlaw.Fit(d)
    p=result.distribution_compare('power_law', 'exponential',normalized_ratio=True)
    exp_compare_list.append(p)

wei_compare_list = []
for g in IBD_graph_list:
    d= g.degree()
    result = powerlaw.Fit(d)
    p=result.distribution_compare('power_law', 'stretched_exponential',normalized_ratio=True)
    wei_compare_list.append(p)



fo = open('ibd_cmp.txt','w')
fo.write('compare,result\n')
for ele in compare_list:
    if ele[0] < 0:
        fo.write('power_law vs lognormal_positive,<0\n')
    else:
        fo.write('power_law vs lognormal_positive,>0\n')
for ele in exp_compare_list:
    if ele[0] < 0:
        fo.write('power_law vs exponential,<0\n')
    else:
        fo.write('power_law vs exponential,>0\n')
for ele in wei_compare_list:
    if ele[0] < 0:
        fo.write('power_law vs Weibull,<0\n')
    else:
        fo.write('power_law vs Weibull,>0\n')
fo.close()


scale_cd = []
scale_uc = []
scale_non = []
fo = open("IBD_alpha.txt",'w')
for i in range(len(IBD_alpha_list)):
    sample_id = IBD_file_list[i].split('.')[0][:5]
    ele = IBD_alpha_list[i]
    if sample_label_dict[sample_id] == 'CD':
        fo.write(str(ele)+','+'CD'+'\n')
        scale_cd.append(ele)
    if sample_label_dict[sample_id] == 'UC':
        fo.write(str(ele)+','+'UC'+'\n')
        scale_uc.append(ele)
    if sample_label_dict[sample_id] == 'non':
        fo.write(str(ele)+','+'Non-IBD'+'\n')
        scale_non.append(ele)
fo.close()


graph_list = []
file_list = []
for file_name in os.listdir("/public/lichen/hgt_network/infant_hgt_network"):
    if fnmatch.fnmatch(file_name, "*csv"):
        df = pd.read_csv("/public/lichen/hgt_network/infant_hgt_network/"+file_name)
        ref_name_list = df.columns.values.tolist()
        if len(ref_name_list) < 100:
            continue
        file_list.append(file_name)
        event_value = []
        for i in range(1, len(ref_name_list)):
            ref_name = ref_name_list[i]
            row_list=df[ref_name].tolist()
            event_value.append(row_list)
        g = ig.Graph.Adjacency(event_value)
        graph_list.append(g)


alpha_list = []
label_list = []
infant_compare_list = []
infant_exp_compare_list = []
infant_compare_dict = {'power_law vs truncated_power_law':[],'power_law vs stretched_exponential':[]}
sigma_list = []
for i in range(len(file_list)):
    flag = file_list[i].split('_')[0][5]
    g = graph_list[i]
    d = g.degree()
    result = powerlaw.Fit(d)
    if math.isnan(result.power_law.alpha):
        continue
    sigma_list.append(result.power_law.sigma)
    p=result.distribution_compare('power_law', 'lognormal_positive')
    infant_compare_list.append(p[0])
    p=result.distribution_compare('power_law', 'exponential',normalized_ratio=True)
    infant_exp_compare_list.append(p[0])
    p=result.distribution_compare('power_law', 'stretched_exponential',normalized_ratio=True)
    infant_compare_dict['power_law vs stretched_exponential'].append(p[0])
    p=result.distribution_compare('power_law', 'truncated_power_law')
    infant_compare_dict['power_law vs truncated_power_law'].append(p[0])

fo = open('infant_powerlaw_fit_alpha.txt','w')
for i in range(len(file_list)):
    flag = file_list[i].split('_')[0][5]
    g = graph_list[i]
    d = g.degree()
    result = powerlaw.Fit(d)
    d_str = ''
    for i in d:
        d_str += str(i)+','
    d_str += str(result.xmin)+','
    d_str += str(result.power_law.alpha)
    fo.write(d_str+'\n')

fi = open('infant_powerlaw_fit_alpha.txt','r')
infant_degree_list = []
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    info_part = buf.split(',')
    d_tmp = info_part[:-2]
    d = []
    for ele in d_tmp:
        d.append(int(ele))
    d_array=np.asarray(d)
    x_min = float(info_part[-2])
    alpha = float(info_part[-1])
    infant_degree_list.append([d_array,x_min,alpha])

fo = open('infant_cmp.txt','w')
fo.write('compare,ratio\n')
for ele in infant_compare_list:
    if ele < 0:
        fo.write('power_law vs lognormal_positive,<0\n')
    else:
        fo.write('power_law vs lognormal_positive,>0\n')
for ele in infant_exp_compare_list:
    if ele < 0:
        fo.write('power_law vs exponential,<0\n')
    else:
        fo.write('power_law vs exponential,>0\n')
fo.close()
fo = open('infant_cmp.txt','a')
for ele in infant_compare_dict['power_law vs stretched_exponential']:
    if ele < 0:
        fo.write('power_law vs Weibull,<0\n')
    else:
        fo.write('power_law vs Weibull,>0\n')
fo.close()