# 这个脚本用于统计 每一对直系同源组共有的物种数目，用于构建pmsa
import pickle
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--output_3', type=str, help='第三个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--output_4_csv', type=str, help='输出每一对直系同源组的信息')
args = parser.parse_args()

with open(args.output_3,"rb") as f:
    x=pickle.load(f)

#   x是这样的[{},{},{}]
#   x[i]是这样的{"queryseqid":"seq","specie1id":"seq","specie2id":"seq",}

results=[]
for i in range(len(x)):
    for j in range(i+1,len(x)):
        keys_i = list(x[i].keys())
        keys_j = list(x[j].keys())
        query_i = keys_i[0]
        query_j = keys_j[0]
        query_i_len=len(x[i][query_i])
        query_j_len=len(x[j][query_j])
        species_i = set(keys_i[1:])
        species_j = set(keys_j[1:])
        if len(species_i)!=len(keys_i[1:]) or len(species_j)!=len(keys_j[1:]):
            print('错误')
            exit(0)
        intersection_count = len(species_i & species_j)
        result = [query_i,query_i_len,query_j,query_j_len,len(species_i),len(species_j),intersection_count]
        results.append(result)

df = pd.DataFrame(results,columns=['seq1','seq1_len','seq2','seq2_len','seq1_species','seq2_species','intersection'])
df.to_csv(args.output_4_csv,index=False)

        