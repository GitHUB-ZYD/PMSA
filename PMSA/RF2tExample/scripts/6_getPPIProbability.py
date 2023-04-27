import subprocess
import numpy as np
import pandas as pd
import uuid 
import os
import argparse

# 用于矫正
def convert_value(x, row_sum, col_sum, total_sum):
    if pd.isna(x):
        return x
    elif x >= 0.9:
        return x + x - (row_sum * col_sum) / total_sum
    else:
        return x - (row_sum * col_sum) / total_sum

parser = argparse.ArgumentParser()
parser.add_argument('--input_csv', type=str, help='要计算的ppi的信息文件')
parser.add_argument('--RTF2track', type=str, help='RoseTTAFold的predict_msa.py的位置')
parser.add_argument('--output_3', type=str, help='第三个脚本输出的比对后的直系同源组文件，pickle')
parser.add_argument('--save_ppi', type=str, help='保存RosseTTAFold计算后的未矫正的互作概率，每成功运行就会flush，即使程序失败也能保存计算好的结果')
parser.add_argument('--save_ppi_correction', type=str, help='保存矫正后的互作概率矩阵')
parser.add_argument('--cpu', type=int, help='0:使用gpu计算;1:使用cpu计算')
args = parser.parse_args()
# df = pd.read_csv("./run_long.csv",header=0)
df = pd.read_csv(args.input_csv,header=0)
#用于保存数据转化为矩阵
seq1s=[]
seq2s=[]
probs=[]
# with open("/ppi/qkj/ppi_propertity_long.csv","w") as f:
with open(args.save_ppi,"w") as f:
    f.write("seq1,seq2,property\n")
    
    for i in range(len(df)):
        aln = "/dev/shm/tmp/"+str(uuid.uuid4())+".aln"
        npz = "/dev/shm/tmp/"+str(uuid.uuid4())+".npz"
        seq1 = str(df.loc[i,'seq1'])
        print("seq1:{}".format(seq1))
        seq1s.append(seq1)
        seq2 = str(df.loc[i,'seq2'])
        print("seq2:{}".format(seq2))
        seq2s.append(seq2)
        L1 = df.loc[i,'seq1_len']
        cmd = "time python scripts/5_getPMSA.py --output_3 {} --seq1 {} --seq2 {} --output_5_aln {}".format(args.output_3,seq1,seq2,aln)
        subprocess.run(cmd,shell=True)
        if args.cpu==0:
            cmd = "time python {} -msa {} -npz {} -L1 {}".format(args.RTF2track,aln,npz,L1)
        elif args.cpu==1:
            cmd = "time python {} -msa {} -npz {} -L1 {} --cpu".format(args.RTF2track,aln,npz,L1)
        else:
            print("请指定cpu或者gpu计算")
            os.remove(aln)
            os.remove(npz)
            exit(0)
        subprocess.run(cmd,shell=True)
        data = np.load(npz)
        dist = pd.DataFrame(data['dist'])
        # 第一次矫正，去除最后10行和前10列
        dist=dist.iloc[:-10, 10:]
        dist.columns = range(dist.shape[1])
        max_value = dist.max().max()
        probs.append(max_value)
        f.write("{},{},{}\n".format(seq1,seq2,max_value))
        f.flush()
        os.remove(aln)
        os.remove(npz)
        print("完成度:{}/{}".format(i+1,len(df)))
        

# 将所有结果转化为矩阵，并进行矫正
data = {"rows":seq1s,"cols":seq2s,"values":probs}
df = pd.DataFrame(data).pivot(index="rows", columns="cols", values="values")
# 计算总和，跳过空值
total_sum = df.sum().sum(skipna=True)
# 对每一行应用转换函数
for index, row in df.iterrows():
    row_sum = row.sum(skipna=True)
    for col in df.columns:
        col_sum = df[col].sum(skipna=True)
        df.loc[index, col] = row.apply(convert_value, args=(row_sum, col_sum, total_sum))[col]

df.to_csv(args.save_ppi_correction,index=True, header=True)