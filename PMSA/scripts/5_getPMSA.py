#从run_4.py中可以得到哪些的序列比较多，现在我们将传入比对后的结果和要比对的序列id,生成pmsa文件

import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output_3', type=str, help='第三个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--seq1', type=str, help='用于合并比对的第一条蛋白')
parser.add_argument('--seq2', type=str, help='用于合并比对的第二条蛋白')
parser.add_argument('--output_5_aln', type=str, help='合并后的比对文件')
args = parser.parse_args()
with open(args.output_3,"rb") as f:
    x=pickle.load(f)
#   x是这样的[{},{},{}]
#   x[i]是这样的{"queryseqid":"seq","specie1id":"seq","specie2id":"seq",}
seq1 = args.seq1
seq2 = args.seq2
# pairs = [d for d in x if list(d.keys())[0] == seq1 or list(d.keys())[0] == seq2]
oth1 = [d for d in x if list(d.keys())[0] == seq1][0]
oth2 = [d for d in x if list(d.keys())[0] == seq2][0]
print(type(oth1))
intersection = list(set(oth1) & set(oth2))

with open(args.output_5_aln,"w") as f:
        # i 是交集的键的名字
        l1=">"+seq1+"_"+seq2+"\n"
        l2=oth1[seq1]+oth2[seq2]+"\n"
        f.write(l1)
        f.write(l2)
        for i in intersection:
            q1id=i+"_1"
            q2id=i+"_2"
            l1=">"+q1id.split('/')[-1]+"_"+q2id.split('/')[-1]+"\n"
            l2=oth1[i]+oth2[i]+"\n"
            f.write(l1)
            f.write(l2)

# with open("q1_q2.fasta","w") as f:
#         q1keys=list(x[seq1].keys())
#         q2keys=list(x[seq2].keys())
#         q1id=q1keys[0]
#         q1seq=x[seq1][q1id]
#         print(len(q1seq))
#         q2id=q2keys[0]
#         q2seq=x[seq2][q2id]
#         print(len(q2seq))
#         l1=">"+q1id+"_"+q2id+"\n"
#         l2=q1seq+q2seq+"\n"
#         f.write(l1)
#         f.write(l2)
#         for k in intersection:
#             q1id=k+"_1"
#             q2id=k+"_2"
#             q1seq=x[i][k]
#             q2seq=x[j][k]
#             l1=">"+q1id.split('/')[-1]+"_"+q2id.split('/')[-1]+"\n"
#             l2=q1seq+q2seq+"\n"
#             f.write(l1)
#             f.write(l2)