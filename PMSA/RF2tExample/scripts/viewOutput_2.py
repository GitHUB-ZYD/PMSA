import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output_2', type=str, help='第二个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--savetxt', type=str, help='输出保存到该文件')
args = parser.parse_args()
with open(args.output_2,"rb") as f:
    x=pickle.load(f)
# 它的输出是这样的
# {
#    "queryseq1":{"specie1":"seq1","specie2":"seq2","specie3":"seq3"}
# }
with open(args.savetxt,"w") as f:
    for i in x: # i是x的键
        f.write(i+":\n")
        for j in x[i]:
            f.write("\t\t{}:{}\n".format(j,x[i][j]))
        

