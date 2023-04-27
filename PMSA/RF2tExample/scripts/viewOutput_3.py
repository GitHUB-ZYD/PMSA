import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output_3', type=str, help='第三个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--savetxt', type=str, help='输出保存到该文件')
args = parser.parse_args()
with open(args.output_3,"rb") as f:
    x=pickle.load(f)
#   x是这样的[{},{},{}]
#   x[i]是这样的{"queryseqid":"seq","specie1id":"seq","specie2id":"seq",}
#

with open(args.savetxt,"w") as f:
    for i in x:
        for j in i: #j 是key
            f.write("{}:\n{}\n".format(j,i[j]))
        f.write("\n\n\n")
        

        

