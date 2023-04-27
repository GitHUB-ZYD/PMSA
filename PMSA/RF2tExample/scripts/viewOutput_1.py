import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output_1', type=str, help='第一个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--savetxt', type=str, help='输出保存到该文件')
args = parser.parse_args()
with open(args.output_1,"rb") as f:
    x=pickle.load(f)
print(type(x))
"""
    x = [
        {"querySeqID":[{"speciesID":"seqID"},{"speciesID":"seqID"},{"speciesID":"seqID"}]},
        {"querySeqID":[{"speciesID":"seqID"},{"speciesID":"seqID"},{"speciesID":"seqID"}]},
        {"querySeqID":[{"speciesID":"seqID"},{"speciesID":"seqID"},{"speciesID":"seqID"}]}
    ]
"""
with open(args.savetxt,"w") as f:
    for i in x:
        queryseqid = list(i.keys())[0]
        f.write(queryseqid+":\n")
        for j in i[queryseqid]:
            speciesID=list(j.keys())[0]
            f.write("\t"+speciesID+": "+j[speciesID]+"\n")
    
# with open("output_2.txt","w") as f:
#     for i in x: # i是x的键
#         f.write(i+":\n")
#         for j in x[i]:
#             f.write("\t\t{}:{}\n".format(j,x[i][j]))
        

