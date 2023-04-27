import subprocess
import pickle
import argparse
import uuid
import os
from concurrent.futures import ThreadPoolExecutor
#获取某个fasta文件中所有序列id
def getEcoliSeqID(EcoliFilePath):
    cmd = 'grep "^>" {} | cut -d " " -f 1 | tr -d ">"'.format(EcoliFilePath)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # 获取标准输出和标准错误输出
    out, err = p.communicate()

    # 将输出按行分割并保存到列表
    output_lines = out.decode().split('\n')
    error_lines = err.decode().split('\n')
    #去除空字符串
    output_lines = [x.strip() for x in output_lines if x.strip()]
    return output_lines

# 根据fasta文件和序列id获取序列
def getseq(fasta,seqid):
    cmd = "samtools faidx {} '{}'".format(fasta,seqid)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # 获取标准输出和标准错误输出
    out, err = p.communicate()
    # 将输出按行分割并保存到列表
    seq = out.decode("utf-8")
    # print("decode的seq:{}".format(seq))
    # print("split的：{}".format(str(seq.split("\n")[1:])))
    seq = "".join(seq.split("\n")[1:])
    # print("getseq的seq:{}".format(seq))
    return seq

def getseq_aln(fasta,seqid):
    with open(fasta,"r") as f:
        seq=""
        x=f.readlines()
        seqcount=0
        for i in x:
            if seqcount==2:
                return seq
            if seqcount==1:
                if i.startswith(">"):
                    seqcount=seqcount+1
                else:
                    seq=seq+i.rstrip()
            if seqid in i:
                seqcount=seqcount+1
    return seq


# 对一条ecoli蛋白的直系同源组进行比对处理，从{"queryseqid1":{"物种1id":"seq1","物种2id":"seq2","物种3id":"seq3"},"queryseqid2":{"物种1id":"seq1","物种2id":"seq2"}} 提取出"queryseqid1":{"物种1id":"seq1","物种2id":"seq2","物种3id":"seq3"} 进行处理
def dealwith_one(queryseqid,queryseq,uniques_all):#queryseqid是当前处理的ecoli蛋白的id,queryseq是它的序列；uniques_all就是所有的直系同源组字典
    # temp_alignment="./aln_0.fasta"
    aln_0 = "/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    aln_1 = "/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    aln_2 = "/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    aln_3 = "/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    err = "/dev/shm/tmp/"+str(uuid.uuid4())+".txt"
    # 当前处理的ecoli蛋白的直系同源组
    seqs=uniques_all[queryseqid]
    # 将所有待比对的序列输出到aln_0.fasta中，包括queryseq和其他物种的seq，queryseq放到第一条
    with open(aln_0,"w") as f:
        seqid = ">"+queryseqid+"\n"
        seq = queryseq+"\n"
        f.write(seqid)
        f.write(seq)
        for i in seqs:
            seqid=">"+i+"\n"
            seq=seqs[i]+"\n"
            f.write(seqid)
            f.write(seq)
    # mafft进行比对，输出到aln_1.fasta
    cmd = "mafft --thread 13 --localpair --inputorder {} > {} 2>{}".format(aln_0,aln_1,err)
    subprocess.run(cmd,shell=True)
    # print(cmd)
    # 从aln_1中提取第一条蛋白，用于去除第一条蛋白上所有的gap
    # 有个bug，mafft生成的aln_1.fasta不要用getseq()方法获取序列，因为末尾会得到其他的东西
    # seq = getseq("aln_1.fasta",queryseqid)
    # print("从aln_1提取第一条的seq:{}".format(seq))
    with open(aln_1,"r") as f:
        seq=""
        x=f.readlines()
        seqcount=0
        for i in x:
            if seqcount==2:
                break
            if i.startswith(">"):
                seqcount=seqcount+1
            else:
                seq=seq+i.rstrip()
    # print("从aln_1提取第一条的seq:{}".format(seq))
    # 记录gap位置
    pos=[]
    # print("pos:\n{}".format(str(pos)))
    for i in range(len(seq)):
        if seq[i]=='-':
            pos.append(i)
    # print("gap的位置：{}".format(str(pos)))
    # print("用于排除gap的seq:{}\n111".format(seq))
    # print("len(pos):{}".format(len(pos)))
    if len(pos)==0:
        # print("没有gap的seq:{}".format(seq))
        # print("进行判断 if")
        cmd = "cat {} > {}".format(aln_1,aln_2)
    else:
        # print("进行判断 else")
        cmd = "trimal -in {} -out {} -selectcols {} ".format(aln_1,aln_2,"{")
        for i in range(len(pos)):
            if i==len(pos)-1:
                cmd=cmd+str(pos[i])+" }"
            else:
                cmd=cmd+str(pos[i])+","
    # print(cmd)
    subprocess.run(cmd,shell=True)
    #移除gap fraction大于0.5的序列
    cmd="trimal -in {} -out {} -gt 0.5".format(aln_2,aln_3)
    subprocess.run(cmd,shell=True)
    seqids=getEcoliSeqID(aln_3)
    result={}
    for i in seqids:
        result[i]=getseq_aln(aln_3,i)
    # print(result)
    os.remove(aln_0)
    os.remove(aln_1)
    os.remove(aln_2)
    os.remove(aln_3)
    os.remove(err)
    # print("dealwithone的result:\n{}".format(result))
    return result

def work(queryProt,i,results,x,j):
    queryseq = getseq(queryProt,i)
    # print("第一次调用queryseq得到的:{}".format(queryseq))
    # print("当前seqid:{}\n{}".format(i,queryseq))
    result = dealwith_one(i,queryseq,x)
    print("{}/{}".format(j,len(x)))
    # print(result)
    results.append(result)
    
    


parser = argparse.ArgumentParser()
parser.add_argument('--output_2', type=str, help='第二个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--queryProteome', type=str, help='查询的蛋白组 或者 fasta序列文件的文件地址')
parser.add_argument('--output_3', type=str, help='当前脚本输出的比对后的pickle文件')
parser.add_argument('--num_thread', type=int, help='最多多少个线程一起运行')

args = parser.parse_args()

with open(args.output_2,"rb") as f:
    x=pickle.load(f)
    #x中保存的是 {“ecoli1seqid”:{specie1id:seq,specie2id:seq,...},“ecoli2seqid”:{specie1id:seq,specie2id:seq,...},}}
# 查询的蛋白文件
queryProt=args.queryProteome
results = []
executor = ThreadPoolExecutor(max_workers=args.num_thread)
tasks = []
j=0
for i in x:
    j=j+1
    future = executor.submit(work,queryProt,i,results,x,j)
    tasks.append(future)
    # queryseq = getseq(queryProt,i)
    # results.append(dealwith_one(i,queryseq,x))
    # j=j+1
    # print("{}/{}".format(j,len(x)))
executor.shutdown(wait=True)
# 所有任务完成之后执行其他代码
print("比对完成，开始输出到文件")
# print(results)

with open(args.output_3,"wb") as f:
    pickle.dump(results,f)