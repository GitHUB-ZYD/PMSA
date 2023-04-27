import pickle
import subprocess
import uuid
import os
import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor
# 1. 读取数据
def getAllCandidateOth(filepath="./Ecoli_2_specie_1w.pickle"):
    with open(filepath,"rb") as f:
        x=pickle.load(f)
    return x

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
    seq = "".join(seq.split("\n")[1:])
    return seq

# 2. 传入一个蛋白的直系同源的列表添加seq信息，就在原来的数组上进行操作
def appendSeqs_one(oth_one_dict,count):
    ecoliseqid = list(oth_one_dict.keys())[0]
    for i in range(len(oth_one_dict[ecoliseqid])):
        specieid = list(oth_one_dict[ecoliseqid][i].keys())[0]
        seqid = oth_one_dict[ecoliseqid][i][specieid]
        seq = getseq(specieid,seqid)
        temp = {seqid:seq}
        oth_one_dict[ecoliseqid][i][specieid]=temp
    print("appendSeqs_{}:finished\n".format(count))
        
# 3. 传入所有ecoli的直系同源蛋白，添加seq信息
def appendSeqs_All(oth_all,num_thread):
    executor = ThreadPoolExecutor(max_workers=num_thread)
    tasks = []
    for i in range(len(oth_all)):
        future = executor.submit(appendSeqs_one, oth_all[i],i+1)
        tasks.append(future)
        # 关闭线程池并等待所有任务完成
    executor.shutdown(wait=True)
        # 所有任务完成之后执行其他代码
    print("已经成功添加所有seq序列")
        # appendSeqs_one(oth_one_dict=oth_all[i])
        # print("添加seq信息：{}/{}".format(i,len(oth_all)))  
    #     p = threading.Thread(target=appendSeqs_one,args=(oth_all[i],))
    #     p.start()
    #     process.append(p)
    #     print("添加seq信息：{}/{}".format(i+1,len(oth_all)))  
    # for i in process:
    #     i.join()
    # print("已经成功添加所有seq序列")
    
        

# 4. 将重复的和不重复的分开
def seperate(oth_one_dict,duplicates,uniques):
    # 一个直系同源组的列表
    oth_one=oth_one_dict[list(oth_one_dict.keys())[0]]
    all_keys = []
    duplicate_keys = set()
    for d in oth_one:
        keys = d.keys()
        for k in keys:
            if k in all_keys:
                duplicate_keys.add(k)
            else:
                all_keys.append(k)
    # 将重复的字典移动到另一个列表
    # duplicates = []
    # uniques = []
    for d in oth_one:
        keys = d.keys()
        for k in keys:
            if k in duplicate_keys:
                duplicates.append(d)
            else:
                uniques.append(d)
    # print("这里是duplicate:{}".format(duplicates))
    # unique的格式[{"speciesID1":{'seqid':'seq'}},{"speciesID2":{'seqid':'seq'}},{"speciesID3":{'seqid':'seq'}} ]
    # duplicates格式类似，不过，每一个speciesID是有重复的，比如[{"speciesID1":{'seqid':'seq'}},{"speciesID1":{'seqid':'seq'}},{"speciesID3":{'seqid':'seq'}},{"speciesID3":{'seqid':'seq'}} ] 
    #


# 5. duplicates按照同一个物种的序列进行分组
def duplicateseperate(duplicates):
    output_dict = {}
    for d in duplicates:
        for key, value in d.items():
            if key not in output_dict:
                output_dict[key] = []
            output_dict[key].append(value)
    # 遍历分组后的重复，一个直系同源组下的一个重复物种的序列id及其序列信息
    # keys = list(output_dict.keys())
    # for i in range(len(keys)):
    #     print(output_dict[keys[i]])  
    return output_dict  #{'speciesid1':[{'seqid1':"seq"},{'seqid2':"seq"},{'seqid3':'seq'}],
                        #  'speciesid2':[{'seqid1':"seq"},{'seqid2':"seq"},{'seqid3':'seq'}]      }

# 6. cd-hit去重复
#  cd-hit -i cdhittest.fasta -o new.fasta -c 0.95 -aS 0.8 -d 0
#传入的是duplicates分组以后的一个物种的序列集  [{seqid1:seq1},{seqid2:seq2}]
def runCdHit(duplicates_one):
    temp_input = str(uuid.uuid4())+"_in.fasta"
    temp_output= str(uuid.uuid4())+"_out.fasta"
    with open(temp_input,"w")as f:
        for i in duplicates_one:
            key = list(i.keys())[0]
            seqid=">"+key+"\n"
            seq=i[key]+"\n"
            f.write(seqid)
            f.write(seq)
    cmd="cd-hit -i {} -o {} -T 14 -c 0.95 -aS 0.8 -d 0 1>/dev/null".format(temp_input,temp_output)
    subprocess.run(cmd,shell=True)
    aftercdhit=getEcoliSeqID(temp_output)
    print("aftercdhit:\n"+str(aftercdhit))
    os.remove(temp_input)
    os.remove(temp_output)
    os.remove(temp_output+".clstr")
    duplicates_one = [d for d in duplicates_one if list(d.keys())[0] in aftercdhit]
    return duplicates_one


# 7. 检验序列是否能够合并 只能aftercdhit后是两条的可以使用
# queryfile是对ecoli序列遍历的时候生成的一个一条序列的fasta文件
# mergeloc是一个字典，保存左右的信息
def ifmerge(queryfastafile,duplicate_one_aftercdhit,mergeloc):
    # temp_queryfile = str(uuid.uuid4())+"_query.fasta"
    temp_sequences = "/dev/shm/tmp/"+str(uuid.uuid4())+"_sequence.fasta"
    temp_output= "/dev/shm/tmp/"+str(uuid.uuid4())+"_output.csv"
    with open(temp_sequences,"w") as f:
        for i in duplicate_one_aftercdhit:
            key = list(i.keys())[0]
            seqid=">"+key+"\n"
            seq=i[key]+"\n"
            f.write(seqid)
            f.write(seq)
    cmd="blastp -query {} -subject {} -outfmt '6 delim=, qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident ppos positive' -evalue 0.01 -num_threads=14 -out {}".format(queryfastafile,temp_sequences,temp_output)
    subprocess.run(cmd,shell=True)
    columns = ["qseqid","qlen","sseqid","slen","qstart","qend","sstart","send","evalue","length","pident","nident","ppos","positive"]
    df = pd.read_csv(temp_output, sep=',', header=None, names=columns)
    print(df)
    os.remove(temp_sequences)
    os.remove(temp_output)
    
    #选取符合条件的e-value最小的两行
    seqinfo=[]
    seqlength=[]
    
    for i in duplicate_one_aftercdhit:
            key = list(i.keys())[0]
            seqinfo.append(df.loc[df['sseqid'] == key].nsmallest(1, 'evalue'))
            seqlength.append(df.loc[df['sseqid'] == key].nsmallest(1, 'evalue')['slen'].values[0])
    #如果没有比对到两条序列
    if(len(seqinfo)!=2):
        print("虽然传入两条，但是比对到queryseq上evalue<0.01的比对没有两条")
        return False
    #区分left和right
    if(seqinfo[0]['qstart'].values[0]<seqinfo[1]['qstart'].values[0]):
        left = seqinfo[0]
        right=seqinfo[1]
    else:
        left=seqinfo[1]
        right=seqinfo[0]
    #找出最长的长度
    max_seqlen=max(seqlength)
    print("最长的长度：{}".format(max_seqlen))
    
    #如果左边的序列比对到queryseq上的末尾比右边的序列的开头还要前，可以合并
    if(left['qend'].values[0]<=right['qstart'].values[0] | left['qend'].values[0]-right['qstart'].values[0]<0.25*max_seqlen):
        mergeloc['left']=left['sseqid'].values[0]
        mergeloc['right']=right['sseqid'].values[0]
        return True
    else:
        return False
    # os.remove(temp_queryfile)

def getUnique_one(queryseqid,result,queryProteome,x,i):
    # queryseqid = list(x[i].keys())[0]
    temp_queryseq="/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    cmd="samtools faidx {} '{}' >> {}".format(queryProteome,queryseqid,temp_queryseq)
    subprocess.run(cmd,shell=True)
    duplicates=[]
    uniques=[]
    #将重复的和不重复的分开
    seperate(x[i],duplicates,uniques)
    #按照同一个物种的序列进行分组
    y=duplicateseperate(duplicates)
    # print(y)
    # #获取所有重复物种的id
    duplicatespecieskey = list(y.keys())
    for j in range(len(duplicatespecieskey)):
        currentspecies=duplicatespecieskey[j]
        duplicate_one_cdhit = runCdHit(y[currentspecies])
        # print(duplicate_one_cdhit)
        #将cd-hit结果转化为字典 便于获取序列
        duplicate_one_cdhit_dict={}
        for d in duplicate_one_cdhit:
            seqid=list(d.keys())[0]
            duplicate_one_cdhit_dict[seqid]=d[seqid]
        # print("dict  {}".format(duplicate_one_cdhit_dict))
        if len(duplicate_one_cdhit)==1:
            uniques.append({currentspecies:duplicate_one_cdhit[0]})
            print("cddit处理后只剩一个:{}/{}".format(j,len(duplicatespecieskey)))
        elif len(duplicate_one_cdhit)==2:
            loc={}
            if(ifmerge(temp_queryseq,duplicate_one_cdhit,loc)):
                print("能够合并")
                print("左右位置：{}".format(loc))
                temp={currentspecies:{str(uuid.uuid4()):duplicate_one_cdhit_dict[loc['left']]+duplicate_one_cdhit_dict[loc['right']]}}
                print(temp)
                uniques.append(temp)
                print("重复物种拼接:{}/{}".format(j,len(duplicatespecieskey)))
                # exit(0)
        else:
            print("无法拼接")
            print("重复物种:{}/{}".format(j,len(duplicatespecieskey)))
    #处理好合并以后，将unique中的保存到总字典
    temp={}
    for u in uniques:
        print("u是什么：{}".format(u))
        specieID=list(u.keys())[0]
        print("specieID是什么：{}\n\n".format(specieID))
        seqid=list(u[specieID].keys())[0]
        print("seqid是什么：{}\n\n".format(seqid))
        temp[specieID]=u[specieID][seqid]
    result[queryseqid]=temp
    os.remove(temp_queryseq)
    print("\n\n当前进度:{}/{}\n\n".format(i+1,len(x)))
    
    


parser = argparse.ArgumentParser()
parser.add_argument('--output_1', type=str, help='第一个脚本输出的pickle文件用作当前脚本的输入')
parser.add_argument('--queryProteome', type=str, help='查询的蛋白组 或者 fasta序列文件的文件地址')
parser.add_argument('--output_2', type=str, help='当前脚本输出的去除重复的pickle文件')
parser.add_argument('--num_thread', type=int, help='最多多少个线程一起运行')

args = parser.parse_args()
# filepath="./result/1_test_5seq_200species.pickle"
filepath=args.output_1
x = getAllCandidateOth(filepath)
appendSeqs_All(x,args.num_thread)
result={}   #{"EcoliSeqID1":{"sprciesID1":"seq1","speciesID2":"seq2"},"EcoliSeqID2":{"sprciesID1":"seq1","speciesID2":"seq2"}}
executor = ThreadPoolExecutor(max_workers=args.num_thread)
tasks = []
for i in range(len(x)):
    future = executor.submit(getUnique_one,list(x[i].keys())[0],result,args.queryProteome,x,i)
    tasks.append(future)
# 关闭线程池并等待所有任务完成
executor.shutdown(wait=True)
# 所有任务完成之后执行其他代码
print("去除重复的候选直系同源:finished\n开始保存到文件")
with open(args.output_2,"wb") as f:
    pickle.dump(result,f)
# for i in range(len(x)):
#     queryseqid = list(x[i].keys())[0]
    # temp_queryseq="/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    # cmd="samtools faidx {} '{}' >> {}".format(args.queryProteome,queryseqid,temp_queryseq)
    # subprocess.run(cmd,shell=True)
    # duplicates=[]
    # uniques=[]
    # #将重复的和不重复的分开
    # seperate(x[i],duplicates,uniques)
    # #按照同一个物种的序列进行分组
    # y=duplicateseperate(duplicates)
    # # print(y)
    # # #获取所有重复物种的id
    # duplicatespecieskey = list(y.keys())
    # for j in range(len(duplicatespecieskey)):
    #     currentspecies=duplicatespecieskey[j]
    #     duplicate_one_cdhit = runCdHit(y[currentspecies])
    #     # print(duplicate_one_cdhit)
    #     #将cd-hit结果转化为字典 便于获取序列
    #     duplicate_one_cdhit_dict={}
    #     for d in duplicate_one_cdhit:
    #         seqid=list(d.keys())[0]
    #         duplicate_one_cdhit_dict[seqid]=d[seqid]
    #     # print("dict  {}".format(duplicate_one_cdhit_dict))
    #     if len(duplicate_one_cdhit)==1:
    #         uniques.append({currentspecies:duplicate_one_cdhit[0]})
    #         print("cddit处理后只剩一个:{}/{}".format(j,len(duplicatespecieskey)))
    #     elif len(duplicate_one_cdhit)==2:
    #         loc={}
    #         if(ifmerge(temp_queryseq,duplicate_one_cdhit,loc)):
    #             print("能够合并")
    #             print("左右位置：{}".format(loc))
    #             temp={str(uuid.uuid4()):duplicate_one_cdhit_dict[loc['left']]+duplicate_one_cdhit_dict[loc['right']]}
    #             print(temp)
    #             uniques.append(temp)
    #             print("重复物种:{}/{}".format(j,len(duplicatespecieskey)))
    #             # exit(0)
    #     else:
    #         print("无法拼接")
    #         print("重复物种:{}/{}".format(j,len(duplicatespecieskey)))
    # #处理好合并以后，将unique中的保存到总字典
    # temp={}
    # for u in uniques:
    #     specieID=list(u.keys())[0]
    #     seqid=list(u[specieID].keys())[0]
    #     print("seqid是什么：{}\n\n".format(seqid))
    #     temp[specieID]=u[specieID][seqid]
    # result[queryseqid]=temp
    # os.remove(temp_queryseq)
    # print("\n\n当前进度:{}/{}\n\n".format(i+1,len(x)))
# with open(args.output_2,"wb") as f:
#     pickle.dump(result,f)