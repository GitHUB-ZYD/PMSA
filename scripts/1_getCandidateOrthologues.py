import subprocess
import pandas as pd
import pickle
import threading
import uuid
import os
import argparse

#获取query Proteome所有的序列id
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

#queryFilePath: query序列的来自哪一个蛋白组文件
#querySeqID: 要查询的query seq的序列ID
#dbpath: query seq和哪一个蛋白组比对
def getBestHits(queryFilePath="UP000000625.fasta",querySeqID="sp|P00350|6PGD_ECOLI",dbpath="GCF_000308975.1_ASM30897v2_protein.faa"):
    temp_file = "/dev/shm/tmp/"+str(uuid.uuid4()) + ".csv"
    # print(temp_file)
    #运行blastp
    temp_fasta = "/dev/shm/tmp/"+str(uuid.uuid4())+".fasta"
    cmd = "samtools faidx {} '{}' > {}".format(queryFilePath,querySeqID,temp_fasta)
    subprocess.run(cmd,shell=True)
    # 读取blastp结果到Pandas DataFrame
    cmd = "blastp -query {} -db {} -outfmt '6 delim=, qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident ppos positive' -out {} -evalue 0.01 -num_threads=14".format(temp_fasta,dbpath,temp_file)
    subprocess.run(cmd, shell=True)
    columns = ["qseqid","qlen","sseqid","slen","qstart","qend","sstart","send","evalue","length","pident","nident","ppos","positive"]
    df = pd.read_csv(temp_file, sep=',', header=None, names=columns)
    #计算query_coverage和hit_coverage并排序按照identity和evalue
    df['query_coverage'] = df.apply(lambda row:row['length']/row['qlen'],axis=1)
    df['hit_coverage'] = df.apply(lambda row:row['length']/row['slen'],axis=1)
    # print("\n\ndf的长度:",len(df))
    # print(df)
    filtered_df = df.loc[((df['slen']>=50) | (df['query_coverage']>=0.5) ) & ((df['query_coverage']>=0.5) | (df['hit_coverage']>=0.5))]
    # print("filtered_df的长度",len(filtered_df)) 

    sorted_identity =  filtered_df.sort_values(by='pident', ascending=False, inplace=False)

    sorted_evalue = filtered_df.sort_values(by='evalue',ascending=True,inplace=False)
    
    #identity
    bestidentity = sorted_identity['pident'].max()
    # print("bestidentity:{}".format(bestidentity))
    sorted_identity_filter = sorted_identity.loc[((sorted_identity['positive']>=0.25*sorted_identity['length']) & (sorted_identity['pident']>bestidentity-10)) | ((sorted_identity['positive']<0.25*sorted_identity['length']) & (sorted_identity['pident']>bestidentity-20))]

    #evalue
    bestevalue = sorted_evalue['evalue'].min()
    #如果是0,那就10的-10次方，不然无论如何能通过
    if(bestevalue==0):
        bestevalue=1e-10
    # print("bestevale",bestevalue)
    sorted_evalue_filter = sorted_evalue.loc[((sorted_evalue['positive']>=0.25*sorted_evalue['length'])&(sorted_evalue['evalue']<bestevalue * 1e5))|((sorted_evalue['positive']<0.25*sorted_evalue['length']))]
    # print("evalue过滤长度:",len(list(sorted_evalue_filter['sseqid'])))
    # print("identity过滤长度",len(list(sorted_identity_filter['sseqid'])))

    #取交集
    intersection = list(sorted_evalue_filter[sorted_evalue_filter['sseqid'].isin(sorted_identity_filter['sseqid'])]['sseqid'])

    # print(df)
    # df.to_csv('output.csv', index=False)
    os.remove(temp_file)
    os.remove(temp_fasta)
    # print("交集长度:",len(intersection))
    return intersection
    #返回列表

# #输出线程名字
def getOneSeqOth(threadname,queryFilePath="UP000000625.fasta",querySeqID="sp|P00350|6PGD_ECOLI",proteome_species=['GCF_000006685.1_ASM668v1_protein.faa','GCF_000308975.1_ASM30897v2_protein.faa']):
    #获得这条querySeq的直系同源组
    candiadateOrtho_one=[]
    i=0
    for proteome_specie in proteome_species:
        
        #一个物种的foreard best hits
        fbhs = getBestHits(queryFilePath=queryFilePath,querySeqID=querySeqID,dbpath=proteome_specie)
        # print(fbhs,"fbh success")
        for fbh in fbhs:
            #一个物种的forward best hits中的一条forward best hit的reverse best hits
            rbhs = getBestHits(queryFilePath=proteome_specie,querySeqID=fbh,dbpath=queryFilePath)
            # print(rbhs,"rbh success")
            #如果，reverseBestHits中有这条querySeqID,说明当前这条fbh属于这条querySeqID的直系同源蛋白序列
            if querySeqID in rbhs:
                candiadateOrtho_one.append({proteome_specie:fbh})
        i=i+1
        print("{}当前species:{}".format(threadname,i))
    #if(len(candiadateOrtho_one)==0):
    #    return None
    #else:
    #    return candiadateOrtho_one
    return candiadateOrtho_one

def getOneSeqOth_mul(num_threads=14,queryFilePath="UP000000625.fasta",querySeqID="sp|P00350|6PGD_ECOLI",proteome_species=['GCF_000006685.1_ASM668v1_protein.faa','GCF_000308975.1_ASM30897v2_protein.faa']):
    # 将数据切分成 num_threads 份
    chunk_size = len(proteome_species) // num_threads
    data_chunks = [proteome_species[i:i + chunk_size] for i in range(0, len(proteome_species), chunk_size)]
    # for chunk in data_chunks:
    #     print(chunk,sep='\n')
    # 创建线程来处理数据
    threads = []
    results = []
    i=0
    for chunk in data_chunks:
        i=i+1
        threadname = "getOneSeqOth_Thread"+str(i)
        thread = threading.Thread(target=lambda threadname,queryFilePath,querySeqID,proteome_species: results.append(getOneSeqOth(threadname,queryFilePath,querySeqID,proteome_species)), args=(threadname,queryFilePath,querySeqID,chunk))
        thread.start()
        threads.append(thread)

    # 等待所有线程完成
    for thread in threads:
        thread.join()

    # 合并结果
    return [item for sublist in results for item in sublist]


#多线程实现
#该函数用于处理已经被线程分割后的列表querySeqs
def process_querySeq(querySeqs,threadname,queryProteome="UP000000625.fasta",species=['GCF_000006685.1_ASM668v1_protein.faa', 'GCF_000308975.1_ASM30897v2_protein.faa']):
    candiadateOrtho_all=[]
    i=0
    for querySeq in querySeqs:
        #这条Ecoli蛋白的直系同源组
        candidateOrtho_one =getOneSeqOth_mul(200,queryFilePath=queryProteome,querySeqID=querySeq,proteome_species=species)
        if(candidateOrtho_one!=None):
            candiadateOrtho_all.append({querySeq:candidateOrtho_one})  
        i=i+1
        print("{}: {}/{}".format(threadname,i,len(querySeqs)))
        # print("candiadateOrtho_all的长度",len(candiadateOrtho_all))
        #  candiadateOrtho_all的数据结构：
            # [
            #   {
            #     Ecoli蛋白id:
            #     [candidateOrtho_one,{物种id,fbh序列},...]
            #   },
            #   {
            #     Ecoli蛋白id:
            #     [{物种id,fbh序列},{物种id,fbh序列},...]
            #   }
            # ]
        #
        #用于测试，只能每个线程3个蛋白
        # if i==3:
        #     print(threadname,"finish success")
        #     return candiadateOrtho_all
    print(threadname,"finish success")
    return candiadateOrtho_all


def process_querySeq_multithreaded(querySeqs, num_threads,queryProteome,species):
    # 将数据切分成 num_threads 份
    chunk_size = len(querySeqs) // num_threads
    data_chunks = [querySeqs[i:i + chunk_size] for i in range(0, len(querySeqs), chunk_size)]
    # for chunk in data_chunks:
    #     print(chunk,sep='\n')
    # 创建线程来处理数据
    threads = []
    results = []
    i = 0
    for chunk in data_chunks:
        i=i+1
        threadname="thread"+str(i)
        thread = threading.Thread(target=lambda querySeqs,threadname,queryProteome,species: results.append(process_querySeq(querySeqs,threadname,queryProteome,species)), args=(chunk,threadname,queryProteome,species))
        thread.start()
        threads.append(thread)

    # 等待所有线程完成
    for thread in threads:
        thread.join()

    # 合并结果
    return [item for sublist in results for item in sublist]


if __name__=="__main__":
    
    # 要查询的蛋白组文件 或者 从蛋白组中抽出的fasta序列文件
    # queryProteome="UP000000625.fasta"
    parser = argparse.ArgumentParser()
    parser.add_argument('--queryProteome', type=str, help='查询的蛋白组 或者 fasta序列文件的文件路径')
    parser.add_argument('--species_loc', type=str, help='所有用于提取直系同源蛋白的蛋白组路径')
    parser.add_argument('--output_1', type=str, help='输出的候选直系同源蛋白的pickle文件')
    parser.add_argument('--num_thread', type=int, help='运行线程数目: 1-4')
    args = parser.parse_args()
    if not os.path.exists('/dev/shm/tmp'):
    # 如果目录不存在，创建目录
        os.makedirs('/dev/shm/tmp')
    if not os.path.exists('./result'):
    # 如果目录不存在，创建目录
        os.makedirs('result')
    # queryProteome="./queryprot/test.fasta"
    queryProteome=args.queryProteome
    # 该蛋白组 或者 fasta序列文件的所有序列id
    querySeqs = getEcoliSeqID(EcoliFilePath=queryProteome)
    # 所有用于提取直系同源蛋白的序列
    species=[]
    # with open("./species_loc_200.txt","r") as f:
    #     x=f.readlines()
    with open(args.species_loc,"r") as f:
        x=f.readlines()
    for i in x:
        species.append(i.rstrip())
    print(len(species))
    allOth = process_querySeq_multithreaded(querySeqs=querySeqs,num_threads=args.num_thread,queryProteome=queryProteome,species=species)
    # with open('./result/1_test_5seq_200species.pickle',"wb") as f:
    #     pickle.dump(allOth,f)
    with open(args.output_1,"wb") as f:
        pickle.dump(allOth,f)
    