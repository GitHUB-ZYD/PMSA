# PMSA的构建

使用RoseTTAFold进行蛋白互作预测需要构建蛋白序列的PMSA文件。该项目编写了一套生成PMSA(Paired Multiple Sequences)的脚本，将PMSA的构建流程化。**仅适用于细菌的蛋白序列**

# 使用流程

1. ### 下载细菌蛋白组

   allProteomeUrl.csv文件中保存有Assembly数据库中标记为representative的细菌蛋白组的下载路径；

2. ### 使用PMSA.yaml创建脚本运行环境

   ```
   conda env create -f PMSA.yaml
   ```

3. ### 使用makeblastdb为每一个细菌蛋白组建立比对库

4. ### 使用samtools faidx为每一个细菌蛋白组建立fai索引文件

5. ### 生成PMSA文件，使用RF2t模型进行蛋白互作预测

   RF2tExample目录中为脚本使用示例。在该目录中，queryprot/test.fasta保存着两条蛋白序列，已经使用makeblastdb对这两条序列建库，并且使用samtools faidx为其建立索引，该脚本将为这两条序列构建PMSA文件；scripts目录保存脚本源码；ref.loc保存下载到本地的细菌蛋白组的文件全路径；results目录保存脚本所有的输出文件

   1. ```python
      #输出每一条蛋白序列的候选直系同源组
      time python scripts/1_getCandidateOrthologues.py --queryProteome queryprot/test.fasta --species_loc ./ref.loc --output_1 result/output_1.pickle --num_thread 4
      
      #使用viewOutput_1.py可以产看output_1.pickle中包含的信息
      python scripts/viewOutput_1.py --output_1 result/output_1.pickle --savetxt result/output_1.txt
      ```

   2. ```python
      #对于一个直系同源蛋白组，该组中每一个物种只能有一条蛋白序列，这一步就是对output_1.pickle中重复的物种的序列进行处理
      time python scripts/2_removeMultipleCandidates.py --output_1 result/output_1.pickle --queryProteome queryprot/test.fasta --output_2 result/output_2.pickle --num_thread 50
      
      #使用viewOutput_2.py可以产看output_2.pickle中包含的信息，该文件保证每一个同源组中每一个物种只有一条蛋白，并且还保存了序列的信息
      python scripts/viewOutput_2.py --output_2 result/output_2.pickle --savetxt result/output_2.txt
      
      ```

   3. ```python
      #使用mafft分别对每一个直系同源蛋白组进行局部比对，并使用trimal处理gap
      time python scripts/3_getAlignment.py --output_2 result/output_2.pickle  --queryProteome queryprot/test.fasta --output_3 result/output_3.pickle  --num_thread 2
      
      #使用viewOutput_3.py可以产看output_3.pickle中包含的信息，该pickle文件保存有每一个直系同源组的比对后的蛋白序列，第一条是query，后面的蛋白用蛋白组的路径来命名
      python scripts/viewOutput_3.py --output_3 result/output_3.pickle --savetxt result/output_3.txt
      
      ```

   4. ```python
      #统计每一个直系同源组的各种信息，包括序列长度，并且将每一个直系同源组两两配对，输出具有的相同物种的数目，输出为csv
      python scripts/4_gatherOrthoStatic.py --output_3 result/output_3.pickle --output_4_csv result/output_4.csv
      
      ```

   5. ```python
      #当从output_4.csv获得要配对的蛋白id以后，使用该脚本进行配对，生成pmsa。输入seq1的ID和seq2的ID，就能够将相同物种的两条序列合并
      python scripts/5_getPMSA.py --output_3 result/output_3.pickle --seq1 seq_218 --seq2 seq_310 --output_5_aln result/seq1_seq2.aln
      
      ```

      完成以上步骤以后，生成的PMSA文件(seq1_seq2.aln)就可以传入RoseTTAFold进行互作预测

   6. ```python
      #如果蛋白序列的数目较多，一个个的构建pmsa来计算互作概率显得很麻烦，因此可以对output_4.csv进行过滤，选出那些想要用于预测的蛋白对，然后传入6_getPPIProbability.py中批量计算互作概率
      
      #运行这一步需要安装RoseTTAFold的环境依赖，并且下载好RF2t所需要的权重文件。
      ## https://github.com/RosettaCommons/RoseTTAFold
      python scripts/6_getPPIProbability.py  --input_csv output_4.csv --RTF2track RoseTTAFold脚本的路径/network_2track/predict_msa.py --output_3 result/output_3.pickle --save_ppi result/ppiprobability.csv --cpu 0
      
      ```

      **通过-h参数能够了解到每一个参数的具体含义**

# 注意事项

1. 多序列比对的比对质量影响PMSA的质量，PMSA的质量对RF2t的准确性有很大的影响。
2. 不能直接构建长蛋白序列(>300AA)的PMSA文件。当蛋白序列较长，序列数目较大的时候，mafft多序列比对质量较差，导致生成的PMSA的质量较差，可以在 3_getAlignment.py 给调用mafft的命令指定最大迭代次数，如 --maxiterate 1000，提高长序列的多序列比对质量，但是相应会需要过多的计算资源。在测试过程中，指定最大迭代次数为500，在比对过程直接吃完64g内存和32g SWAP，导致死机无法得到测试结果，因此脚本没有给出设置最大迭代次数的参数，而选择让mafft自己选择合适的迭代次数。

# 参考文献

Baek M, Dimaio F, Anishchenko I, et al. Accurate prediction of protein structures and interactions using a three-track neural network[J]. Science, 2021,373(6557):871-876.

Humphreys I, Pei J, Baek M, et al. Computed structures of core eukaryotic protein complexes[J]. Science, 2021,374(6573):m4805.



