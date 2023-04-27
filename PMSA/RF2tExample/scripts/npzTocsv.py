#将 rosetta fold 生成的npz文件转化为csv表格
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--output_npz', type=str, help='RoseTTAFold输出的npz文件作为输入')
parser.add_argument('--output_6_csv', type=str, help='两条蛋白质各个残基的距离矩阵')
args = parser.parse_args()

data = np.load(args.output_npz)
df = pd.DataFrame(data['dist'])
df.to_csv(args.output_6_csv,index=False)