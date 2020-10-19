import math
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import sys

def sigmoid(x):
    return 1 / (1 + math.exp(-x))

# Filter out irrelevant columns
def genes_only(df):
    return df.drop(['Cell_ID','Animal_ID','Animal_sex','Behavior','Bregma','Centroid_X','Centroid_Y',
        'Neuron_cluster_ID'], axis = 1)

# Drop null rows
def drop_empty_rows(df):
    return df.dropna(axis = 1)

# # Calculate correlation coefficient from dataframes
# def correlation_coeff(col1, col2):
#     if type(col1) != list or type(col2) != list:
#         print("Columns must be type list to calculate correlation coefficient.")
#         return None
#     x_bar = np.mean(col1)
#     y_bar = np.mean(col2)
#     num = 0
#     denom_x = 0
#     denom_y = 0
#     for i in col1:
#         for j in col2:
#             num += (i - x_bar) * (j - y_bar)
#             denom_x += (i - x_bar) ** 2
#             denom_y += (j - y_bar) ** 2
#     print(num / (math.sqrt(denom_x * denom_y)))
#     print(pearsonr(col1, col2))
#     sys.exit()

# def corr_cells(df1, df2):
#     out = pd.DataFrame(index = df1.columns, columns = df2.columns)
#     for row in df1.columns:
#         for col in df2.columns:
#             out[col][row] = round(correlation_coeff(df1[row].tolist(), df2[col].tolist()), 4)
#     return out