import pandas as pd
import keras
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
import helper

######################################################################################
# User input
######################################################################################

# Will want to change the way we do this for the user interface
# Probably have a drop down menu of all genes
temp_gene = 'Adora2a'
temp_compare_gene = 'Amigo2'
file_name = 'Sample_Data.csv'
temp_cell = 'Inhibitory'
temp_compare_cell = "Excitatory"
file_save_location = '' # include?

# Import data - note that I currently am using shorter sample data, so if there is ever
# a null value for a standard error, that's probably because there was only one cell that
# was classified as that cell type
df = pd.read_csv(file_name, sep=',')
clean_df = helper.drop_empty_rows(helper.genes_only(df))

######################################################################################
# Single gene expression data
######################################################################################

# Find highest average gene expression by cell type
def groupBy_avg_expression(df, gene_name):
    avgs = df.groupby("Cell_class").mean()
    sems = df.groupby("Cell_class").sem().rename(columns = {gene_name: gene_name + '_sem'})
    counts = df.groupby("Cell_class").count().rename(columns = {gene_name: gene_name + '_count'})
    grouped_df = pd.merge(avgs, sems, left_index = True, right_index = True, how='inner')
    grouped_df = pd.merge(grouped_df, counts, left_index = True, right_index = True, how='inner').sort_values(by = [gene_name], ascending = False)
    return grouped_df

# Visualize: bar plot 
def single_bar_plot(df, gene_name):
    # Establish variables
    X = df.index
    y = df[gene_name]
    sem = df[gene_name + '_sem']
    labels = df[gene_name + '_count'].tolist()

    # Set up graph
    x_pos = np.arange(len(X))
    fig, ax = plt.subplots()
    ax.bar(x_pos, height = y, yerr = sem, align = 'center', alpha = 0.5, ecolor='black', capsize=10)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(X)
    for i, v in enumerate(y.tolist()):
        ax.text(i, v + 0.01, str(labels[i]), fontweight='bold')
    plt.xlabel('Cell Type')
    plt.ylabel('Average ' + gene_name + ' expression')
    plt.title(gene_name + ' expression by cell type with SEM error bars and number of cells in cell type')
    # plt.show()
    plt.savefig(gene_name + '_expression.png')

######################################################################################
# Single gene expression data DEMO
######################################################################################

grouped = groupBy_avg_expression(clean_df, temp_gene)
# print(grouped[[temp_gene, temp_gene + '_sem', temp_gene + '_count']])
single_bar_plot(grouped, temp_gene)
grouped = groupBy_avg_expression(clean_df, temp_compare_gene)
# print(grouped[[temp_compare_gene, temp_compare_gene + '_sem', temp_compare_gene + '_count']])
single_bar_plot(grouped, temp_compare_gene)

######################################################################################
# Multiple genes expression data
######################################################################################

def groupBy_avg_expression_correl(df, gene_names):
    if type(gene_names) != list or len(gene_names) != 2:
        print("Second argument must be a list of two genes.")
        return None
    # Get the cell types
    cell_types = df.Cell_class.unique()
    # Set up output
    correls = {}
    for cell_type in cell_types:
        # Filter out the undesired cell types
        grouped = df[df["Cell_class"] == cell_type]
        # Get the correlation between the two genes
        corr = grouped[gene_names[0]].corr(grouped[gene_names[1]])
        correls[cell_type] = corr
    out = pd.DataFrame.from_dict(correls, orient = 'index').dropna().rename(columns = {0: 'correl'})
    return out.sort_values(by = 'correl', ascending = False)

# Visualize: bar plot with 
def multi_bar_plot(df, gene_names):
    # Establish variables
    X = df.index
    y = df['correl']

    # Set up graph
    x_pos = np.arange(len(X))
    fig, ax = plt.subplots()
    ax.bar(x_pos, height = y, align = 'center', alpha = 0.5, ecolor='black', capsize=10)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(X)
    plt.xlabel('Cell Type')
    plt.ylabel('Correlation')
    plt.title('Correlation between ' + gene_names[0] + ' and ' + gene_names[1] + ' by cell type')
    plt.show()
    # plt.savefig(gene_names[0] + '_' + gene_names[1] + '_correlation.png')

######################################################################################
# Multiple gene expression data DEMO
######################################################################################

temp_genes = [temp_gene, temp_compare_gene]
grouped = groupBy_avg_expression_correl(clean_df, temp_genes)
multi_bar_plot(grouped, temp_genes)









# ######################################################################################
# # Find genes with the highest correlation in cell types
# ######################################################################################

# def highest_gene_correl(df, cell_types):
#     if type(cell_types) != list or len(cell_types) != 2:
#         print("Second argument must be a list of two cell types.")
#         return None
#     out = {}
#     for cell_type in cell_types:
#         # Filter out the undesired cell types
#         grouped = df[df["Cell_class"] == cell_type]
#         out[cell_type] = grouped.drop("Cell_class", axis = 1)
#     corr = helper.corr_cells(out[cell_types[0]], out[cell_types[1]])
#     corr2 = out[cell_types[0]].corrwith(out[cell_types[1]])
#     # Show figure
#     print(corr)
#     plt.matshow(corr.astype(int))
#     # Plot
#     plt.xticks(range(len(corr.columns)), corr.columns, rotation=90)
#     plt.yticks(range(len(corr.columns)), corr.columns)
#     plt.colorbar()
#     plt.rcParams.update({'font.size':10})
#     plt.gcf().set_size_inches(14.5, 14.5)
#     plt.show()
#     # plt.title('Correlation of All Strains', pad = 60)
#     # plt.savefig('feature_exploration/corr.png')
#     # plt.close()
#     return corr

# ######################################################################################
# # Find genes with the highest correlation in cell types DEMO
# ######################################################################################
# cell_types = [temp_cell, temp_compare_cell]
# out = highest_gene_correl(clean_df, cell_types)
