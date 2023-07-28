import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import KernelPCA

# Frobenius norm function
def frob(real_df, reconstructed_df):
  sub = np.subtract(real_df, reconstructed_df)
  sq = np.square(sub)
  summed = np.sum(sq)
  final = np.sqrt(summed)
  return final

# Adding labels onto 
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],y[i].round(3), va = 'top')

def sorting(pair):
    return abs(pair[1])

x = pd.read_csv('pathways.csv')
data_set = pd.read_csv('pathways_names.csv').T
x = x.to_numpy()
# x = (x - np.mean(x, axis=0)) / np.std(x, axis=0)
data_set = data_set.to_numpy()

# Testing different types of kernels
kernels = ['rbf', 'cosine', 'linear', 'poly']
for ker in kernels:
    exvarsum = 0
    pc_count = 0
    kpca = KernelPCA(kernel=ker, fit_inverse_transform=True)
    # x_kpca = kpca.fit_transform(x)
    # kpca_transform = kpca.fit_transform(x)
    where = kpca.fit_transform(x.T)
    what = kpca.eigenvectors_
    help = kpca.eigenvalues_
    huh = kpca.inverse_transform(where)

    ev = help / np.sum(help)

    pc_c = 0
    sum = 0
    while sum < 0.8:
        vari = help[pc_c] / np.sum(help)
        sum = sum + vari
        pc_c = pc_c + 1
        if sum >= 0.8 or pc_c == 60:
            break
        

    print(ker + ": ", sum.round(4), "number of pc: ", pc_c)

    # explained_variance = np.var(kpca_transform, axis=0)
    # ev = explained_variance / np.sum(explained_variance)
    # while exvarsum < 0.8:
    #     exvarsum = exvarsum + ev[pc_count]
    #     pc_count = pc_count + 1

    # print(ker + " explained variance: ", exvarsum.round(4), "number of pc: ", pc_count)

    #--------- Bar Graph for Explained Variance Ratio ------------
    plt.bar(list(range(1, 16)),list(ev[:15]),label='Principal Components',color="#0072BD")
    plt.legend()
    plt.xlabel('Principal Components ')
    #----------------------
    n=list(ev[:15])
    pc=[]
    for i in range(len(n)):
            n[i]=round(n[i],4)
            pc.append('PC-'+str(i+1)+'('+str(n[i])+')')

    #----------------------
    plt.xticks(list(range(1, 16)),pc, fontsize=7, rotation=30)
    plt.ylabel('Explained Variance')
    plt.title('Explained Variance Using Kernel:'+str(ker))
    addlabels(list(range(1, 16)),list(ev[:15]))
    plt.savefig(ker + '.csv' + '.png')
    plt.clf()

# hehes = ['rbf', 'poly']
# gammas = [1/10000, 1/1000, 1/100]
# for hehe in hehes:
#     for gamma in gammas:
#         exvarsum = 0
#         pc_count = 0
#         kpca = KernelPCA(kernel=hehe, fit_inverse_transform=True, gamma= gamma)
#         # x_kpca = kpca.fit_transform(x)
#         # kpca_transform = kpca.fit_transform(x)
#         where = kpca.fit_transform(x.T)
#         what = kpca.eigenvectors_
#         help = kpca.eigenvalues_
#         huh = kpca.inverse_transform(where)

#         ev = help / np.sum(help)

#         pc_c = 0
#         sum = 0
#         while sum < 0.8:
#             vari = help[pc_c] / np.sum(help)
#             sum = sum + vari
#             pc_c = pc_c + 1
#             if sum >= 0.8 or pc_c == 60:
#                 break
            

#         print(hehe + ": ", sum.round(4), "number of pc: ", pc_c, "gamma: ", gamma)


# fkpca = KernelPCA(kernel='linear', gamma=15, fit_inverse_transform=True)
# fx_kpca = fkpca.fit_transform(x.T)

# hello = fkpca.eigenvectors_
# what = fkpca.eigenvalues_

# pc_c = 1
# sum = 0
# while sum < 0.8:
#     sum = np.sum(np.square(hello[:,:pc_c])) / np.sum(np.square(hello))
#     if sum >= 0.8:
#         break
#     pc_c = pc_c + 1

# print(pc_c)

# setx = set()
# repeats = set()

# assoc_genes = []
# for i in range(0, 8):
#     gn = 0
#     each_pc = []
#     gene_numbers = []
#     for col in x.T:
#         correlation = np.corrcoef(fx_kpca[...,i], col)
#         gene = (gn, correlation[0][1])
#         # print(correlation[0][1])
#         gene_numbers.append(gene)
#         gn = gn + 1
#     gene_numbers.sort(key = sorting, reverse=True)
#     for j in range(0, 10):
#         each_pc.append(gene_numbers[j])
#         # print(gene_numbers[j])

#     assoc_genes.append(each_pc)

# for i in assoc_genes[2]:
#     print(data_set[1][i[0]])

# from scipy.cluster.hierarchy import dendrogram, linkage

fkpca = KernelPCA(kernel='linear', fit_inverse_transform=True)
fx_kpca = fkpca.fit_transform(x.T)
pcs = fkpca.eigenvectors_



np.random.seed(0)
from sklearn.cluster import AgglomerativeClustering
# from collections import defaultdict

fig, axs = plt.subplots(nrows=7, ncols=7)
plt.suptitle('Hierarchical Clustering of Gene Pathways on PCs 1-7')
fig.supxlabel('PC 2') 
fig.supylabel('PC 1')

print(axs.shape)
axs
gene_clusters = [[] for _ in range(1283)]
first_gene_cluster = []
second_gene_cluster = []
title = 0

for i in range(0, 7):
    for j in range(0, 7):
        group = np.empty((1283, 2))
        for k in range(0, 1283):
            # x.append(pcs[i, 3])
            # y.append(pcs[i, 4])
            vector = np.array([pcs[k, i], pcs[k, j]])
            group[k] = vector       

            # Perform hierarchical clustering
            

        ax = axs[i][j]
        clustering = AgglomerativeClustering(n_clusters=10, linkage='ward')
        clustering.fit(group)

        # Get cluster labels
        labels = clustering.labels_
        for p in range(0, 1283):
            gene_clusters[p].append(labels[p])

        ax.scatter(group[:, 0], group[:, 1], c=labels, cmap='inferno')
        ax.label_outer()
        
        plt.show()


        # max_gene = ''
        # max_inter = 0
        # for m in range(1, 1282):
        #     for n in range(m + 1, 1283):
        #         if labels[m] == labels[n]:
        #             max_inter = max_inter + 1
            # second_gene_cluster.append(labels[p])
            # intersection = set(first_gene_cluster) & set(second_gene_cluster)
            # if len(intersection) > max_inter:
            #     max_inter = len(intersection)
            #     max_gene = data_set[1][p]

        # print(max_gene, max_inter)



        # dict = defaultdict(list)
        # # dict = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[]}
        # k = 0
        # for name in data_set[1]:
        #     key = labels[k]
        #     dict[key].append(name)
        #     k = k + 1

        # list_8_8_first_cluster = dict[0]
        # print(list_8_8_first_cluster)

        # Plot the data points with different colors for each cluster

def intersection(list1, list2):
    total = 0
    for i in range(0, len(list1)):
        if list1[i] == list2[i]:
            total = total + 1
    return total

comparisons = np.zeros((1283, 1283))
for i in range(0, 1283):
    for j in range(0, 1283):
        inter = intersection(gene_clusters[i], gene_clusters[j])
        prop = round(inter / 64, 3)
        comparisons[i][j] = prop

row_num = 0
most_corr = [[] for _ in range(1283)]
most_corr_num = [[] for _ in range(1283)]
for row in comparisons:
    most_corr[row_num].append(data_set[1][row_num])
    most_corr_num[row_num].append(row_num)
    for i in range(0, len(row)):
        if row[i] > 0.6:
            most_corr[row_num].append(data_set[1][i])
            most_corr_num[row_num].append(row[i])

    row_num = row_num + 1

mcdf = pd.DataFrame(most_corr)
mcndf = pd.DataFrame(most_corr_num)
mcdf.fillna('', inplace=True)
mcdf.to_csv('most_corr_cosine.csv')

mcndf.fillna('', inplace=True)
mcndf.to_csv('most_corr_num_cosine.csv')


# np.savetxt('most_corr.csv', most_corr)
# np.savetxt('most_corr_num.csv', most_corr_num)

print(most_corr[0])







# Perform hierarchical clustering
# Z = linkage(group, method='ward')

# # Plot the dendrogram
# plt.figure(figsize=(8, 6))
# dendrogram(Z)
# plt.xlabel('Data Point')
# plt.ylabel('Distance')
# plt.title('Hierarchical Clustering Dendrogram')
# plt.show()


# plt.scatter(x, y)

# # Set the x and y axis labels
# plt.xlabel('X')
# plt.ylabel('Y')

# # Set the title of the plot
# plt.title('Scatter Plot')

# # Show the plot
# plt.show()

# fig, ax = plt.subplots()
# # colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'pink', 'black', 'grey', 'magenta']
# labels = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15']

# # Plot the vectors
# for vector, label in zip(group,labels):
#     ax.quiver(0, 0, vector[0], vector[1], angles='xy', scale_units='xy', scale=1, color = 'blue')
#     ax.text(vector[0], vector[1], label, ha='left', va='bottom', color='black')

# # Set the x and y axis limits
# # ax.set_xlim([-.1, .1])
# # ax.set_ylim([-.1, .1])

# # Add labels and grid
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.grid(True)

# # Show the plot
# plt.show()

reconstruct = fkpca.inverse_transform(fx_kpca)
# print(x[0][0], reconstruct[0][0])
# print(x[0][1], reconstruct[0][1])
# print(x[0][2], reconstruct[0][2])
reconstruction_error = np.mean(np.square(x.T - reconstruct))
print(reconstruction_error)

third = KernelPCA(kernel='poly', fit_inverse_transform=True)
thirdm = third.fit_transform(x)
reconstruct_third = third.inverse_transform(thirdm)
print("poly err:", np.mean(np.square(x - reconstruct_third)))

# print("poly:", frob(x, reconstruct_third))

a = KernelPCA(kernel='rbf', fit_inverse_transform=True, gamma = 1/1283)
am = a.fit_transform(x)
reconstruct_a = a.inverse_transform(am)
print("rbf err:", np.mean(np.square(x - reconstruct_a)))

# print("rbf:", frob(x, reconstruct_a))

b = KernelPCA(kernel='cosine', fit_inverse_transform=True, gamma = 1/1283)
bm = b.fit_transform(x)
reconstruct_b = b.inverse_transform(bm)
print("cosine err:", np.mean(np.square(x - reconstruct_b)))

# print("cosine:", frob(x, reconstruct_b))

# import pandas as pd
# import plotly.express as px

# n_components = 4

# labels = {str(i): f"PC {i+1}" for i in range(n_components)}
# labels['color'] = 'Median Price'

# fig = px.scatter_matrix(
#     fx_kpca,
#     dimensions=range(n_components),
#     labels=labels,
# )
# fig.update_traces(diagonal_visible=False)
# fig.show()


