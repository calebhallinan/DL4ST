

import scanpy as sc

# Import libraries
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from PIL import Image
import matplotlib.patches as patches




# NOTE: needed to change "spacial/tissue_positions.csv" --> "spacial/tissue_positions_list.csv"

# file name
file_name = "breastcancer_xenium_sample1_rep1"
# resolution
resolution = 12
# read in the data
adata_xenium = sc.read_10x_h5('/Users/calebhallinan/Desktop/jhu/fanlab/projects/histology_to_gene_prediction/data/breastcancer_xenium_sample1/outs/cell_feature_matrix.h5')

# Load the full-resolution spatial data
cell_centers = pd.read_csv(f"/Users/calebhallinan/Desktop/jhu/fanlab/projects/histology_to_gene_prediction/data/{file_name}/{file_name}_fullresolution_STalign.csv.gz", index_col=0)
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/xenium_cell_centroids_visium_high_res.csv")
# cell_centers = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_visium/scaled_spots_for_xenium_image.csv", index_col=0)
# cell_centers.columns = ["x_centroid", "y_centroid"]

# Load the full-resolution image
Image.MAX_IMAGE_PIXELS = None
img_name = "Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image"
img = np.array(Image.open("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/" + file_name + "/" + img_name + ".tif"))
# img = np.load("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/janesick_nature_comms_2023_companion/visium_high_res_image.npy")
# plt.imshow(img)

# add .obs
adata_xenium.obs = cell_centers
# add .obsm
adata_xenium.obsm["spatial"] = adata_xenium.obs[["x_centroid", "y_centroid"]].to_numpy().astype(int)
# add image
adata_xenium.uns['spatial'] = img
# need to add this for subsetting
adata_xenium.obs.index = adata_xenium.obs.index.astype(str)

# adata_xenium.X = np.arcsinh(adata_xenium.X).toarray()

# scale genes with cpm
# sc.pp.normalize_total(adata_xenium, target_sum=1e6)

# log transform the data
# sc.pp.log1p(adata_xenium)

# sc.pp.normalize_total(adata_xenium, target_sum=1e6)

# get rid of genes that aren't in visium
# gene_list = pd.read_csv("/home/caleb/Desktop/projects_caleb/histology_to_gene_prediction/data/breastcancer_xenium_sample1_rep1/rastGexp_df.csv", index_col=0)
# gene_list = [gene for gene in gene_list.index if "BLANK" not in gene and "Neg" not in gene and  "antisense" not in gene]
# gene_list = [gene for gene in gene_list if gene not in ['AKR1C1', 'ANGPT2', 'APOBEC3B', 'BTNL9', 'CD8B', 'POLR2J3', 'TPSAB1']]
# # subset the data
# adata_xenium = adata_xenium[:, gene_list]

# make an array of the gene expression data
adata_xenium.X_array = pd.DataFrame(adata_xenium.X.toarray(), index=adata_xenium.obs.index)

# plot
plt.scatter(adata_xenium.obsm["spatial"][:,0], adata_xenium.obsm["spatial"][:,1], s=.01, c="black")
plt.axis("off")
# plt.savefig("xenium1_spots.svg", dpi=300)
