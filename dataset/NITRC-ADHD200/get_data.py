from nilearn import datasets, plotting, image
from nilearn.connectome import ConnectivityMeasure
from nilearn.input_data import NiftiLabelsMasker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import csv
import os
import pickle

# Download the ADHD200 data
adhd_dataset = datasets.fetch_adhd(n_subjects=30)
phenotypic_data = adhd_dataset.phenotypic
print(pd.DataFrame(phenotypic_data).head())

# Load the Harvard-Oxford atlas
n_rois = 100
schaefer_atlas = datasets.fetch_atlas_schaefer_2018(
    n_rois=n_rois
)  # up to 1000 (ROI, number nodes) for this atlas
atlas_filename = schaefer_atlas.maps
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)

connectivity_measure = ConnectivityMeasure(kind="correlation")

node_coords_path = f"dataset/NITRC-ADHD200/brains/node_coords_{n_rois}.pkl"
if os.path.exists(node_coords_path):
    with open(node_coords_path, "rb") as pkl_file:
        node_coords = pickle.load(pkl_file)
else:
    node_coords = plotting.find_parcellation_cut_coords(labels_img=atlas_filename)
    os.makedirs(os.path.dirname(node_coords_path), exist_ok=True)
    with open(node_coords_path, "wb") as pkl_file:
        pickle.dump(node_coords, pkl_file)


def create_connectivity_matrix(f):
    m = connectivity_measure.fit_transform([masker.fit_transform(f)])[0]
    np.fill_diagonal(m, 0.0)  # no self-loops
    return m


cm = [create_connectivity_matrix(f) for f in np.array(adhd_dataset.func)]

# Building the adjacency matrices
Gs = []
adj_matrices = []
edge_cutoff_threshold = 0.5
for connectivity_matrix in cm:
    adj_matrix = (connectivity_matrix >= edge_cutoff_threshold).astype(int)
    adj_matrices.append(adj_matrix)
    G = nx.from_numpy_array(adj_matrix)

    Gs.append(G)
    print(f"number of nodes: {G.number_of_nodes()}, number of edges: {G.number_of_edges()}")

# Save the adjacency matrices as edge lists in a csv file
for i, G in enumerate(Gs):
    subject_id = phenotypic_data["Subject"][i]
    G.graph["subject_id"] = subject_id

    edge_list = list(nx.generate_edgelist(G, data=False))
    edge_list = [edge.split(" ") for edge in edge_list]
    edge_list_df = pd.DataFrame(edge_list, columns=["source", "target"])
    edge_list_df["subject_id"] = subject_id

    edge_list_path = f"dataset/NITRC-ADHD200/brains/edge_list_{subject_id}.csv"

    edge_list_df.to_csv(header=True, path_or_buf=edge_list_path, index=False)
    print(f"Edge list saved to {edge_list_path}")

# Saving phenotypic data as csv
phenotypic_data_df = pd.DataFrame(phenotypic_data)
phenotypic_data_df.to_csv(
    header=True, path_or_buf="dataset/NITRC-ADHD200/phenotypic_data.csv", index=False
)
