import numpy as np
from sklearn.manifold import TSNE
import h5py

# Data
tcga = h5py.File('../data/tcga_all.h5', mode='r') # Adjust for correct file
cancers = list(tcga['tcga/train'])
tcga_stack = np.vstack(list([tcga['tcga/train/'+c] for c in cancers]))

# Save cancer labels
labels = np.array([c for c in cancers for i in range(tcga['tcga/train/'+c].shape[0])])
np.save('labels.npy',labels)

# TSNE
Z = TSNE(n_components=2, verbose=1).fit_transform(tcga_stack)

np.save('Z.npy',Z)