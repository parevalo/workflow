import numpy as np
import os

# This script finds all the training cache files in folder, ending with 
# npz extension, and stacks them together to create a single training cache
# file, used to train a single classifier to classify all of the scenes

files = [f for f in os.listdir('.') if f.endswith(".npz")]
X = np.empty([0, 35])

for cache in (files):
    with np.load(cache) as f:
        X_ = f['X']
        y_ = f['y']
        row_ = f['row']
        col_ = f['col']
        labels_ = f['labels']
        
    X = np.vstack((X, X_))  
    y = np.hstack((y, y_))
    row = np.hstack((row, row_))
    col = np.hstack((col, col_))
    labels = np.hstack((labels, labels_))

np.savez('M1_full_traincache.npz', X=X, y=y, row=row, col=col, labels=labels)
    
