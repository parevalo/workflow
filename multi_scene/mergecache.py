import numpy as np
import os

# This script finds all the training cache files in folder, ending with 
# npz extension, and stacks them together to create a single training cache
# file, used to train a single classifier to classify all of the scenes

# cd to folder
os.chdir("/projectnb/landsat/projects/Colombia/Training/train_cached")

# Find npz files
files = [f for f in os.listdir('.') if f.endswith(".npz")]

# Create empy lists to store the data
Xl, yl, rowl, coll, labelsl = ([] for i in range(5))

# Stack!
for cache in (files):
    with np.load(cache) as f:
        X_ = f['X']
        y_ = f['y']
        row_ = f['row']
        col_ = f['col']
        labels_ = f['labels']
    
    Xl.append(X_)
    yl.append(y_)
    rowl.append(row_)
    coll.append(col_)
    labelsl.append(labels_)

X = np.vstack(Xl)  
y = np.hstack(yl)
row = np.hstack(rowl)
col = np.hstack(coll)
labels = np.hstack(labelsl)

np.savez('M1_full_traincache.npz', X=X, y=y, row=row, col=col, labels=labels)
    
