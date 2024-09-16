import sys
import numpy as np
from kmeans import kmeans
import symnmfmodule
from sklearn.metrics import silhouette_score

from symnmf import initialize_H
from symnmf import read_file

#Retruns the closest centroid for given vector
def closest_centroid(vector, centroids):
    distances = np.linalg.norm(centroids - vector, axis=1)
    return np.argmin(distances)

def calculate_kmeans_clusters(X, k):
    centroids = kmeans(X, k)
    return [closest_centroid(vector, centroids) for vector in X]

#Calculate the symnmf clusters as mentioned in 1.5
def calculate_symnmf_clusters(X, k):
    X = X.tolist()
    W = symnmfmodule.norm(X)
    H = initialize_H(W,k)
    symnmf_result = symnmfmodule.symnmf(H,W)
    H = np.array(symnmf_result)
    return H.argmax(axis=1)


def main():
    np.random.seed(0) 
    #Get args
    k = int(sys.argv[1]) 
    filename = str(sys.argv[2])

    X,rows,cols= read_file(filename) 
    
    X = np.array(X) #read_file returns a List therefore parse to np.array to be competable to othen methods 
    if X is None:
        print("An Error Has Occurred")
        return

    kmeans_cluster = calculate_kmeans_clusters(X, k)
    symnmf_cluster = calculate_symnmf_clusters(X, k)

    print("nmf: %.4f" % silhouette_score(X, symnmf_cluster))
    print("kmeans: %.4f" % silhouette_score(X, kmeans_cluster))

if __name__ == "__main__":
    main()