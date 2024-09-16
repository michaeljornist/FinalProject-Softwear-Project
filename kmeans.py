import numpy as np
import sys
class Cluster:
    EPSILON = 0.0001
    def __init__(self, centroid,d):
        self.centroid = np.copy(centroid) 
        self.prevCentroid = np.copy(centroid) 
        self.d = d
        self.points = np.empty((0,self.d)) # saves the indexs of the vectors from the matrix
        
    def update_centroid(self):
        # Calculate the mean of the points along each dimension to get the new centroid
        self.prevCentroid = np.copy(self.centroid)
        self.centroid = np.mean(self.points, axis=0)
        
        self.points = np.empty((0,self.d))  # Reset points for the next iteration

    def getCentroid(self):
        return self.centroid
    

    def add_point(self, vector):
        self.points = np.vstack([self.points, np.copy(vector)])

    def has_converged(self):
        return np.linalg.norm(self.centroid - self.prevCentroid) < Cluster.EPSILON
    
    def __str__(self):
        return f'the centroid is \n {self.centroid} \n the prevCentroid is \n {self.prevCentroid} \n ,the points are \n  {self.points}'

def checkInput(argv):
    if len(argv)>6 or len(argv)<5:
        raise Exception('An Error Has Occurred')
    k=int(argv[1])if argv[1].isdigit() else None
    N=int(argv[2])if argv[2].isdigit() else None
    d=int(argv[3])if argv[3].isdigit() else None
    
    if N==None or N<=1:
        raise Exception('Invalid number of points!')
    if k==None or (k<=1 or k>=N):
        raise Exception('Invalid number of clusters')
    if d==None or d<1:
        raise Exception('Invalid dimension of point')
    if(len(argv)!=5):
        iter=int(argv[4])if argv[4].isdigit() else None
        if iter==None or (iter<=1 or iter>=1000):
            raise Exception('Invalid maximum iteration')
        if argv[5].lower().endswith('.txt')==False:
            raise Exception('NA')
    else:
        if argv[4].lower().endswith('.txt')==False:
            raise Exception('NA')

def addPointsToClusters(matrix,clusterList,N,k):
    for i in range(0,N):
        bestDis = np.linalg.norm(matrix[i,:] - clusterList[0].getCentroid())
        bestIndxOfCluster = 0
        for j in range(0,k):
            distance = np.linalg.norm(matrix[i,:] - clusterList[j].getCentroid())
        
            if distance < bestDis:
                bestDis = distance
                bestIndxOfCluster = j
        
        clusterList[bestIndxOfCluster].add_point(matrix[i,:])

def printMat(mat):
    i=0
    j=0
    for i in range(mat.shape[0]):  
        for j in range(mat.shape[1]):
            if j!=(mat.shape[1]-1):
                formatted = f"{mat[i,j]:.4f}"
                print(formatted,end=',')
            else:
                formatted = f"{mat[i,j]:.4f}"
                print (f'{formatted}')
 


# in HM1 this was the main , slightly changed to work with analysis.py
#Receives the data matrix and number of centroids (k) , returns the centroieds.
def kmeans(matrix,k):
    N, d = matrix.shape
    # initialition 
    clusterList = []
    for i in range(0,k):
        cluster = Cluster(matrix[i,:],d)
        clusterList.append(cluster)

    iter_number = 0 

    while iter_number < 300:
        addPointsToClusters(matrix,clusterList,N,k)

        flag = True
        for i in range(0,k):
            clusterList[i].update_centroid()
            if not clusterList[i].has_converged():
                flag = False
        
        if flag:
            break

    centroidMatrix = np.empty((0,d))
    for i in range(0,k):
        centroidMatrix = np.vstack([centroidMatrix, clusterList[i].getCentroid()])
    
    return centroidMatrix


