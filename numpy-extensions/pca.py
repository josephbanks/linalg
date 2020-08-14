import numpy as np
import matplotlib.pyplot as plt

class PCA:
    def __init__(self, X):
        self.data = X
        self.cov = self.covar_matrix()
        self.eigvals, self.pcs = self.pca()
        self.importance = self.importance()
    def pca(self):
        return self.eig()
    def eig(self):
        eigvals, eigvecs = np.linalg.eigh(self.cov)
        index = eigvals.argsort()[::-1]
        eigvals = eigvals[index]
        eigvecs = eigvecs[:,index]
        return eigvals.real, eigvecs.real
    def covar_matrix(self):
        B = self.mean_dev()
        return np.matmul(B.transpose(),B)/(self.data.shape[0])
    def mean_dev(self):
        M = np.mean(self.data, axis=0)
        M = np.tile(M,(self.data.shape[0],1))
        return np.subtract(self.data, M)
    def closest_approx(self, y, U):
        UUt = np.dot(U, U.transpose())
        return np.dot(UUt, y)
    def importance(self):
        tr = np.sum(self.eigvals)
        return np.divide(self.eigvals,tr)
    def necessary_dims(self, entropy):
        entr = 0
        dim = 0
        per_var = self.importance
        while entr < entropy:
            entr += per_var[dim]
            dim += 1
        return dim
    def weights(self,dim):
        W = np.zeros((self.data.shape[0],dim))
        for j in range(self.data.shape[0]): 
            weights = np.zeros(dim)
            for i in range(dim):
                weights[i] = np.dot(self.data[j], self.pcs[:,i])
            W[j] = weights
        return W
    def reduce_dim(self, dim):
        U = self.pcs[:,:dim]
        return np.dot(U, U.transpose())

class PCA_svd:
    def __init__(self, X):
        self.X = X
        self.eigvals, self.pcs = self.pca()
        self.importance = self.importance()
    def pca(self):
        U, S, Vh = np.linalg.svd(self.mean_dev())
        eigvals = np.square(S)/(self.X.shape[0]-1)
        pcs = Vh.transpose()
        return eigvals, pcs
    def mean_dev(self):
        M = np.mean(self.X, axis=0)
        M = np.tile(M,(self.X.shape[0],1))
        return np.subtract(self.X, M)
    def closest_approx(self, y, U):
        UUt = np.dot(U, U.transpose())
        return np.dot(UUt, y)
    def importance(self):
        tr = np.sum(self.eigvals)
        return np.divide(self.eigvals,tr)
    def necessary_dims(self, entropy):
        entr = 0
        dim = 0
        per_var = self.importance
        while entr < entropy:
            entr += per_var[dim]
            dim += 1
        return dim
    def weights(self,dim):
        W = np.zeros((self.data.shape[0],dim))
        for j in range(self.data.shape[0]): 
            weights = np.zeros(dim)
            for i in range(dim):
                weights[i] = np.dot(self.data[j], self.pcs[:,i])
            W[j] = weights
        return W
    def reduce_dim(self, dim):
        U = self.pcs[:,:dim]
        return np.dot(U, U.transpose())

def weights(y,Basis,size):
    weights = np.zeros(size)
    for i in range(size):
        weights[i] = np.dot(y, Basis[:,i])
        print(weights[i])
    return weights

def show_eigenimages(height,width,PCA):
    for i in range(height*width):
        plt.subplot(height,width,i+1)
        plt.imshow(PCA.pcs[:,i].reshape(62,47))