'''linalg/Dumbpy is a simple linear algebra libray built by Joseph Banks'''

import math
import random

def mat_mult(A, B):

    '''Takes two matrices and multiplies them and returns the resulting matrix.

        Parameters:
        A = mxn matrix (type list)
        B = nxr matrix (type list)

        Output:
        AB = mxr matrix (type list)'''

    if len(A[0]) == len(B):
        AB = []
        for i in range(len(A)):
            row = []
            for j in range(len(B[0])):
                element = 0
                for k in range(len(A[0])):
                    element += A[i][k]*B[k][j]
                row.append(element)
            AB.append(row)
        return AB
    else:
        return 'These matrices cannot be multiplied'

def mat_addition(A, B):

    '''Takes two matrices and adds them together.

        Parameters:
        A = mxn matrix (type list)
        B = mxn matrix (type list)

        Output:
        AplusB = mxn matrix (type list)'''

    AplusB = []
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        return 'These matrices cannot be added'
    for i in range(len(A)):
        row = []
        for j in range(len(A[0])):
            element = A[i][j]+B[i][j]
            row.append(element)
        AplusB.append(row)
    return AplusB

def transpose(A):
        
    '''Transposes matrix.

        Parameter:
        A = mxn matrix (type list)

        Output:
        At = nxm matrix (type list)'''

    At = []
    for i in range(len(A[0])):
        row = []
        for j in range(len(A)):
            element = A[j][i]
            row.append(element)
        At.append(row)
    return At

def scale(A, c):
        
    '''Multiplies all of the elements by c.

        Parameter:
        A = mxn matrix (type list)
        c = a scalar (type float or int)

        Output:
        returns the result of multiplying A by an mxm matrix full of zeros with c in the diagonals. This is equivalent to A being
        left multiplied by c times the identity matrix (type list).'''

    cI = []
    for i in range(len(A)):
        cI.append([0]*len(A))
        cI[i][i] = c
    return mat_mult(cI, A)

def vectorize(A):

    '''Takes a list and creates a column vector. This column vector is a list of lists containing a single element within the 
        passed list. This allows the list to be treated as a vector and thus multiplied and augmented with matrices.

        Parameter:
        A = list containing the ordered elements of the desired vector (type list)

        Output:
        vec = list with the elements of A wrapped in individual lists (type list)'''

    vec = []
    for i in A:
        vec.append([float(i)])
    return vec

def unvectorize(A):

    '''Takes a vector and transforms it into a simple python list.
    
        Parameter:
        A = vector that must be in the form outputted by the vectorize function (type list)
        
        Output:
        vecT = list containing the elements in the passed vector (type list)'''
    
    vecT = []
    for i in range(len(A)):
        vecT.append(A[i][0])
    return vecT

def normalize_vector(x):

    '''Take a vector and returns a vector in the same direction with a magnitude of 1

        Parameter:
        x = vector (type list)

        Output:
        returns a new vec that is equivalent to x scaled by 1 divided by the magnitude of x (type list)'''

    norm = 0
    for i in range(len(x)):
        norm += x[i][0]*x[i][0]
    if norm == 0:
        return x
    return scale(x, 1/math.sqrt(norm))

def mat_print(A):

    '''Prints the rows of a matrix on individual lines and an empty line at the end.
    
        Parameter:
        A = mxn matrix (type list) 
        
        Output:
        None'''

    for i in range(len(A)):
        print(A[i])
    print()

def row_swap(A, i, j):

    '''Takes a matrix and the indices of two rows and returns a new matrix with the rows swapped.

        Parameter:
        A = mxn matrix (type list)
        i = index of row (starting from 0) (type int)
        j = index of row (starting from 0) (type int)

        Output:
        B = mxn matrix identical to A except row i and j swapped (type list)'''

    B = A[:]
    temp = B[i]
    B[i] = B[j]
    B[j] = temp
    return B

def row_addition(A, i, j, c):

    '''Takes a matrix and the indices of two rows and a scalar value adds row j scaled by c to row i.

        Parameter:
        A = mxn matrix (type list)
        i = index of row (starting from 0) (type int)
        j = index of row (starting from 0) (type int)
        c = scalar (type float or int)
        
        Output:
        B = mxn matrix with row j scaled by c and added to row i(type list)'''

    B = A[:]
    for k in range(len(B[0])):
        B[i][k] += c*B[j][k]
    return B

def row_elimination(A, i, j):

    '''Takes a matrix and eliminates the equivalent term in i to the leading term of j using row addition.

        Parameter:
        A = mxn matrix (type list)
        i = index of row (starting from 0) (type int)
        j = index of row (starting from 0) (type int)

        Output:
        B = mxn matrix with row i minus j by the ratio of the term in i over the leading term in j (type list)'''
        
    mult = 1
    for k in range(len(A[j])):
        if A[j][k] != 0:
            mult = -A[i][k]/A[j][k]
            row_addition(A, i, j, mult)
            return A
    return A

def row_echelon(A):

    '''Takes a matrix and uses gaussian elimination to reduce it to row echelon form without swap rows. This preserves
        the determinant.

        Parameter:
        A = mxn matrix (type list)

        Output:
        B = mxn matrix in row echelon form (type list)'''

    B = A[:]
    for i in range(len(B)):
        for j in range(i+1, len(B), 1):
            row_elimination(B, j, i)
    # coefficent_to_one(B)
    return B

def rearrange(A):

    '''Rearranges the rows of matrix so that leading terms only have zeros below them in their column. Also known
        as an upper triangular matrix.

        Parameter:
        A = mxn matrix (type list)

        Output:
        A with rows swapped to put it in upper triangular form (type list)'''

    for i in range(len(A)):
        if A[i][i] == 0 and i < len(A)-1:
            A = row_swap(A, i, i+1)
    return A

def coefficent_to_one(A):

    '''Scales all of the rows of the matrix so that the leading coefficent is 1 and then returns the matrix.

        Parameter:
        A = mxn matrix (type list)

        Output:
        A with rows scale (type list)
    '''

    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j] != 0:
                A[i] = scale(vectorize(A[i]), 1/A[i][j])
                A[i] = unvectorize(A[i])
                break
    return A

def rref(A):

    '''Performs guassian elimination to get the matrix into reduced row echelon form and then scales the rows so
        the leading coefficient is 1 and rearranges to terms so they are on the diagonals.

        Parameter:
        A = mxn matrix (type list)

        Output:
        B = mxn matrix in reduced row echelon form (type list)'''

    B = A[:]
    row_echelon(B)
    for i in range(len(B)):
        for j in range(i+1, len(B), 1):
            row_elimination(B, i, j)

    coefficent_to_one(B)
    return rearrange(B)

def determinant(A):

    '''Finds and then returns the determinant of the matrix passed.

        Parameter:
        A = mxn matrix (type list)

        Output:
        determinant = determinant of the matrix (type float or int)'''

    if len(A) != len(A[0]):
        return 'Matrix must be square to have a defined determinant'
    B = row_echelon(A)
    determinant = 1
    for i in range(len(B)):
        determinant *= B[i][i]
    return determinant

def trace(A):

    '''Finds and then returns the trace of the matrix passed.

        Parameter:
        A = mxn matrix (type list)

        Output:
        trace = trace of the matrix (type float or int)'''

    if len(A) != len(A[0]):
        return 'Matrix must be square to have a defined trace'
    trace = 0
    for i in range(len(A)):
        trace += A[i][i]
    return trace

def augment(A, B):

    '''Takes two matrices and returns the 'augmented' matrix. Equivalent to A sitting next to B.

        Parameter:
        A = mxn matrix (type list)
        B = mxr matrix (type list)

        Output:
        AaugB = mx(n+r) matrix (type list)'''

    AaugB = A[:]
    Blen = len(B)
    Alen = len(A)
    if len(A) != len(B):
        return 'They must be row equivalent'
    for i in range(len(B)):
        AaugB[i].extend(B[i])
        B[i] = B[i][:len(B)]
        A[i] = A[i][:len(A)]
    return AaugB

def identity(n):

    '''Creates and returns an nxn Identity matrix.

        Parameter:
        n = the size (type int)

        Output:
        I = nxn identity matrix'''

    I = []
    for i in range(n):
        I.append([0]*n)
        I[i][i] = 1
    return I

def inverse(A):

    '''Finds the inverse of a matrix by augmenting it with the identity matrix and returns the left side
        of the matrix which is the inverse of the original matrix.

        Parameter:
        A = nxn invertible matrix (type list)

        Output:
        inverse = nxn invertible matrix that is the inverse of A (type list)'''

    B = A[:]
    n = len(B)
    if determinant(B) == 0 or len(B) != len(B[0]):
        return 'This matrix is not invertible'
    IaugInverse = augment(B, identity(len(B)))
    for i in range(len(B)):
        A[i] = A[i][:n]
    IaugInverse = rref(IaugInverse)
    inverse = []
    for i in range(len(B)):
        inverse.append(IaugInverse[i][n:])
        A[i] = A[i][:n]
    return inverse

def change_of_basis_matrices(A, B):

    '''Takes two sets of basis vectors in the form of matrices and returns the change of basis matrix 
        from B to A and from A to B.
        
        Parameter:
        A =  nxn invertible matrix (type list)
        B = nxn invertible matrix (type list)
        
        Output:
        BtoA = nxn matrix that is the change of basis matrix from B to A and is invertible (type list)
        inverse(BtoA) = nxn matrix that is the change of basis matrix from A to B and is invertible (type list)'''

    BtoA = rref(augment(B, A))
    for i in range(len(A)):
        BtoA[i] = BtoA[i][len(A):]
    # mat_print(B)
    # mat_print(A)
    return BtoA, inverse(BtoA)

def change_of_coordinates(P, x):

    '''Takes the basis P of a coordinate system and the coordinate vector x returns a new coordinate vector
        in the new coordinate system.

        Parameter:
        P = nxn invertible matrix is basis of new coordinate system (type list)
        x = nx1 coordinate vector (type list)

        Output:
        mat_mult(inverse(P), x) = new nx1 coordinate vector in the new coordinate system (type list)'''

    return mat_mult(inverse(P), x)

def largest_element(x):

    '''Takes a vector and returns the largest element.
    
        Parameter:
        x = nx1 vector (type list)
        
        Output:
        largest = the largest element of the passed vector (type list)'''

    for i in x:
        largest = 0
        if abs(i[0]) > largest:
            largest = i[0]
    return float(largest)

def power_method_eigval(A, x, k):

    '''This method strictly converges to the eigenvalue with the largest absolute value of a matrix and a corresponding eigenvector. 
        http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf is a good resourse to find out about the method used.

        Parameter:
        A = nxn matrix with orthogonal eigenvectors (type list)
        x = nx1 coordinate vector that serves as an initial guess (type list)
        k = the number of iterations the method should use (type int)
        
        Output:
        eigval = approximately the dominant eigenvalue (type float)
        x = a corresponding eigenvector to the found eigenvalue (type list)'''

    for i in range(k):
        x = mat_mult(A, x)
        eigval = largest_element(x)
        u = 1/eigval
        x = scale(x, u)
    return eigval, x

def eig(A, x, k):

    '''This method finds all of the eigenvalues and the corresponding eigenvectors if the eigenvectors are orthogonal (ex. symmetric matrix).
        The eigenvectors must be orthogonal so that we can remove the largest eigenvalues as we find them and replace them with 0 to be
        able to use the power method to find the next eigenvalue and vector.
        http://www.robots.ox.ac.uk/~sjrob/Teaching/EngComp/ecl4.pdf is a is a good resourse to find out about the method used.
        
        Parameter:
        A = nxn matrix with orthogonal eigenvectors (type list)
        x = nx1 coordinate vector that serves as an initial guess (type list)
        k = the number of iterations the method should use (type int)
        
        Output:
        eigvals = a list of the eigenvalues of the matrix from largest to smallest (type list)
        eigvecs = a list of corresponding eigenvectors to the eigenvalues found in eigvals (type list)'''

    B = A[:]
    eigvals = []
    eigvecs = []
    for i in range(len(B)):
        eigval, eigvec = power_method_eigval(B, x, k)
        eigvec = normalize_vector(eigvec)
        eigvals.append(eigval)
        eigvecs.append(eigvec)
        hotel = scale(mat_mult(eigvec, transpose(eigvec)), -eigval)
        B = mat_addition(B, hotel)
    return eigvals, eigvecs

def random_mat(m, n):

    '''Takes the number of desired rows and columns as input and returns a mxn matrix with random integer elements
        in the range of -100 to 100.

        Parameter:
        m = the number of rows in the matrix (type int)
        n = the number of columns in the matrix (type int)

        Output:
        A = mxn matrix with random integer elements in the range of -100 t0 100 (type list)'''

    A = []
    for i in range(m):
        row = []
        for j in range(n):
            element = random.randint(-100,100)
            row.append(element)
        A.append(row)
    return A

def random_sym(n):

    '''Takes the desired size of the matrix and returns a nxn symmertic matrix with random integer elements.
    
        Parameter:
        n = the number of rows/columns in the matrix (type int)
        
        Output:
        mat_mult(transpose(A), A) = nxn symmetric matrix with random integer elements (type list)'''

    A = random_mat(n, n)
    return mat_mult(transpose(A), A)

def random_vec(n):

    '''Takes the length (n) of the desired vector and returns a n-dimensional vector with random integer elements in the 
        range of -5 to 5.
        
        Parameter:
        n = the length or dimension of the vector (type int)
        
        Output:
        vec = n-dimensional vector with random integer elements between -5 and 5 (type list)'''

    vec = []
    for i in range(n):
        element = random.randint(-5,5)
        vec.append([element])
    return vec

def diag(vals):

    '''Takes a list and creates a nxn matrix with the elements of the list as the diagonal entries in the matrix.
    
        Parameter:
        vals = a list of floats or ints (type list)
        
        Output:
        diag = nxn diagonal matrix (type list)'''

    n = len(vals)
    diag = []
    for i in range(n):
        diag.append([0]*n)
        diag[i][i] = vals[i]
    return diag

def vecs_to_mat(vecs):

    '''Takes a set of equal length vectors and returns a matrix where the columns of the matrix are the vectors in the set.
    
        Parameter:
        vecs = set of equal length vectors like that returned by eig() (type list)
        
        Output:
        A = mxn matrix with the vectors as its column (type list)'''

    A = vecs[:]
    for i in range(len(vecs)):
        A[i] = unvectorize(vecs[i])
    return transpose(A)

def eigdecomp(A, k):

    '''Takes a diagonalizable matrix and performs eigendecomposition and returns an orthogonal matrix of eigenvectors, a diagonal
        matrix of the eigenvalues and the transpose/inverse of the eigenvectors matrix.
        
        Parameter:
        A = nxn diagonalizable matrix (type list)
        k = the number of iterations used to find the eigenvalues and eigenvectors (type int)
        
        Output:
        Q = nxn orthogonal matrix of eigenvectors (type list)
        kappa = nxn diagonal matrix of the eigenvalues (type list)
        Qtranspose = transpose/inverse of Q (type list)'''

    eigvals, eigvecs = eig(A, random_vec(len(A)), k)
    kappa = diag(eigvals)
    Q = vecs_to_mat(eigvecs)
    Qtranspose = transpose(Q)
    return Q, kappa, Qtranspose

def singularvals(A):

    '''Takes a mxn matrix and returns its singular values.
    
        Parameter:
        A = mxn matrix (type list)
        
        Output:
        singvals = the singular values of the matrix (type list)'''

    AtA = mat_mult(transpose(A),A)
    eigvals, eigvecs = eig(AtA)
    singvals = []
    for val in eigvals:
        singvals.append(math.sqrt(val))
    
    return singvals

def vec_add(U, V):

    '''Takes two lists that are the same size and adds to the components.
    
        Parameter:
        U = size n standard python list (type list)
        V = size n standard python list (type list)
        
        Output:
        UplusV = size n standard python list with elements equals to U + V element at the index (type list)'''

    UplusV = []
    for u, v in zip(U, V):
        UplusV.append(u + v)
    return UplusV

def mean_dev(X):

    '''Given a observation matrix it will return a matrix in mean deviation form.
    
        Parameter:
        X = mxn matrix (type list)
        
        Output:
        mxn matrix (type list)'''

    M = [0]*len(X)
    for x in transpose(X):
        M = vec_add(M, x)
    M = unvectorize(scale(vectorize(M), -1/len(X[0])))
    return transpose([vec_add(x, M) for x in transpose(X)])

def covar_matrix(X):

    '''Given a observation matrix it will return the covariance matrix of the data.
        
        Parameter:
        X = mxn matrix (type list)
        
        Output:
        mxm symmetrix covariance matrix (type list)'''

    B = mean_dev(X)
    return scale(mat_mult(B, transpose(B)), 1/(len(B[0])-1))

def pca(X, k):

    '''Given a observation matrix it will return the eigvals in descending order and the principal components.
        
        Parameter:
        X = mxn matrix (type list)
        k = number of iterations used to find the eigenvalues of the covariance matrix (type int)
        
        Output:
        eigvals = eigenvalues of the covariance matrix in descending order (type list)
        eigvecs = eigenvectors of the covariance matrix corresponding to the eigvals (type list)'''

    return eig(covar_matrix(X), random_vec(len(X)), k)