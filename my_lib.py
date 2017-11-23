import numpy as np
import graph_lib

def get_G_H_matrixes():
    # returns G and H matrixes in a tuple

    A = np.eye(12, dtype=np.int);
    B = np.array([[1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1],
                  [0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],
                  [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1],
                  [1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0],
                  [1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1],
                  [1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0],
                  [1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1],
                  [1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0],
                  [0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],
                  [0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0],
                  [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                  [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1]], dtype=np.int);
    G = np.concatenate((A, B), axis=1);
    H = np.concatenate((B, A), axis=0);
    return (G, H);

def bin_sum(a, b):
    n = max(a.shape[0], b.shape[0]);
    c = a+b;
    for i in np.arange(n-1, 0, -1):
        if(c[i]==2):
            c[i]=0;
            c[i-1]+=1;
        if(c[i]==3):
            c[i]=1;
            c[i-1]+=1;
    return c;

def bin_or(a, b):
    n = max(a.shape[0], b.shape[0]);
    c = np.ndarray(shape = (n,), dtype=np.int);
    for i in np.arange(0, n, 1):
        c[i] = (a[i] or b[i])
    return c;

def bin_multipl(A, B):
    shape_a = A.shape;
    shape_b = B.shape;
    C = np.ndarray(shape = (shape_a[0], shape_b[1]), dtype=np.int);

    if(shape_a[1]!=shape_b[0]):
        return 0;

    for i in np.arange(0, shape_a[0], 1):
        for j in np.arange(0, shape_b[1], 1):
            C[i][j] = 0;
            for k in np.arange(0, shape_a[1], 1):
                C[i][j] += A[i][k]*B[k][j];
            C[i][j] = C[i][j] % 2;
    return C;

def bin_weight(a):
    weight = 0;
    for i in np.arange(0, a.shape[0], 1):
        if (a[i] == 1):
              weight += 1;
    return weight;

def hamming_dist(a, b):
    dist = 0;
    for i in np.arange(0, a.shape[0], 1):
        if(a[i] != b[i]):
            dist += 1;
    return dist;

def get_bin_nums_matrix():
    # matrix of all possible binary words with length 12
    bin_one = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype=np.int);

    B = np.ndarray(shape=(4096, 12), dtype=np.int);
    B[0] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    for i in np.arange(0, 4095, 1):
        next = bin_sum(B[i], bin_one);
        B[i+1] = next;
    return B;

def get_codewords(Din, G):
    Code_words = np.ndarray(shape = (Din.shape[0], 24), dtype=np.int);
    Code_words = bin_multipl(Din, G);
    return Code_words;

def get_words_of_weight_eight(codewords):
    codewords_w_eight = np.ndarray(shape=(759,24), dtype=np.int);
    j = 0;
    for i in np.arange(0, codewords.shape[0], 1):
        if(bin_weight(codewords[i]) == 8):
            codewords_w_eight[j] = codewords[i];
            j += 1;
    return codewords_w_eight;

def get_connectivity_matrix(codewords):
    connectivity_matrix = np.zeros((codewords.shape[0], codewords.shape[0]), dtype=np.int);
    for i in np.arange(0, codewords.shape[0], 1):
        for j in np.arange(0, codewords.shape[0], 1):
            if((hamming_dist(codewords[i], codewords[j]) >= 12) or (hamming_dist(codewords[i], codewords[j]) == 0)):
                connectivity_matrix[i][j] = 1;
    return connectivity_matrix;

def get_connectivity_list(codewords):
    connectivity_list = list([]);

    for i in np.arange(0, codewords.shape[0], 1):
        line = list([]);
        for j in np.arange(0, codewords.shape[0], 1):
            if(hamming_dist(codewords[i], codewords[j]) >= 12):
                line.append(j);
        connectivity_list.append(line);
    return connectivity_list;
