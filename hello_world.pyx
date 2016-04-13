import numpy as np
cimport numpy as np
from cpython cimport array



def hello_world():
    print "hello world"

def add(int a, int b):
    return a+b

def sum_array(np.ndarray arr):
    cdef array.array arr_c = array.array('d',arr)

    return sum(arr_c.data.as_doubles,arr.size)

def sum_array_2D(np.ndarray arr):

    #cdef array.array arr_c = array.array('d',arr)
    cdef double [:,:] arr_c = arr
    return sum_2D(arr_c,arr.shape[0],arr.shape[1])

def sum_dict(cdef dict d)

    return sum_dict_c(d['a'], d['b'])

cdef int sum_dict_c(int a, int b)

    return a+b

cdef double sum(double* arr, int N):

    cdef int i
    cdef double sum_arr = 0.0

    for i in range(0,N):
        sum_arr +=arr[i]

    return sum_arr

cdef double sum_2D(double[:,:] arr, int M, int N):

    cdef int i
    cdef int j
    cdef double sum_arr = 0.0

    for i in range(0,M):
        for j in range(0,N):
            sum_arr +=arr[i][j]

    return sum_arr