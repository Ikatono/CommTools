from collections import deque
from collections import Iterable
import numpy as np
from math import log2
#import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_factor

#basic structure of objects used in this package
#feed() is used to send data to the object
#it stores whatever state data it needs to determine future bits
#feed() returns everything new up to (not including) the first undetermined bit
#for example, a block encoder will return every newly completed block and store part of the next
#these bits must be used or stored immediately and cannot be received again
#clear() removes all state data fed to the object and returns it to initialization state
#All CommObjects allow real numerics and Booleans. Decoders will return real numerics.
#Hard decoders return 1 and -1. True and False are interpretted as 1 and -1 respectively
#some objects allow other data types, but they may behave strangely if you use Iterables as single data points
#because of this, DO NOT USE STRINGS AS DATA POINTS
class CommObject:
    def feed(self):
        raise NotImplementedError()
    
    def clear(self):
        raise NotImplementedError()
    
    def getinverse(self):
        raise NotImplementedError()

#turns a sequence of CommObjects into a single CommObject
class Chain(CommObject):
    def __init__(self, comms):
        #validate entries
        for i in comms:
            if not isinstance(comms, CommObject):
                raise ValueError('comms must be an iterable of CommObjects')
        
        self.comms = list(comms)
    
    #feeds data through each CommObject in sequence
    def feed(self, data):
        if not isinstance(data, Iterable):
            out = (data,)
        else:
            out = data
        
        for i in self.comms:
            out = i.feed(out)
        
        return out
    
    def clear(self):
        [i.clear() for i in self.comms]
    
    #creates inverses for each component and places them in reverse order
    def getinverse(self):
        ln = len(self.comms)
        inv = [None] * ln
        for i in range(ln):
            inv[i] = self.comms[ln-i-1]
        return Chain(inv)

class LinearBlockEncoder(CommObject):
    def __init__(self, mat):
        if not isinstance(mat, np.matrix):
            self.gen = np.matrix(mat)
        else:
            self.gen = mat
        self.datalen = self.gen.shape[0]
        self.codelen = self.gen.shape[1]
        self.buff = Queue()
        self.dist = self.distance(self.gen)
    
    def feed(self, data):
        self.buff.add(data)
        out = []
        while len(self.buff) >= self.datalen:
            datavect = np.matrix(self.buff.pop(self.datalen))
            codevect = mul(datavect,self.gen)
            out.append(codevect.tolist()[0])
        return out
    
    def clear(self):
        self.buff = Queue()
            
    def checkmat(self):
        if not self.standard(self.gen):
            raise ValueError('Currently only standard matrices can automatically generate check matrices')
            raise ValueError('Currently only standard matrices can automatically generate check matrices')
        (h,w) = self.gen.shape
        return np.concatenate((self.gen[:,h:].T, np.identity(w-h)), 1)
    
    def getinverse(self):
        return LinearBlockDecoder(self.checkmat())
    
    @staticmethod
    def distance(gen):
        datalen = gen.shape[0]
        datawords = [tobin(i,datalen) for i in range(2**datalen)]
        codewords = [mul(i.T, gen) for i in datawords]
        return min(int(sum(i.T)) for i in codewords if int(sum(i.T))>0)
    
    @staticmethod
    def standard(mat):
        if not isinstance(mat, np.matrix):
            mat = np.matrix(mat)
        h = mat.shape[0]
        return np.array_equal(mat[:,:h],np.identity(h))
        
    #literally useless, why did I make this?
    #LINEAR ALGEBRA MAKES LINEAR CODES!
    @staticmethod
    def validate(mat):
        if not isinstance(mat, np.matrix):
            mat = np.matrix(mat)
        #list of columns in mat
        #cols = [i.T for i in mat.T]
        #verify that every sum of 2 columns yields another column
        #TODO see if there's a more efficient way to do this
        datawords = [tobin(i,mat.shape[0]) for i in range(2**mat.shape[0])]
        #print(datawords)
        codewords = [mul(i.T, mat) for i in datawords]
        for c1 in range(len(codewords)):
            for c2 in range(c1, len(codewords)):
                if not arreq_in_list(add(codewords[c1],codewords[c2]),codewords):
                    return False
        return True

#return behavior is only defined for syndromes with weight <= (d-1)/2
#where d is the min Hamming distance of the generator matrix
class LinearBlockDecoder(CommObject):
    def __init__(self, chk, error=False):
        if not isinstance(chk, np.matrix):
            chk = np.matrix(chk)
        self.chk = chk
        self.buff = Queue()
        self.codelen = self.gen.shape[1]
        self.datalen = self.codelen - self.gen.shape[0]
    
    def feed(self, data):
        self.buff.append(data)
        while len(self.buff) >= self.codelen:
            code = np.matrix(self.buff.pop(self.codelen)).T
            synd = self.chk * code
            #data = 
#creates a cyclic encoder, a type of linear block code
#can be used in two different ways:
#CyclicEncoder(n, k) creates an [n,k] code from a code polynomial of max distance
#raises a ValueError if no such polynomial exists
#CyclicEncoder(poly [, n]) creates a code using poly as the code polynomial
#poly is a list of 1 and 0, where poly[i] is the coefficient of X**i
#if n is not defined then it is infered as len(poly)-1
#raises a ValueError if poly is invalid for this order of code
class CyclicEncoder(LinearBlockEncoder):
    def __init__(self, *args, **kwargs):
        #TODO validate inputs better
        if 'standard' in kwargs:
            tostand = kwargs['standard']
        else:
            tostand = True
        if isinstance(args[0], list):
            #TODO implement polynomial loading
            pass
        else:
            n, k = args
            outerpoly = [1] + [0]*(n-1) + [1]
            factors = gf_factor(ZZ.map(outerpoly), 2, ZZ)
            doms = [len(factor[0])-1 for factor in factors[1]]
            #print(factors)
            #print(doms)
            if n-k not in doms:
                raise ValueError('there is no cyclic code with these values')
            #list of all factors of the correct order
            facts = [i for i, x in enumerate(doms) if x == n-k]
            facts = [factors[1][i][0] for i in facts]
            dists = [0] * len(facts)
            gens = [None] * len(facts)
            for i in range(len(facts)):
                gen = np.matrix(np.zeros((k,n)))
                l = len(facts[i])
                for j in range(k):
                    poly = np.matrix([0]*j + facts[i][::-1] + [0]*(n-j-k))
                    gen[j] = poly
                gens[i] = gen
                dists[i] = self.distance(gen)
            #chooses arbitrary factor of maximum distance
            choice = dists.index(max(dists))
            self.poly = facts[choice]
            self.dist = dists[choice]
            self.gen = gens[choice]
            self.buff = Queue()
            self.datalen = self.gen.shape[0]
            self.codelen = self.gen.shape[1]
            

class HammingEncoder(LinearBlockEncoder):
    def __init__(self, paritybits=3):
        pass
        

#basic queue structure for CommObject buffers
class Queue:
    def __init__(self, data=()):
        self.data = deque(data)
    
    #TODO add __blah__ functions
    #or don't, this shouldn't be used outside this module
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, key):
        return self.data[key]
    
    def add(self, data):
        if isinstance(data, Iterable):
            for i in data:
                self.data.append(i)
        else:
            self.data.append(data)
    
    def pop(self, num=1):
        if num > len(self.data):
            raise IndexError('Cannot pop more values than are in queue')
            
        out = [0] * num
        for i in range(num):
            out[i] = self.data.popleft()
        return out

def add(x1, x2):
    if not isinstance(x1, np.matrix):
        x1 = np.matrix(x1)
    if not isinstance(x2, np.matrix):
        x2 = np.matrix(x2)
    return (x1+x2) % 2

def mul(x1, x2):
    if not isinstance(x1, np.matrix):
        x1 = np.matrix(x1)
    if not isinstance(x2, np.matrix):
        x2 = np.matrix(x2)
    return x1.dot(x2) % 2

#converts each value to 1 (if True or >0) or 0
#for use with hard decoders
#TODO finish this
def bininput(x):
    if isinstance(x, np.matrix):
        #TODO make this less bad
        return np.matrix([[1 if elem else 0 for elem in row] for row in x])
    
    

def tobin(num, bits=0):
    if num == 0:
        return np.matrix([0] * max(bits,1)).T
    if not isinstance(num, int):
        raise TypeError('tobin only takes integers')
    bitsmax = int(log2(num))+1
    if not bits:
        bits = bitsmax
    if bitsmax > bits:
        raise ValueError('insufficient bits')
    out = [0] * bits
    for i in range(bitsmax-1,-1,-1):
        if num >= 2**i:
            out[i] = 1
            num -= 2**i
    return np.matrix(out).T

#I stole this from stackoverflow
#https://stackoverflow.com/questions/5488307/numpy-array-in-python-list/
#thanks Sven Marnach!
def arreq_in_list(myarr, list_arrays):
    return any((myarr == x).all() for x in list_arrays)

#adds row a to row b, mod 2
def rowadd(mat, a, b):
    mat[a] += mat[b]
    mat %= 2

#TODO finish this
#left=True if identity matrix is on the left side
def standard(mat, left):
    k = mat.shape[1]
    for i in mat.T[:k]:
        if not i.any():
            pass
            
