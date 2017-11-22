from collections import deque
from collections import Iterable
import numpy as np
from math import log2
#import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_factor
import itertools as it

#TODO rewrite this, no longer accurate
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
class CommObject(object):
    def feed(self, data=()):
        raise NotImplementedError()
    
    def clear(self):
        raise NotImplementedError()
    
    def getinverse(self):
        raise NotImplementedError()
    
    def load(self, data):
        self.buff.add(data)

class Decoder(object):
    def geterrors(self):
        return self.errors.pop(-1)

#turns a sequence of CommObjects into a single CommObject
class Chain(CommObject):
    def __init__(self, comms):
        #validate entries
        for i in comms:
            if not isinstance(comms, CommObject):
                raise ValueError('comms must be an iterable of CommObjects')
        
        self.comms = list(comms)
    
    #feeds data through each CommObject in sequence
    def feed(self, data=()):
        #putsdata through a buffer first because all CommObjects should support load()
        self.load(data)
        out = self.buff.pop(-1)
        
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
    
    def feed(self, data=()):
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
    #TODO rework this to test linear independence of rows (copy, tostandard()?)
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

#implements single-shift register convolutional encoder
#each row of mat represents an output
#columns represent input, delay1, delay2... in that order
#recursive follows the same structure, always start it with a 1
class ConvolutionalEncoder(CommObject):
    def __init__(self, mat, recursive=None):
        if isinstance(mat, np.matrix):
            self.mat = mat
        else:
            self.mat = np.matrix(mat)
        self.shiftlen = len(self.mat.T)-1
        if recursive is None:
            self.recursive = np.matrix([1] + [0]*self.shiftlen)
        else:
            self.recursive = np.matrix(recursive)
        self.reg = deque([0]*self.shiftlen,self.shiftlen)
    
    def feed(self, data=()):
        if not isinstance(data, Iterable):
            data = (data,)
        #TODO replace with empty list of correct size
        out = deque([])
        for i in data:
            state = np.matrix([i]+list(self.reg)).T
            outs = mul(self.mat, state)
            self.reg.appendleft(mul(self.recursive, state)[0,0])
            #turn vertical vector into list
            outs = [i[0] for i in outs.tolist()]
            out.extend(outs)
        return list(out)
    
    def flush(self):
        return self.feed([0]*self.shiftlen)

#relatively slow maximum-likelyhood decoder for convolutional codes
#ViterbiDecoder( convEncoder ) automatically creates a decoder for convEncoder
#ViterbiDecoder( matrix [, recursive [, puncturing ]])
class ViterbiDecoder(CommObject):
    def __init__(self, *args, **kwargs):
        self.buff = Queue()
        if len(args) < 1:
            raise TypeError('class ViterbiDecoder requires at least one argument, received 0')
        #TODO change this logic for new init format
        if isinstance(args[0], ConvolutionalEncoder):
            self.mat = args[0].mat
            self.recursive = args[0].recursive
            self.shiftlen = args[0].shiftlen
            self.block = args[0].block
        else:
            self.mat = np.matrix(args[0])
            if len(args) >=2:
                self.recursive = np.matrix(args[1])
                if not self.recursive:
                    self.recursive = np.matrix([1] + [0]*self.mat.shape[1])
            else:
                self.recursive = np.matrix([1] + [0]*self.mat.shape[1])
        self.shiftlen = self.mat.shape[1]-1
        #trellis forward and reverse links
        self.tf = [None]*self.shiftlen
        self.tr = [[]]*self.shiftlen
        self.datawords = [None] * 2**self.shiftlen
        self.codewords = list(self.datawords)
        for i in range(2**self.shiftlen):
            imat = tobin(i)
            datawords[i] = tobin(i, self.shiftlen+1).T[0]
            codewords[i] = mul(self.mat, datawords[i])
            #bit entering shift register if input is 0
            shiftin = mul(self.recursive[0,1:],imat)
            if shiftin:
                self.tf[i] = [i//2, i//2 + 2**(self.shiftlen-1)]
                self.tr[i//2].append(i)
                self.tr[i//2 + 2**(self.shiftlen-1)].append(i)
            else:
                self.tf[i] = [i//2 + 2**(self.shiftlen-1), i//2]
                self.tr[i//2].append(i)
                self.tr[i//2 + 2**(self.shiftlen-1)].append(i)
        self.dists = [(0,None)] + [(float('inf'),None)]*(2**self.shiftlen - 2)
        #use deque for efficient appends, plus it works without specified block length
        self.trell = deque()
        #self.trellc = deque()
    
    #TODO make all data go through same trellis-building code
    #TODO make single function to build maximum-likelyhood list
    #TODO make trellv only record links forward, store a SINGLE column of net distances elsewhere
    def feed(self, data):
        self.buff.add(data)
        out = []
        while self.block and len(self.buff) + len(self.trellc) >= self.block*self.mat.size[0]:
            l = self.buff.pop(self.block*self.mat.size[0] - len(self.buff) - len(self.trellc))
            blockdata = [l[i:i + n] for i in range(0, len(l), n)]
            for dat in len(blockdata):
                trellcol = [None]*len(self.datawords)
                #TODO implement more efficient euclidean distance (calculate distance from 1/0 beforehand)
                for num in range(len(datawords)):
                    #chooses reverse link with least distance
                    choices = [self.trell[-1][i] for i in self.tr[num]]
                    #TODO implement random decision if reverse links have equal distance
                    choice = self.tr[num][1 * (choices[1] > choices[0])]
                    dist = edist(tobinl(num), blockdata)
                    dist = dist + self.trell[-1][choice]
                    trellcol[num] = (dist, choice)
                self.trell.append(trellcol)
            #TODO make this better
            links = [self.trell[-1][i] for i in self.tr[self.trell[-1].index(min(self.trell[-1]))]]
            statelist = [None]*len(self.block)
            statelist[0] = index(min(self.trell[-1]))
            ind = min(self.trell[-1])[1]
            
            for i, col in it.izip(it.count(), reversed(self.trell)):
                #TODO find a better way to skip first index of trell
                if not i:
                    continue
                statelist[i] = col[ind][1]
                ind = col[ind][1]
            
    
    def clear(self):
        self.buff = Queue()
        trellfirst = [0] + [float('inf')]*(2**self.shiftlen - 2)
        self.trell = deque((trellfirst,))
        #self.trellc = deque()
            

#return behavior is only defined for syndromes with weight <= (d-1)/2
#where d is the min Hamming distance of the generator matrix
class LinearBlockDecoder(CommObject):
    def __init__(self, gen=None, chk=None, error=False):
        if not isinstance(chk, np.matrix):
            chk = np.matrix(chk)
        self.chk = chk
        self.buff = Queue()
        self.codelen = self.gen.shape[1]
        self.datalen = self.codelen - self.gen.shape[0]
    
    def feed(self, data=()):
        self.buff.append(data)
        out = Queue()
        while len(self.buff) >= self.codelen:
            code = np.matrix(self.buff.pop(self.codelen)).T
            synd = self.chk * code
            #TODO finish this

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
        if 'left' in kwargs:
            left = kwargs['left']
        else:
            left = False
        if isinstance(args[0], list):
            #TODO implement polynomial loading
            pass
        else:
            n, k = args
            outerpoly = [1] + [0]*(n-1) + [1]
            factors = gf_factor(ZZ.map(outerpoly), 2, ZZ)
            doms = [len(factor[0])-1 for factor in factors[1]]
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
            if tostand:
                tostandard(self.gen, left)
    
    def getinverse(self, errors=False):
        pass
    
    #TODO more efficiect, cycleic-specific feee() algorithm

#TODO make this
class HammingEncoder(LinearBlockEncoder):
    def __init__(self, paritybits=3):
        pass
        

#basic deque-based FIFO structure for CommObject buffers
class Queue(object):
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
        #-1 pops all, -2 all but one, etc.
        if num < 0:
            num = max(0, len(self.data)+num+1)
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

#DEPRECATED: I'm actually better off not indexing so a deque is fine
#very basic shift register for use in convolutional codes
#size is determined by the size of the iterable passed at initialization
#supports indexing arbitrary values by list or tuple
#pushing data is O(1), indexing is O(1)
class ShiftRegister(object):
    def __init__(self, data):
        self.length = len(data)
        self.data = [i for i in data]
        self.ind = 0
    
    def __len__(self):
        return self.length
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self.data[(key+self.ind)%len(self)]
        elif isinstance(key, (tuple, list)):
            return [self.data[(i+self.ind)%len(self)] for i in key]
        #TODO implement slice indexing
        
    
    def __iter__(self):
        yield [self.data[(i+self.ind)%len(self)] for i in range(len(self))]
    
    def push(self, datum):
        self.ind -= 1
        if self.ind == -1:
            self.ind = len(self)-1
        self.data[self.ind] = datum

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

def tobinl(num, bits=0):
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
    return out

#I stole this from stackoverflow
#https://stackoverflow.com/questions/5488307/numpy-array-in-python-list/
#thanks Sven Marnach!
def arreq_in_list(myarr, list_arrays):
    return any((myarr == x).all() for x in list_arrays)

#adds row a to row b, mod 2, in place
def rowadd(mat, a, b):
    mat[a] += mat[b]
    mat %= 2

#calculates euclidean distance for two iterables
def edist(x1, x2):
    if not len(x1) == len(x2):
        raise ValueError('iterables must have same length')
    dist = 0
    for i in range(x1):
        dist += (x1[i]-x2[i])**2
    return dist**.5

#TODO finish this
#left=True if identity matrix is on the left side
def tostandard(mat, left):
    k, n = mat.shape
    #TODO validate input
    #for i in mat.T:
        #if not i.any():
            #pass
    #first column of identity matrix
    s = (n-k) * (not left)
    for i in range(k):
        #places a 1 in the appropriate cell if necessary
        if not mat[i,s+i]:
            for j in [l for l in range(k) if l != i]:
                if mat[i,s+j]:
                    rowadd(mat, j, i)
                    break
        #places 0s in the rest of the column
        for j in [l for l in range(k) if l != i]:
            if mat[i,s+j]:
                rowadd(mat, i, j)

#0 no; 1 left; 2 right
def isstandard(mat):
    n, k = mat.shape
    #test left
    left = True
    for i in range(k):
        if not (mat[i,i] and mat[:,i].sum() == 1):
            left = False
            break
    if left:
        return 1
    s = n-k
    for i in range(k):
        if not (mat[i,s+i] and mat[:,s+i].sum() == 1):
            return 0
    return 2
