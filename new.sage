#from multiprocessing import Pool
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ring import Ring
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.matrix.constructor import Matrix
#from sage.matrix.constructor import matrix
from sage.matrix.special import block_diagonal_matrix
from sage.misc.prandom import randrange
from sage.rings.ideal import Ideal
from sage.combinat.partition import Partitions
from itertools import permutations, accumulate

n=4
q=257

allHess={}
allHess[1]=[[1]]
for k in range(2,n+1):
    allHess[k]=[]
    for func in allHess[k-1]:
        allHess[k].extend([[i+1]+[j+1 for j in func] for i in range(1+func[0])])

filename='data.txt'

#for n in range(2,N):
#    exec(open("./hess.sage").read())
#    file=open(filename,'a')
#    file.write('-'*80+'\n\n')
#    file.close()

#exec(open("./hess.sage").read())
#varlist=[]
#for i in range(1,n+1):
#   for j in range(1,n+1):
#       varlist.append('x'+str(i)+str(j))
varlist=['x%d%d'%(i+1,j+1) for i in range(n) for j in range(n)]
A = PolynomialRing(FiniteField(q),varlist)

g=Matrix(A,n,A.gens())
det=g.det()

l=list(A.gens())
l.reverse()
#x=dict({})
x=dict({(i+1,j+1):l.pop() for i in range(n) for j in range(n)})

def semisimple(lam): 
    return matrix.diagonal(sum([[i]*lam[i] for i in range(lam.length())],[]))

def oneBlock(k):
    l=[]
    for i in range(k-1): l.extend([0]*(i+1)+[1]+[0]*(k-i-2))
    l.extend([0]*k)
    return Matrix(k,l)

def nilpotent(partition): return Matrix(list(block_diagonal_matrix([oneBlock(k) for k in partition])))

Matrix(list(block_diagonal_matrix([oneBlock(k) for k in lam])))

def randLin(): return randrange(q)+homRandLin()

def homRandLin(): return sum(randrange(q)*y for y in A.gens())

def genMult(I):
'''
only works for I a primary ideal
'''
    d=I.dimension()
    J=I+[homRandLin() for i in range(d-1)]
    while not(J.dimension()==1): J=I+[homRandLin() for i in range(d-1)]
    K=J+randLin()
    while not(K.dimension()==0): K=J+randLin()
    deg=K.radical().vector_space_dimension()
    return (K.vector_space_dimension()/deg,deg)

def HessComps(x,func):
    I=Ideal(A,0)
    for k in range(1,n):
        h=Matrix(A,[g.columns()[i] for i in range(func[k-1])]+[(x*g).column(k-1)]) 
        I=I+h.minors(func[k-1]+1)
    l=[(J,J.dimension()) for J in I.primary_decomposition() if det not in J.radical()]
    maxDim=max(b for (_,b) in l)
    l=[x for x in l if x[1]=maxDim]
    return l
#    return([(J,J.dimension()) for J in I.primary_decomposition() if det not in J.radical()])

def springerComps(x): return HessComps(x,[1..x.ncols()])
#    I=Ideal(A,0)
#    g=Matrix(A,n,A.gens())
#    for k in range(1,n):
#        h=Matrix(A,[g.columns()[i] for i in range(k)]+[(s*g).column(k-1)]) 
#        I=I+h.minors(k+1)
#    return([J for J in I.primary_decomposition() if det not in J.radical()])

#file=open(filename,'a')
#file.write('n='+str(n)+'\n')
#file.close()
#
#for lam in Partitions(n)[1:]:
#    file=open(filename,'a')
#    file.write('Nilpotent partition='+str(lam.conjugate()))
#    file.write(', Semisimple partition='+str(lam))
#    file.write('\n')
#    file.close()
#    for func in allHess[n]:
#        file=open(filename,'a')
#        file.write('Hessenberg function = '+str(func)+'\n')
#        file.close()
#        e=nilpotent(lam.conjugate())
#        file=open(filename,'a')
#        le=[genMult(J) for J in HessComps(e,func)]
#        te=sum([a*b for (a,b) in le])
#        file.write('(Mult, degree) for nilpotent components: ')
#        file.write(str(le)+', '+str(te))
#        file.write('\n')
#        file.close()
#        ls=[genMult(J) for J in HessComps(semisimple(lam),func)]
#        ts=sum([a*b for (a,b) in ls])
#        file=open(filename,'a')
#        file.write('(Mult, degree) for semisimple components: '+ str(ls)+', '+str(ts)+'\n')
#        file.write(str(te==ts)+'\n\n')
#        file.close()
#    file=open(filename,'a')
#    file.write('-'*80+'\n')
#    file.close()

def printDimensions(func,lam):
    e=nilpotent(lam)
    I=Ideal(A,0)
    for k in range(1,n):
        h=Matrix(A,[g.columns()[i] for i in range(func[k-1])]+[(e*g).column(k-1)]) 
        I=I+h.minors(func[k-1]+1)
#    print(func,lam,I.dimension()-n*(n+1)/2,conjDimension(func,lam))
    return (I.dimension()-(n*(n+1)/2)==conjDimension(func,lam))
    
def springerDimension(lam):
#    e=nilpotent(lam)
    func=[i for i in range(1,n+1)]
    return max(J.dimension() for J in HessComps(nilpotent(lam),list(range(1,n+1))))
#    I=Ideal(A,0)
#    for k in range(1,n):
#        h=Matrix(A,[g.columns()[i] for i in range(k)]+[(e*g).column(k-1)]) 
#        I=I+h.minors(k+1)
#    return I.dimension()-(n*(n+1)/2)
    
def conjDimension(func,lam):
    acc=[list(accumulate(comp))[:-1] for comp in permutations(lam)]
    return springerDimension(lam)+ max(sum(func[i-1]-i for i in range(1,n+1) if i not in partial) for partial in acc)

for func in allHess[n]:
    for lam in Partitions(n):
        if not printDimensions(func,lam): print (func,lam)

def myFunc(val):
    (i,(func,lam))=val
    e=nilpotent(lam.conjugate())
    le=[genMult(J) for (J,_) in HessComps(e,func)]
    te=sum([a*b for (a,b) in le])
    ls=[genMult(J) for (J,_)in HessComps(semisimple(lam),func)]
    ts=sum([a*b for (a,b) in ls])
    f=open("data.txt",'a')
    f.write('i = '+str(i)+'\n')
    f.write('Semisimple Partition = '+str(lam)+', Nilpotent Partition = '+str(lam.conjugate())+'\nHessenberg Function = '+str(func)+'\n')
    f.write('(Mult, degree) for nilpotent components: '+ str(le)+', '+str(te)+'\n')
    f.write('(Mult, degree) for semisimple components: '+ str(ls)+', '+str(ts)+'\n')
    f.write(str(ts==te)+'\n\n')
    f.close()
#    print('i = '+str(i)+'\n')
#    print('Semisimple Partition = '+str(lam)+', Nilpotent Partition = '+str(lam.conjugate())+'\nHessenberg Function = '+str(func)+'\n')
#    print('(Mult, degree) for nilpotent components: '+ str(le)+', '+str(te)+'\n')
#    print('(Mult, degree) for semisimple components: '+ str(ls)+', '+str(ts)+'\n')
#    print(str(ts==te)+'\n\n')
    
#pool=Pool(processes=2)
funcPartPairs=[(a,b) for a in allHess[n] for b in Partitions(n)[1:-1]]
numFuncPartPairs=[(i,funcPartPairs[i]) for i in range(len(funcPartPairs))]
#numFuncPartPairs=[(i,funcPartPairs[i]) for i in range(13,18)]
#pool.map(myFunc,numFuncPartPairs[1:])
for x in (numFuncPartPairs[1:]): myFunc(x)

f=open("data.txt",'a')
f.write('_'*80+'\n\n')
f.close()
