def reducedSpringerIdeal(s):
    I=Ideal(A,0)
    g=Matrix(A,n,A.gens())
    for k in range(1,n):
        h=Matrix(A,[g.columns()[i] for i in range(k-1)]+[(s*g).column(k-1)]) 
        I=I+h.minors(k)
    return([J for J in I.primary_decomposition() if det not in J.radical()])


X = SchubertPolynomialRing(ZZ)
X.one()

s1 = Permutation([2,1,3,4])
s2 = Permutation([1,3,2,4])
s3 = Permutation([1,2,4,3])
w0 = Permutation([4,3,2,1])
p1 = Permutation([2,3,4,1])  #s1s2s1

p2 = Permutation([2,3,1,4])  #1213
p2 = Permutation([2,3,1,4])  #1213


R = RootSystem(['A',3]).root_space()
W = WeylGroup(R,prefix="s") #,implementation='coxeter3')
s = W.simple_reflections()

X([2,3,4,1])+ X([1,2,4,3])*X([2,3,1,4])+X([2,1,3,4])*X([1,4,2,3])+X([4,1,2,3])

n=8
y=0
def wL(x):
    if (x<n):return (n-x)
    else : return(x)

for i in [0..n-1]:
    wi=[1..i]+[n]+[i+1..n-1]
    w0wLwi=[n+1-i for i in [wL(j) for j in wi]]
    print (wi)
    print(w0wLwi,'\n')
    y=y+X(w0wLwi)*X(wi)
print (y)

