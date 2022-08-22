from multiprocessing import Pool
n=5
q=251

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

exec(open("./hess.sage").read())

def myFunc(val):
    (i,(func,lam))=val
    e=nilpotent(lam.conjugate())
    le=[genMult(J) for J in HessComps(e,func)]
    te=sum([a*b for (a,b) in le])
    ls=[genMult(J) for J in HessComps(semisimple(lam),func)]
    ts=sum([a*b for (a,b) in ls])
    f=open("data.txt",'a')
    f.write('i = '+str(i)+'\n')
    f.write('Semisimple Partition = '+str(lam)+', Nilpotent Partition = '+str(lam.conjugate())+'\nHessenberg Function = '+str(func)+'\n')
    f.write('(Mult, degree) for nilpotent components: '+ str(le)+', '+str(te)+'\n')
    f.write('(Mult, degree) for semisimple components: '+ str(ls)+', '+str(ts)+'\n')
    f.write(str(ts==te)+'\n\n')
    f.close()
    
pool=Pool(processes=4)
funcPartPairs=[(a,b) for a in allHess[n] for b in Partitions(n)[1:-1]]
numFuncPartPairs=[(i,funcPartPairs[i]) for i in range(len(funcPartPairs))]
pool.map(myFunc,numFuncPartPairs)
f=open("data.txt",'a')
f.write('_'*80+'\n\n')
f.close()
