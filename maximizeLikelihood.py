###sort
def pivotSort(a,first,last):
#sorting first,last,and middle and choosing the middle as the pivot
    middle = int((first+last+1)/2)
    pivotNumber = middle
    if a[first] > a[last]: a[first], a[last] = a[last], a[first]
    if a[first] > a[middle]: a[first], a[middle] = a[middle], a[first]
    if a[middle] > a[last]: a[last], a[middle] = a[middle], a[last]
    #print('start')
    #print('first %s,last %s'%(first,last))
    #print(a)
    if last-first == 1:return #if the length is of the array is 2, it is already sorted by performing the above process, so no need to continue.
    xP = a[middle]
    #2.looping from the segments first element to its last element, for i, j at the same time.
    i=first
    j=last
    while i<last+1 and j>=first:
        if a[i]<xP:i = i+1
        if a[j]>xP:j = j-1
        if a[i]>=xP and a[j]<=xP:
            if j<=i:break #i,j crossed or saw each other at the same element(pivot), so break.
            else:
                starti = i
                startj = j
                startPivot = pivotNumber
                a[i],a[j]=a[j],a[i]
                #print('startingPivot:%s'%startPivot)
                #print('swapping %s:%s'%(i,j))
                #print(a)
                #print('newi,j:(%s,%s)'%(i,j))
                #if the element being swaped is pivot, take care of the pivot number.
                if starti == startPivot:
                    pivotNumber = j
                    i = i+1
                    #print('i was pivot')
                #if the element being swaped is not pivot, continue.
                elif startj != startPivot:
                    j = j-1
                #if the element being swaped is pivot, take care of the pivot number.
                if startj == startPivot:
                    pivotNumber = i
                    j = j-1
                    #print('j was pivot')
                #if the element being swaped is not pivot, continue.
                if startj != startPivot and starti != startPivot:
                    i = i+1
                #print('newi,j:(%s,%s)'%(i,j))
                #print('pivotNumberNew: %s'%pivotNumber)
    #print('ended the loop,results')
    #print(a)
    #print('pivotNumber: %s'%pivotNumber)
    #3.make sure that we still need to continue the algorithm (we haven't reached the first or last element)
    #sort left part
    if pivotNumber-1>0 and pivotNumber-1>first:
        pivotSort(a,first,pivotNumber-1)
    #sort right part
    if pivotNumber<last-1 and pivotNumber+1<last:
        pivotSort(a,pivotNumber+1,last)







################
###read the data, which will return a list of x

f = open('data/satgals_m14.txt','r')
lines = f.readlines()
f.close()

haloNumber = int(lines[3].split('\n')[0])
haloCount = numpy.zeros(haloNumber)
Xlist = []

haloID = 0
for i in range(4,len(lines)):
    firstRead = lines[i]
    if firstRead == '#\n' and i+1 != len(lines):
        secondRead = lines[i+1]
        if secondRead != '#\n':
            haloCount[haloID] = haloCount[haloID] + 1
            Xlist.append(float(secondRead.split(' ')[0]))
        else: haloID = haloID + 1
    else:
        if i-1 != 3:
            secondRead = lines[i-1]
            if secondRead != '#\n':
                haloCount[haloID] = haloCount[haloID] + 1
                Xlist.append(float(secondRead.split(' ')[0]))
        if i+1 != len(lines):
            nextLine = lines[i+1]
            if nextLine == '#\n':haloID = haloID + 1


####################################
###find a,b,c to maximize likelihood


#first I want to define a function which returns the log-likelihood of a given data realization. for that, I keep in mind that from eqiation (2), I have the density profile, and the likelihood of a realization set, is simply the product of the n(x) for each one of the X = (x_i), so if I take Ln() from both sides, the log likelihood would be the product of each one of the Ln(n(x_i)) for all X = (x_i).

def negativeLogLikelihoodFunction(a,b,c,A,X):
    result = 0.
    for i in range(len(X)):
        l = numpy.log(A)+(a-3.)*numpy.log(X[i]/b)-c*X[i]/b
        result = result - l
    return result


#now I want to calculate this function, for those 6240 point that already calculated the A(a,b,c) for

negLogLikeList = numpy.zeros(((15),(16),(26)))
for i in range(15):
    for j in range(16):
        for k in range(26):
            negLogLikeList[i][j][k] = negativeLogLikelihoodFunction(aArray[i],bArray[j],cArray[k],Alist[i][j][k],Xlist)

#writing a function to minimize another function using Downhill Simplex.
fakeLogLikeList = []
for i in range(15):
    for j in range(16):
        for k in range(26):
            fakeLogLikeList.append(negLogLikeList[i][j][k])
    
pivotSort(fakeLogLikeList,0,len(fakeLogLikeList)-1)
dicX = {}
for n in range(len(fakeLogLikeList)):
    for i in range(15):
        for j in range(16):
            for k in range(26):
                if (fakeLogLikeList[n]-negLogLikeList[i][j][k])**2.<= 0.000000001: dicX[n] = [aArray[i],bArray[j],cArray[k],Alist[i][j][k],negLogLikeList[i][j][k]]

def downHill(fakeLogLikeList,
sumA = 0
sumB = 0
sumC = 0
for i in range(len(X)):
    sumA = sumA + dicA[X[i]]
    sumB = sum
    sumC = sum



