import numpy
from matplotlib import pyplab as plot


####################################
###producing random values for a,b,c


a = XORshift()*1.4 + 1.1
b = XORshift()*1.5 + 0.5
c = XORshift()*2.5 + 1.5


######################################################################
###defining the function that I want to integrate, in order to find a.


def iWantToIntegrateThis(x):
    global a
    global b
    global c
    if x == 0: y = 0
    else:
        y = 4. * numpy.pi * (x/b)**(a-3.) * numpy.exp(-(x/b)**c) * x**2.
    return y


################################################################
###building up the numerical integrator, using Romberg algorithm


#I'll make a table, containg all the steps taken to the final estimation of the integral (process)
#1.simplest estimation, using just one trapezoid.
#2.doubling the number of trapezoids in each step, and since we already know the values at half of these points, calculating the other half. so we are making better and better estimations by increasing the number of trapizoids, but, this is getting better linearly, and very slow. 
#3.we can make a better estimation, and very fast, using romberg integration. for that, solving the romberg equation, we use the previously calculated values, and make a better estimation (practically we are making our estimation better, by using the fact that we already know our error estimations).
#this return the whole table and the best estimated value which is the [-1,-1] member of the table.

def integration(f,a,b,steps):
    process = numpy.zeros((steps,steps),dtype=numpy.float64)
    #1
    process[0,0] = ( f(a)+f(b) )*(b-a)/2.
    #2
    for i in range(steps-1):
        process[i+1,0] = (b-a)*sum( f(a + (b-a)*(2.*j+1)/(2.**(i+1)) ) for j in range(2**i) )
        process[i+1,0] = process[i,0]/2. + process[i+1,0]/(2.**(i+1))
        #3
        for k in range(i+1):
            process[i+1,k+1] = process[i+1,k] + (process[i+1,k]-process[i,k])/(4**(k+1) -1.)

    return process #this return the whole table and the best estimated value which is the [-1,-1] member of the table.


###########################################################################################################################
###now using the defined function and the numerical integrator, I can integrate that function, hence the value of A(a,b,c).


def findA():
    denom = integration(iWantToIntegrateThis,0.,5.,5)
    return 1./denom[-1,-1]


######################################
###defining the p(x) function (part d)


def p(x,A,a,b,c):
    y = A*4.*numpy.pi*(x**2.)*(x/b)**(a-3) * numpy.exp(-(x/b)**c)
    return y

#########################################################
###here I'm using p(x) function to produce 100 satellites 


#I produce 1500 numbers for x, and 1500 numbers for y, of course different numbers! but, since I want to plot a log-log at the end, the idea is to produce x logarithmically, so I give a chance to the smaller x to appear in the plot.
#for each x, I'm producing a y, taking this y as the random probability for that x, if y is smaller than f(x) it can be a member of the final 100 galaxies, else, I throw it away.

n = 0
pointsList = numpy.zeros(100)
while n< 100:
    #generate a x in the logarithmic scale, and then taking it to the linear scale 
    random = XORshift()
    logX = random*4.698970004336019-4.
    x = 10.**logX
    #generating a probability for the y, f(x).
    y = XORshift()
    #if f(x) is smaller than the density function value for that x, keep the x, else, forget it.
    if float(y) <= float(p(x)):
        pointsList[n] = x
        n = n+1

print(pointsList)
#print(pointsList[-1])


#############################################################################################################################
###making 1000 halos each with 100 galaxies, which is the same as repeating the above code but with different random numbers.


#defining bins
binSize = 4.698970004336019/20.
binStarts = numpy.zeros(21)
logBinStarts = numpy.zeros(21)
for i in range(21):
    logBinStarts[i] = -4. + binSize*float(i)
    binStarts[i] = 10.**(logBinStarts[i])
binCounts = numpy.zeros(20,dtype=int)
biggestBinCounts = numpy.zeros(1000)
biggestBinX = []

halos = numpy.zeros((1000,100))
for haloNumber in range(1000):
    n = 0
    #xRandoms = XORshift(10*(haloNumber+1),1500)
    #yRandoms = XORshift(11*(haloNumber+1),1500)
    while n<100:
        logX = XORshift()*4.698970004336019-4.
        x = 10.**logX
        y = XORshift()
        if float(y) <= float(p(x)):
            halos[haloNumber,n] = x
            #checking which bin does this x belong to
            binNumber = int((logX + 4.)/binSize)
            binCounts[binNumber] = binCounts[binNumber]+1
            n = n+1
            #if it belongs to the biggest bin, save the x for that bin to the biggest binX and 
            if binNumber == 17:
                biggestBinCounts[haloNumber] = biggestBinCounts[haloNumber]+1
                biggestBinX.append(x)


###########
###plotting


xvals = numpy.arange(0.0001, 5, 0.0001)
plt.loglog(xvals,100.*p(xvals))
plt.ylabel("p(x)")
plt.xlabel("x")
plt.xlim(.0001,5)

binCountsAverage = numpy.zeros(20)
for i in range(20):
    binCountsAverage[i] = binCounts[i]/1000.

for i in range(20):
    xvals = numpy.arange(binStarts[i],binStarts[i+1],(binStarts[i+1]-binStarts[i])/3.)
    pltFunction = numpy.zeros(len(xvals))
    for j in range(len(xvals)):
        pltFunction[j] = binCountsAverage[i]
    plt.loglog(xvals,pltFunction)

#plt.show()
plt.savefig('binning.pdf')


##########################################################################################################
###this function will sort a given array in ascending order, for this I am using the quich sort algorithm.


#1.get the first, end, and middle element, and sort them. choose middle element as the pivot.
#2.for i and j in the same loop (i from 0 to n, j from n-1 to 0), if a[i]>=pivot, and a[j]<=pivot, if j>i, swap a[i], a[j]. aim of this step is to make sure that all the elements at pivots right are bigger(or equal) than pivot and all the elements at its left are smaller(or equal) than pivot.
#3.pivot is at its right place, so now take the left and rigth arrays excluding pivot and perform the same thing.
#ignore all the prints, they are just for myself in order to debug the code, if necessary.

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

#these are just tests, ignore them.
#myArray = [1012,57,42,63,97,1234,53,41253,112,4,566,123,34,153]
#myArray = [31,42,42,42,42,42,53,41253,112,4,566,123,34,153]
#pivotSort(myArray,0,13)


#######################################################################################################################
