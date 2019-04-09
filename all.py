import numpy
import matplotlib.pylab as plot

######1a
###two functions which together will return the poisson distribution at any give point:


#1.defining the poisson function when k = 0 (no factorial):
def zeroPoisson(lamb):
    result = numpy.float64(numpy.exp(-lamb))
    return result

#2.now that we have the poisson function value at a given lambda (with k = 0); one can derive a relation for p(lambda,k): p(lambda,k) = p(lambda,k-1) * lambda / k. this function will use this relation to find poisson's function value at the give k recursively.
def poisson(lamb,k):
    if k == 0:
        result = numpy.float64(zeroPoisson(lamb))
        return result
    result = numpy.float64(poisson(lamb,k-1) * lamb / k)
    return result

printLines = []
lastRandom = 5761015743107
printLines.append('This is my Seed: %i\n\n1(a):\n'%lastRandom)
printLines.append('P(1,0) = %1.30f\n'%(poisson(1,0)))
printLines.append('P(5,10) = %1.30f\n'%(poisson(5,10)))
printLines.append('P(3,21) = %1.30f\n'%(poisson(3,21)))
printLines.append('P(2.6,40) = %1.30f\n'%(poisson(2.6,40)))
printLines.append('P(101,200) = %1.30f\n'%(poisson(101,200)))
f = open('text/1a.txt','w')
f.writelines(printLines)
f.close()


######1b
###random number generator, which for each time calling it, will produce a new random number for you


printLines = []
lastRandom = 5761015743107 #starts with this seed, but this is actually my storage unit for later stages               
#printLines.append('This is my Seed: %i\n'%lastRandom)
def XORshift():
    a = 21
    b = 43 #XORshift values
    c = 4
    #first inserting the global randomNumber, which is the previous number
    global lastRandom
    #now use MLCG
    lastRandom = (2685821657736338717*lastRandom)%(2**64)
    #combine it with XORshift
    lastRandom ^= (lastRandom>>a)
    lastRandom ^= (lastRandom<<b)
    lastRandom ^= (lastRandom>>c)
    #to specify a range for my random number
    return float(lastRandom%(2**64))/(2**64)


###producing 1000 random numbers between 0 and 1, and then plotting the scatter plot of x_i+1 vs x_i


tousandRandomNumbers = numpy.empty(1000)
for i in range(1000):
    tousandRandomNumbers[i] = XORshift()

nextRandomNumbers = numpy.empty(1000)
nextRandomNumbers[999] = nextRandomNumbers[0]
for i in range(999):
    nextRandomNumbers[i] = tousandRandomNumbers[i+1]

plot.scatter(nextRandomNumbers,tousandRandomNumbers)
plot.title('Scatter Plot of Sequential Random Numbers')
plot.savefig('plots/1bScatter.png')
plot.close()


###generating one million random numbers between 0 and 1, and plotting the results in 20 bins


millionRandomNumbers = numpy.empty(1000000)
for i in range(1000000):
    millionRandomNumbers[i] = XORshift()

arange = numpy.arange(0,1.05,0.05)
plot.hist(millionRandomNumbers,arange)
plot.xlabel('20 bins')
plot.ylabel('counts')
plot.title('20 bins of 1,000,000 random numbers')
plot.savefig('plots/1bHistogram.png')
plot.close()


######2a
###producing random values for a,b,c


a = XORshift()*1.4 + 1.1
b = XORshift()*1.5 + 0.5
c = XORshift()*2.5 + 1.5

###defining the function that I want to integrate, in order to find a.


def iWantToIntegrateThis(x):
    global a
    global b
    global c
    if x == 0: y = 0
    else:
        y = 4. * numpy.pi * (x/b)**(a-3.) * numpy.exp(-(x/b)**c) * x**2.
    return y


###building up the numerical integrator, using Romberg algorithm


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


###now using the defined function and the numerical integrator, I can integrate that function, hence the value of A(a,b,c).


def findA():
    denom = integration(iWantToIntegrateThis,0.,5.,5)
    print(denom)
    return 1./denom[-1,-1]

A = findA()
printLines = []
printLines.append('A(a,b,c) = %f, a = %f, b = %f, c = %f\n'%(A,a,b,c))
f = open('text/2a.txt','w')
f.writelines(printLines)
f.close()


######2b
###defining the function again N = 100


def n(x):
    global a
    global b
    global c
    global A
    y = A * 100 * (x/b)**(a-3.) * numpy.exp(-(x/b)**c)
    return y


###in the upcoming lines, I will define 4 functions, which I need for Natural Cubic Interpolation:


#(first) solve the tridiagonal linear equation (which is the result of Cubic Interpolation method), and by solving I mean: given the 3 diagonal arrays (which I will define later in the defineCubicSpline function), this function will solve the set of linear equations by matrice manipulation methods, and returns the solutions, which are y second derivatives at the given points.

#up is the upper diagonal, diag is the middle, and down in the bottom diagonl, rhs (short for rigth hand side), is the right hand side of each equation. I want to achieve a diagonal matrice with the diagonal values all equal to one (basic linear algebra), which in this case (tridiagonal matrice) is achievable, at the end the rhs values, would be the solutions:

#1.first lets make a two diagonal matrice by getting rid of the bottom diag. so I solve for each down[i] and add the upper row to the current row, to make the down[i] zero, since for each row, down[i] is below the diag[i-1] of the previous row, I need to add -row[i-1]*down[i-1]/diag[i-1] to each row which will result down[i] becoming zero, and since I know that already, Im just keeping track of the rest of the row (up,diag,and rhs).
#2.after the first loop is finished, since the last row contains only diag and rhs, its solvable for rhs, meaning dividing the row by diag value to make diag equal to one (remember, from algebra, this is the goal)
#3.now that we know the solution for the last row, we can use it to find the rest of the solutions, with the same method as above, but this time moving from bottom to top instead, and getting rid of up elements, same as we did for the down elements.
#4.at the end I'm left out with a diagonal matrice with all diag values equal to one, so the rhs of this matrice is simply the solution to the starting linear equation.

def solve(up,diag,down,rhs):
    #1
    n = len(diag)
    for i in range(n-1):
        diag[i+1] = diag[i+1] - up[i]*down[i]/diag[i]
        rhs[i+1] = rhs[i+1] - rhs[i]*down[i]/diag[i]
    #2
    rhs[-1] = rhs[-1]/diag[-1]
    diag[-1] = 1.
    #3
    for i in range(n-1):
        rhs[n-2-i] = (rhs[n-2-i] - up[n-2-i]*rhs[n-1-i])/diag[n-2-i]
        diag[n-2-i] = 1.
        up[n-2-i] = 0.
    #4
    return rhs


#(second) defining the Cubic Spline, from the given data set, simply by the given formula for the cubic spline.

#1.the loop simply goes through all the data set, and calculates the coefficients of the cubic interpolation formula, which together will define a tridiagonal matrice,solvable to find the values of y second derivatives at each point.
#2.at the end I'm returning this coefficients as the three diagonal arrays of the tridiagonal matrice, which will be solved by the above function.

def defineCubicSpline(x,y):
    n = len(x)
    up = numpy.zeros(n-1,dtype=numpy.float64)
    diag = numpy.ones(n,dtype=numpy.float64)
    down = numpy.zeros(n-1,dtype=numpy.float64)
    rhs = numpy.zeros(n,dtype=numpy.float64)
    #1
    for i in range(1,n-1):
        diag[i] = (x[i+1]-x[i-1])/3.
        rhs[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        if i < n-2:
            down[i] = (x[i]-x[i-1])/6.
            up[i] = (x[i+1]-x[i])/6.
    #2
    return up,diag,down,rhs

#(third) finds where does x belong in a give data set (segment).

#tries one point in the middle of the array, if x is bigger than this point, gets rid of the bottom half, and if x is smaller than that point, gets rid of the upper half. do this till left out with an array of length 2. returns that array, and the original array number of the bottom member.

def findX(array,x):
    i = 0
    seq = array
    while len(seq) != 2:
        m = int(len(seq)/2)
        if x >= seq[m]:
            del seq[0:m]
            i = i + m
        elif x < seq[m]:
            del seq[m+1:]
    return seq,i


#(fourth) now that we know where does this point of interest belongs in our grid (by using the function above), and finnaly, we can insert y second derivative values(solved by functions 1 and 2) into the natural cubic interpolator formula, to get back the interpolated value at the point of interest, f(point of interest).


def cubicInterpolation(pointOfInterest,segmentNumber,cubicSolutions,x,y):
    i = segmentNumber
    result = (cubicSolutions[i]/6.)*( (pointOfInterest-x[i+1])**3./(x[i]-x[i+1]) - (pointOfInterest-x[i+1])*(x[i]-x[i+1]) )
    result = result - (cubicSolutions[i+1]/6.)*( (pointOfInterest-x[i])**3./(x[i]-x[i+1])-(pointOfInterest-x[i])*(x[i]-x[i+1]) )
    result = result + ( y[i]*(pointOfInterest-x[i+1])-y[i+1]*(pointOfInterest-x[i]) )/(x[i]-x[i+1])
    return result


#putting four functions together to make the one-dimensional interpolator, which includes:
#1.defineCubicSpline
#2.solve
#3.findX
#4.cubicInterpolator
#this function gets two arrays for x(should be sorted), and f(x) values; and a point of interest, and retturns the interpolated values for f(point of interest) at the end.

def interpolate(xArray,yArray,pointOfInterest):
    #find the tridiagonal matrice from the given values for xArray, and yArray:
    up,diag,down,rhs = defineCubicSpline(xArray,yArray)
    #solve the tridiagonal matrice to find the values for y's second derivative at each one of the given x points:
    solution = solve(up,diag,down,rhs)
    #copy xArray into a new array(xCopy) to keep the xArray safe from the manipulation of the upcoming function(findX):
    xCopy = []
    for i in range(len(xArray)):
        xCopy.append(xArray[i])
    #finding the location of the point of interst in the sorted x array, which returns the closest upper and lower x to the point of interest as the segment, and the array number of the lower x belonging to that segment as lowerX(so, for example, if point of interest lies between the first and second members of the xArray, this function will return [xArray[0],xArray[1]] as the segment, and 0 as lowerX.
    segment,lowerX = findX(xCopy,pointOfInterest)
    #now we know where does this point of interest belongs in our grid, and finnaly, we can insert y second derivative values into the natural cubic interpolator formula, to get back the interpolated value at the point of interest, f(point of interest):
    valueAtPointOfInterest = cubicInterpolation(pointOfInterest,lowerX,solution,xArray,yArray)
    return valueAtPointOfInterest


################################################################################################################################
#now I want to use the above function to make a three dimensional interpolator, but for upper dimensions, the idea would remain the same:


#first, what I want to do? well, I know the values of f at a set of points (xi,yi,zi), and I want to find its value at (x,y,z). for that, I will devide the interpolation to three one-dim interpolations, and apply them as follows, and in each step I attemp to lower the dimension of the function by one, by getting rid of one of the coordinates.
#1.imagine the line f(yConstant,zConstant). for this line, I can interpolate the f(x,yConstant,zConstant) which is a 1-dim function, f_yConstant_zConstant(x)
#2.repeat the above process, but for a different value of yConstant this time. this is like doing the interpolation in the surface of f(zConstant).
#3.do the step 2, till you go through all the yi. at this point you have a set of new grid points which would look like this: f(x,yi,zConstant). which define a line passing through f(x,y,zConstant).
#4.interpolate on this line, and you will get the value at f(x,y,zConstant).
#5.to the above steps for all the zi, and you will end up with a line f(x,y,zi) passong through f(x,y,z).
#6.finally interpolate that last line to get the f(x,y,z).

def threeDInterpolate(xArray,yArray,zArray,fArray,xInterest,yInterest,zInterest):
    zConstantNumber = 0
    zLine = numpy.zeros(len(zArray)) #defining the line of f(xInterest,yInterest,zArray)
    #loop to produce the f(xInterest,yInterest,zArray) for all the zArray
    for zConstant in zArray:
        yConstantNumber = 0
        yLine = numpy.zeros(len(yArray)) #defining the line of f(xInterest,yArray,zConstant)
        #loop to produce the f(xInterest,yArray,zConstant) for all the yArray
        for yConstant in yArray:
            xLine = numpy.zeros(len(xArray)) #defining the line of f(xArray,yConstant,zConstant)
            #loop to produce the f(xArray,yConstant,zConstant) for all xArray
            for i in range(len(xArray)):
                xLine[i] = fArray[i,yConstantNumber,zConstantNumber]
            yLine[yConstantNumber] = interpolate(xArray,xLine,xInterest)
            yConstantNumber = yConstantNumber + 1
        zLine[zConstantNumber] = interpolate(yArray,yLine,yInterest)
        zConstantNumber = zConstantNumber + 1
    result = interpolate(zArray,zLine,zInterest)
    return result


x = [10.**(-4),10.**(-2),10.**(-1),1.,5.]
y = []

for i in range(len(x)):
    y.append(n(x[i]))

plot.scatter(x,y)
plot.xlabel('x')
plot.ylabel('n(x)')
plot.xlim(10**(-4),5)
plot.xscale('log')
plot.yscale('log')
plot.title('interpolation plot')
#plot.show()


xPlot = numpy.arange(-4,0.67,0.01)
xFinal = []
for i in range(len(xPlot)):
    xFinal[i] = 10**(xPlot[i])
yPlot = 0
yFinal = numpy.zeros(len(xPlot))
for i in range(len(xPlot)):
    yPlot = interpolate(x,y,xFinal[i])
    yFinal[i] = yPlot
plot.loglog(xFinal,yFinal)
######

