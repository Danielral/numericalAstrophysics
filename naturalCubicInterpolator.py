import numpy


##################################################################################################
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
    for i in range(len(x)):
        xCopy.append(x[i])
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


########################################################################################################################
###now I want to produce that 6240 values for A(a,b,c)


#all a,b,c values in ascending order
aArray = numpy.arange(1.1,2.6,0.1)
bArray = numpy.arange(0.5,2.1,0.1)
cArray = numpy.arange(1.5,4.1,0.1)

#count = 0
Alist = numpy.zeros(((15),(16),(26)))
for i in range(15):
    for j in range(16):
        for k in range(26):
            def iWantToIntegrateThis(x):
                if x == 0: y = 0
                else:
                     y = 4. * numpy.pi * (x/bArray[j])**(aArray[i]-3.) * numpy.exp(-(x/bArray[j])**cArray[k]) * x**2.
                return y
            Alist[i][j][k] = findA()
            print('A(a,b,c) = %2.7f(%1.1f,%1.1f,%1.1f)'%(Alist[i][j][k],aArray[i],bArray[j],cArray[k]))
            #count = count + 1
            #print(count)        


#########################################
