import numpy

def function(x):
    y = x**2. - 16.
    return y


#########################
###root finding algorithm


#I'm using the bisection method, beacasue it will converge for sure. this method works by linearly connecting two points of the function that have values with different sign. taking the avarage of these two points and calculating f at that point, if its bigger than zero ,then this is the new second guess for the root, and if its smaller than zero, this is the new first guess for the root. repeat this process untill you get close enough to the root.

def bisection(f,firstGuess,secondGuess,desiredError):
     f0 = f(firstGuess)
     f1 = f(secondGuess)
     newGuess = (firstGuess+secondGuess)/2.
     #print(firstGuess,secondGuess,desiredError,newGuess)
     fN = f(newGuess)
     if float(f0*fN) >= 0.:
        firstGuess = newGuess
     if float(f0*fN) < 0.:
        secondGuess = newGuess
     return newGuess if float(fN**2.) < float((desiredError)**2.) else bisection(f,firstGuess,secondGuess,desiredError)


#############################################
#firstTest = bisection(function,0.,8.,0.001)
#print(firstTest)


