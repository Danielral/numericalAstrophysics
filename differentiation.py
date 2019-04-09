import numpy 

def n(x,a,b,c,A,N):
    y = A*N*(x/b)**(a-3) * numpy.exp(-(x/b)**c)
    return y


###########################################################
###Ridder's differentiation method, returns the whole table


#1.caculate the centeral differential
#2.decreasing h and repeating 1
#3.combining pairs using Ridder's formula

def differentiator(f,x,h,d,steps,analAnswer,howManyDigits):
    #1,2
    process = numpy.zeros((steps,steps),dtype = numpy.float64)
    #calculating d(x), for n values of h, in descending oreder
    for n in range(steps):
        process[n,0] = ( f(x+h/d**n)-f(x-h/d**n) )/( 2.*(h/d**n) )
    bestError = float((process[steps-1,0]-analAnswer)**2.)  #choosing the best error
    #3
    for j in range(1,steps):
        for i in range(0,steps-j):
            process[i,j] = ( (d**(2.*j))*process[i+1,j-1]-process[i,j-1] )/((d**(2.*j))-1.)
            error = float((process[i,j]-analAnswer)**2.)
            #if the error is small enough, break it 
            if error < 10.**(-(2.*howManyDigits)):
                return process[i,j],error
                break
        #if the error is getting large, break and output the table
        if float(bestError) < float(error):
            return process[i,j],error
            break
        bestError = error 
    return process,bestError


#############################################################################

def function(x):
    y = numpy.exp(x)/(numpy.sin(x)-x**2.)
    return y

def cube(x):
    y = x**3
    return y
