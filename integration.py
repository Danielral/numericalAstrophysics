import numpy 
from randomNumberGenerator import XORshift

randoms = XORshift(10,3)
print(randoms)

a = randoms[0]/(2**64)*1.4 + 1.1
b = randoms[1]/(2**64)*1.5 + 0.5
c = randoms[2]/(2**64)*2.5 + 1.5
#a = 2
#b = 1
#c = 1

print(a,b,c)

def iWantToIntegrateThis(x):
    if x == 0: y = 0
    else:
        y = 4. * numpy.pi * (x/b)**(a-3.) * numpy.exp(-(x/b)**c) * x**2.
    return y

def Integration(f,a,b,steps):
    process = numpy.zeros((steps,steps),dtype=numpy.float64)
    process[0,0] = ( f(a)+f(b) )*(b-a)/2.
    
    for i in range(steps-1):
        process[i+1,0] = (b-a)*sum( f(a + (b-a)*(2.*j+1)/(2.**(i+1)) ) for j in range(2**i) )
        process[i+1,0] = process[i,0]/2. + process[i+1,0]/(2.**(i+1))
        for k in range(i+1):
            process[i+1,k+1] = process[i+1,k] + (process[i+1,k]-process[i,k])/(4**(k+1) -1.)

    return process

calc = Integration(iWantToIntegrateThis,0.,5.,6)
print(calc)

A = 1. / calc[-1,-1]
print(a,b,c,A)

def sqr(x):
    return x**2

sqrd = Integration(sqr,0.,5.,6)
print(sqrd)


