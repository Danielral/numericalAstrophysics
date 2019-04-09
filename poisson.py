######
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


print('P(1,0) = %1.15f'%(poisson(1,0)))
print('P(5,10) = %1.15f'%(poisson(5,10)))
print('P(3,21) = %1.15f'%(poisson(3,21)))
print('P(2.6,4) = %1.15f'%(poisson(2.6,4)))
#print('P(
######
###
