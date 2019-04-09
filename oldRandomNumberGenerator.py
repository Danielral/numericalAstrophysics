import numpy
import matplotlib.pylab as plot

##########################
###random number generator


#this function is going to generate an overflow warning, which is due to the overflowing in the MLCG part of the code, which is totally fine, beacause that overflowing is the same as multiplying in mode 2**64.

def XORshift(x,howMany):
    a = 21
    b = 35
    c = 4
    xors = numpy.zeros(howMany)
    #first using XOR Shift method
    for i in range(howMany):
        x ^= (x>>a)
        x ^= (x<<b)
        x ^= (x>>c)
        xors[i] = x
    #combining it with MLCG
    for i in range(howMany):
        #print(type(xors[i]))
        xors[i] = (7664345821815920749*xors[i])%(2**64)
    return xors


def randomGenerator(seed,count):
    a = 21
    b = 35
    c = 4
    rendomNumbers = numpy.zeros(count)
    seed


####################################################################################################
###random number generator, which for each time calling it, will produce a new random number for you


lastRandom = 5761015743107 #starts with this seed, but this is actually my storage unit for later stages               
print('This is my Seed: %i'%lastRandom)
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


####################################################################################################
###producing 1000 random numbers between 0 and 1, and then plotting the scatter plot of x_i+1 vs x_i


1000randomNumbers = numpy.empty(1000)
for i in range(1000):
    1000randomNumbers[i] = XORshift()

1001randomNumbers = numpy.empty(1000)
1001randomNumbers[1000] = 1000randomNumbers[0]
for i in range(999):
    1001randomNumbers[i] = 1000randomNumbers[i+1]

plot.scatter(1001randomNumbers,1000randomNumbers)


#############################################################################################
###generating one million random numbers between 0 and 1, and plotting the results in 20 bins


millionRandomNumbers = numpy.empty(1000000)
for i in range(1000000):
    millionRandomNumbers[i] = XORshift()

arange = numpy.arange(0,1.05,0.05)
plot.hist(millionRandomNumbers,arange)


######################################
