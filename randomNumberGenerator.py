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


tousandRandomNumbers = numpy.empty(1000)
for i in range(1000):
    tousandRandomNumbers[i] = XORshift()

nextRandomNumbers = numpy.empty(1000)
nextRandomNumbers[999] = nextRandomNumbers[0]
for i in range(999):
    nextRandomNumbers[i] = tousandRandomNumbers[i+1]

plot.scatter(nextRandomNumbers,tousandRandomNumbers)


#############################################################################################
###generating one million random numbers between 0 and 1, and plotting the results in 20 bins


millionRandomNumbers = numpy.empty(1000000)
for i in range(1000000):
    millionRandomNumbers[i] = XORshift()

arange = numpy.arange(0,1.05,0.05)
plot.hist(millionRandomNumbers,arange)



