import random
n = 200 # n is the number of runs to simulate
ts = 10000 # ts is the number of time steps in a second
MtrStepDecider = random.randrange (ts) #this will be used to determine if try to step
MtrDetachDecider = random.randrange(ts) #This is used to see if falls off
# Numsteps = 0  #this initializes the counter for the number of steps the motor takes
for i in range(n):  #This is a loop which will go n times
    numtimesteps = 0
    Numsteps = 0
    MtrDetachDecider = 1000  #this makes sure you get into the while loop
    while MtrDetachDecider >= 100:
        numtimesteps = numtimesteps + 1
        MtrStepDecider = random.randrange (ts)
        # print('timestep is', numtimesteps, 'MtrStepDecider is', MtrStepDecider, 'MtrDetachDec is', MtrDetachDecider )
        if MtrStepDecider < 100:  #if true, then is attempting to step
            MtrDetachDecider = random.randrange(ts)
            if MtrDetachDecider >= 100:
                Numsteps = Numsteps + 1
              #  print ('entered if loop')
              #     print ('steps taken so far', Numsteps)
            else:
              #  print('got to else')
                print(Numsteps*8, numtimesteps, (Numsteps*8)/(numtimesteps/10000))
                #print(Numsteps)
                
	   
                
