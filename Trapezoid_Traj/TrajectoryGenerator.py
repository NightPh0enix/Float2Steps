
# pr 1 ts=2481961511 DDA: start=[-6.000000 0.000000 0.000000] end=[25.000000 0.000000 0.000000]
# ***Distances***
#  accelStopDistance = 2.500000 steadyDistance = 28.500000 decelStartDistance = 31.000000
# ***Times***
#  accelStopTime = 93750.000000 steadyTime = 534375.000000 decelStartTime = 628125.000000 totalTime = 628125
#  vec=[1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000]
# a=5.6889e-10 d=5.6889e-10 reqv=5.3333e-5 startv=0.0000e+0 topv=5.3333e-5 endv=5.3333e-5 cks=628125 fp=4294967295 fl=0541
# DMX: dir=F steps=2480 next=1 rev=2481 interval=6629 psl=201 A=0.0000e+0 B=0.0000e+0 C=4.3945e+7 dsf=2.5000e+0 tsf=93750.0
# Error: Attempt to move motors when VIN is not in range
# pr 1 ts=2482589636 DDA: start=[25.000000 0.000000 0.000000] end=[40.000000 0.000000 0.000000]
# ***Distances***
#  accelStopDistance = -0.000000 steadyDistance = 12.500000 decelStartDistance = 12.500000
# ***Times***
#  accelStopTime = 0.000000 steadyTime = 234375.000000 decelStartTime = 234375.000000 totalTime = 328125
#  vec=[1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000]
# a=5.6889e-10 d=5.6889e-10 reqv=5.3333e-5 startv=5.3333e-5 topv=5.3333e-5 endv=0.0000e+0 cks=328125 fp=4294967295 fl=0545
# DMX: dir=F steps=1200 next=1 rev=1201 interval=234 psl=1001 A=1.7776e+10 B=0.0000e+0 C=2.3438e+2 dsf=1.2500e+1 tsf=234375.0



import matplotlib.pyplot as plt
import numpy as np 
import math

startTime=0
rounding = 0

clockRate=937500
clocksToSecMultiplier=1/clockRate
acceleration = 5.6889e-10 #mm/cks^2

###################################    Parameters subject to change    ###################################
startVelocity=0
startDist=0
time_interval=0.01

frac, whole = math.modf(time_interval)
frac=str(frac)
for j in frac :
    rounding=rounding+1
    if (j=='.'):
        rounding = 0
###################################    dist in mm     ################################### 

# accelStopDistance = 2.500000 
# steadyDistance = 0.056250 
# decelStartDistance = 2.556250

accelStopDistance=2.500000
steadyDistance=5.000000 # distnce covered while we are moving with constant speed
decelStartDistance=7.500000
totalDistance=10

# delta_time_acc=accelStopTime-startTime
# delta_time_steady=steadyTime
# delta_time_decel=totalTime-decelStartTime
# delta_dist_acc=accelStopDistance-startDist
# delta_dist_steady=steadyDistance
# delta_dist_decel=totalDistance-decelStartDistance

###################################    Times in clocks    ###################################  
# accelStopTime = 93750.000000 
# steadyTime = 1054.684814 
# decelStartTime = 94804.687500 # steadyTime + accelStopTime
# totalTime = 174492




accelStopTime=93750
steadyTime=93750 #for how long are we moving with const speed
decelStartTime=187500 # steadyTime + accelStopTime
totalTime=281250


###################################    Times in seconds    ###################################

totalTime=round(totalTime*clocksToSecMultiplier, rounding)
# print(totalTime)

accelStopTime=round(accelStopTime*clocksToSecMultiplier, rounding)
# print(accelStopTime)

steadyTime=round(steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
# print(steadyTime)

decelStartTime=round(decelStartTime*clocksToSecMultiplier, rounding)
# print(decelStartTime)

acceleration = round(acceleration*clockRate*clockRate,rounding)
deceleration =acceleration*(-1)

startVelocity=round(startVelocity*clockRate,rounding)

t1=accelStopTime
t2=decelStartTime
t3=totalTime

T1=accelStopTime
T2=steadyTime=decelStartTime-accelStopTime
###################################    setting up time list    ################################### 
loop_iterator=totalTime+time_interval
time=[] 
for i in np.arange(0,loop_iterator,time_interval):
    # time.append(round(i,rounding)) # increments in time will be in 0.1 sec step i.e. 0.1 0.2 0.3
    time.append(round(i,rounding))

print(time)

print("Time len",len(time))
###################################    velocity profile    ###################################

f_2e=F=f_1e=startVelocity+acceleration*T1

velocity=[]
for t in time:
    if (t>=0 and t<=t1):
        vel = startVelocity + acceleration * t
        velocity.append(vel)

    elif (t>=t1 and t<=t2):
        vel=F
        velocity.append(vel)

    elif (t>=t2 and t<=t3):
        vel = f_2e + deceleration * (t-t2)
        velocity.append(vel)
print("vello",velocity)
###################################    distance profile    ###################################

s_1e = startDist + startVelocity*T1 + 0.5*acceleration*T1**2
s_2e = s_1e + F*T2
# print("t1 t2 t3",t1,t2,t3)
# print("F",F)
# print(s_1e,s_2e)

distance=[]
for t in time:
    if (t>=0 and t<=t1):
        dist = startDist+startVelocity*t + 0.5*acceleration*t**2
        distance.append(dist)
    elif (t>=t1 and t<=t2):
        dist = s_1e+F*(t-t1)
        distance.append(dist)
    elif (t>=t2 and t<=t3):
        dist = s_2e +f_2e*(t-t2) + 0.5*deceleration*(t-t2)**2
        distance.append(dist)

print(distance)
print("hehe:",distance[-1])
# print(distance[-1],time[-1],velocity[-1])
# ###################################   Plotting Graphs    ###################################

# print("totalTime",totalTime)
# print("totalTime_move1",totalTime_move1)
# print("totalTime_move2",totalTime_move2)
# print("accelStopTime",accelStopTime)
# print("steadyTime",steadyTime)
# print("decelStartTime",decelStartTime)
# print("startVelocity",startVelocity)
# print("t1",t1)
print("t2",t2)
# print("t3",t3)
# print("T1",T1)
# print("T2",T2)
# print("T3",T3)
# print("F, f_1e, f_2e", F)

fig, axs = plt.subplots(2)
fig.suptitle('Motion Profiles')
axs[0].plot(time, velocity)
axs[0].set_title("Velocity Profile")
axs[0].set_xlabel("time (s)")
axs[0].set_ylabel("velocity (mm/s)")

axs[1].plot(time, distance)
axs[1].set_title("Distance Profile")
axs[1].set_xlabel("time (s)")
axs[1].set_ylabel("Distance (mm)")
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()


######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################

######################################################################################### FOR 2 MOVES ############################################################################################################
######################################################################################### FOR 2 MOVES  ############################################################################################################


# import matplotlib.pyplot as plt
# import numpy as np 
# import math

# startTime=0
# rounding = 0

# clockRate=937500
# clocksToSecMultiplier=1/clockRate
# acceleration = 5.6889e-10 #mm/cks^2

# ###################################    Parameters subject to change    ###################################
# startVelocity=0
# startDist=0
# time_interval=0.1

# frac, whole = math.modf(time_interval)
# frac=str(frac)
# for j in frac :
#     rounding=rounding+1
#     if (j=='.'):
#         rounding = 0
# ###################################    dist in mm     ################################### 
# accelStopDistance = 2.500000 
# steadyDistance = 45.056252 
# decelStartDistance = 47.556252

# # accelStopDistance = 2.500000 
# # steadyDistance = 28.500000 
# # decelStartDistance = 31.000000

# ###################################    Times in clocks    ###################################  
# # accelStopTime = 93750.000000 
# # steadyTime = 534375.000000 
# # decelStartTime = 628125.000000 

# accelStopTime = 93750.000000 
# steadyTime = 844804.687500 
# decelStartTime = 938554.687500 


# totalTime_move1 = 1018242
# totalTime_move2 = 1205742

# totalTime = totalTime_move1+totalTime_move2
# ###################################    Times in seconds    ###################################

# totalTime=round(totalTime*clocksToSecMultiplier, rounding)

# totalTime_move1=round(totalTime_move1*clocksToSecMultiplier, rounding)
# totalTime_move2=round(totalTime_move2*clocksToSecMultiplier, rounding)

# accelStopTime=round(accelStopTime*clocksToSecMultiplier, rounding)
# steadyTime=round(steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
# decelStartTime=round(decelStartTime*clocksToSecMultiplier, rounding)


# acceleration = round(acceleration*clockRate*clockRate,rounding)
# deceleration =acceleration*(-1)

# startVelocity=round(startVelocity*clockRate,rounding)

# t1=accelStopTime
# t2=decelStartTime
# t3=totalTime_move1

# T1=accelStopTime
# T2=steadyTime=decelStartTime-accelStopTime
# T3=totalTime_move1-decelStartTime


# ###################################    setting up time list    ################################### 
# loop_iterator=totalTime+time_interval
# # loop_iterator=totalTime_move1+time_interval
# time=[] 
# for i in np.arange(0,loop_iterator,time_interval):
#     # time.append(round(i,rounding)) # increments in time will be in 0.1 sec step i.e. 0.1 0.2 0.3
#     time.append(round(i,rounding))

# ###################################    velocity profile    ###################################

# f_2e=F=f_1e=startVelocity+acceleration*T1


# velocity=[]
# loopcounter=0
# for t in time:
#     loopcounter=loopcounter+1
#     if (t<=totalTime_move1):
#         # print("T value is ",t)
#         if (t>=0 and t<=t1):
#             # velocity.append("AA")
#             vel = startVelocity + acceleration * t
#             velocity.append(vel)
#             # print("A",t)
#         elif (t>=t1 and t<=t2):
#             # velocity.append("BB")
#             vel=F
#             velocity.append(vel)
#             # print("B",t)
#         elif (t>=t2 and t<=t3):
#             # velocity.append("CC")
#             vel = f_2e + deceleration * (t-t2)
#             velocity.append(vel)
#             # print("C",t)


# print(velocity)


# # ###################################    distance profile    ###################################

# s_1e = startDist + startVelocity*T1 + 0.5*acceleration*T1**2
# s_2e = s_1e + F*T2
# print("t1 t2 t3",t1,t2,t3)


# distance=[]
# for t in time:
#     if (t<=totalTime_move1):
#         if (t>=0 and t<=t1):
#             dist = startDist+startVelocity*t + 0.5*acceleration*t**2
#             distance.append(dist)
#         elif (t>=t1 and t<=t2):
#             dist = s_1e+F*(t-t1)
#             distance.append(dist)
#         elif (t>=t2 and t<=t3):
#             dist = s_2e +f_2e*(t-t2) + 0.5*deceleration*(t-t2)**2
#             distance.append(dist)


# print("totalTime",totalTime)
# print("totalTime_move1",totalTime_move1)
# print("totalTime_move2",totalTime_move2)
# print("accelStopTime",accelStopTime)
# print("steadyTime",steadyTime)
# print("decelStartTime",decelStartTime)
# print("startVelocity",startVelocity)
# print("t1",t1)
# print("t2",t2)
# print("t3",t3)
# print("T1",T1)
# print("T2",T2)
# print("T3",T3)
# print("F, f_1e, f_2e", F)


# # ###################################   MOVE 2   ###################################


# move2_startVelocity=velocity[-1]
# move2_startDistance=distance[-1]

# move2_accelStopTime=79687.500000
# move2_steadyTime=1032304.687500
# move2_decelStartTime=1111992.250000
# move2_totaltime=1205742


# move2_accelStopTime=round(move2_accelStopTime*clocksToSecMultiplier, rounding)
# move2_steadyTime=round(move2_steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
# move2_decelStartTime=round(move2_decelStartTime*clocksToSecMultiplier, rounding)
# move2_totaltime=round(move2_totaltime*clocksToSecMultiplier, rounding)

# move2_t1=move2_accelStopTime
# move2_t2=move2_decelStartTime
# move2_t3=move2_totaltime

# move2_T1=move2_accelStopTime
# move2_T2=move2_steadyTime

# move2_loop_iterator=totalTime_move2+time_interval

# move2_F=move2_f_2e=move2_f_1e=move2_startVelocity+acceleration*move2_T1
# move2_s_1e = move2_startDistance + move2_startVelocity*move2_T1 + 0.5*acceleration*move2_T1**2
# move2_s_2e = move2_s_1e + move2_F*move2_T2

# move2_timearray=[] 

# for k in np.arange(0,move2_loop_iterator,time_interval):
#     move2_timearray.append(round(k,rounding))


# for q in move2_timearray:
    
#     if (q>0 and q<=move2_t1):
#         # velocity.append("A")
#         vel = startVelocity + acceleration * q
#         velocity.append(vel)

#     elif (q>=move2_t1 and q<=move2_t2):
#         # velocity.append("B")
#         vel=move2_F
#         velocity.append(vel)

#     elif (q>=move2_t2 and q<=move2_t3):
#         # velocity.append("C")
#         vel = f_2e + deceleration * (q-move2_t2)
#         velocity.append(vel)





# for p in move2_timearray:

#     if (p>0 and p<=move2_t1):
#         # distance.append("A")
#         dist = move2_startDistance+move2_startVelocity*p + 0.5*acceleration*p**2
#         distance.append(dist)

#     elif (p>move2_t1 and p<=move2_t2):
#         # distance.append("B")
#         dist = move2_s_1e + move2_F * (p-move2_t1)
#         distance.append(dist)

#     elif (p>move2_t2 and p<=move2_t3):
#         # distance.append("C")
#         dist = move2_s_2e + move2_f_2e*(p-move2_t2) + 0.5*deceleration*(p-move2_t2)**2
#         distance.append(dist)



# # print("move2_accelStopTime",move2_accelStopTime)
# # print("move2_steadyTime",move2_steadyTime)
# # print("move2_decelStartTime",move2_decelStartTime)
# # print("move2_startVelocity",move2_startVelocity)
# # print("move2_t1",move2_t1)
# # print("move2_t2",move2_t2)
# # print("move2_t3",move2_t3)
# # print("move2_T1",move2_T1)
# # print("move2_T2",move2_T2)
# # print("move2_F",move2_F)
# # print("move2_s_1e,",move2_s_1e)
# # print("move2_s_2e",move2_s_2e)        








# # print(velocity)
# # print(len(time))
# # print(len(velocity))
# # print(len(distance))

# fig, axs = plt.subplots(2)
# fig.suptitle('Motion Profiles')

# axs[0].plot(time, velocity)
# axs[0].set_title("Velocity Profile")
# axs[0].set_xlabel("time (s)")
# axs[0].set_ylabel("velocity (mm/s)")
# axs[0].axvline(x=accelStopTime, color='r', linestyle='dashed')
# axs[0].axvline(x=decelStartTime, color='g', linestyle='dashed')
# axs[0].axvline(x=move2_accelStopTime+totalTime_move1, color='b', linestyle='dashed')
# axs[0].axvline(x=move2_decelStartTime+totalTime_move1, color='c', linestyle='dashed')

# axs[1].plot(time, distance)
# axs[1].set_title("Distance Profile")
# axs[1].set_xlabel("time (s)")
# axs[1].set_ylabel("Distance (mm)")
# axs[1].axvline(x=accelStopTime, color='r', linestyle='dashed')
# axs[1].axvline(x=decelStartTime, color='g', linestyle='dashed')
# axs[1].axvline(x=move2_accelStopTime+totalTime_move1, color='b', linestyle='dashed')
# axs[1].axvline(x=move2_decelStartTime+totalTime_move1, color='c', linestyle='dashed')

# fig.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.show()






######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################
######################################################################################### END ############################################################################################################

######################################################################################### FOR -ve +ve MOVES ############################################################################################################
######################################################################################### FOR -ve +ve MOVES  ############################################################################################################



# import matplotlib.pyplot as plt
# import numpy as np 
# import math

# startTime=0
# rounding = 0

# clockRate=937500
# clocksToSecMultiplier=1/clockRate
# acceleration = 5.6889e-10 #mm/cks^2

# ###################################    Parameters subject to change    ###################################
# startVelocity=0
# startDist=0
# time_interval=0.1

# frac, whole = math.modf(time_interval)
# frac=str(frac)
# for j in frac :
#     rounding=rounding+1
#     if (j=='.'):
#         rounding = 0
# ###################################    dist in mm     ################################### 
# accelStopDistance = 2.500000 
# steadyDistance = 45.056252 
# decelStartDistance = 47.556252

# # accelStopDistance = 2.500000 
# # steadyDistance = 28.500000 
# # decelStartDistance = 31.000000

# ###################################    Times in clocks    ###################################  
# # accelStopTime = 93750.000000 
# # steadyTime = 534375.000000 
# # decelStartTime = 628125.000000 

# accelStopTime = 93750.000000 
# steadyTime = 844804.687500 
# decelStartTime = 938554.687500 


# totalTime_move1 = 1018242
# totalTime_move2 = 1205742

# totalTime = totalTime_move1+totalTime_move2
# ###################################    Times in seconds    ###################################

# totalTime=round(totalTime*clocksToSecMultiplier, rounding)

# totalTime_move1=round(totalTime_move1*clocksToSecMultiplier, rounding)
# totalTime_move2=round(totalTime_move2*clocksToSecMultiplier, rounding)

# accelStopTime=round(accelStopTime*clocksToSecMultiplier, rounding)
# steadyTime=round(steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
# decelStartTime=round(decelStartTime*clocksToSecMultiplier, rounding)


# acceleration = round(acceleration*clockRate*clockRate,rounding)
# deceleration =acceleration*(-1)

# startVelocity=round(startVelocity*clockRate,rounding)

# t1=accelStopTime
# t2=decelStartTime
# t3=totalTime_move1

# T1=accelStopTime
# T2=steadyTime=decelStartTime-accelStopTime
# T3=totalTime_move1-decelStartTime


# ###################################    setting up time list    ################################### 
# loop_iterator=totalTime+time_interval
# # loop_iterator=totalTime_move1+time_interval
# time=[] 
# for i in np.arange(0,loop_iterator,time_interval):
#     # time.append(round(i,rounding)) # increments in time will be in 0.1 sec step i.e. 0.1 0.2 0.3
#     time.append(round(i,rounding))

# ###################################    velocity profile    ###################################

# f_2e=F=f_1e=startVelocity+acceleration*T1


# velocity=[]
# loopcounter=0
# for t in time:
#     loopcounter=loopcounter+1
#     if (t<=totalTime_move1):
#         # print("T value is ",t)
#         if (t>=0 and t<=t1):
#             # velocity.append("AA")
#             vel = startVelocity + acceleration * t
#             velocity.append(vel)
#             # print("A",t)
#         elif (t>=t1 and t<=t2):
#             # velocity.append("BB")
#             vel=F
#             velocity.append(vel)
#             # print("B",t)
#         elif (t>=t2 and t<=t3):
#             # velocity.append("CC")
#             vel = f_2e + deceleration * (t-t2)
#             velocity.append(vel)
#             # print("C",t)





# # ###################################    distance profile    ###################################

# s_1e = startDist + startVelocity*T1 + 0.5*acceleration*T1**2
# s_2e = s_1e + F*T2
# # print("t1 t2 t3",t1,t2,t3)


# distance=[]
# displacement=[]

# for t in time:
#     if (t<=totalTime_move1):
#         if (t>=0 and t<=t1):
#             dist = startDist+startVelocity*t + 0.5*acceleration*t**2
#             distance.append(dist)
#         elif (t>=t1 and t<=t2):
#             dist = s_1e+F*(t-t1)
#             distance.append(dist)
#         elif (t>=t2 and t<=t3):
#             dist = s_2e +f_2e*(t-t2) + 0.5*deceleration*(t-t2)**2
#             distance.append(dist)

# displacement=distance.copy()


# if displacement==distance:
#     print("NONONO")

# print("totalTime",totalTime)
# print("totalTime_move1",totalTime_move1)
# print("totalTime_move2",totalTime_move2)
# print("accelStopTime",accelStopTime)
# print("steadyTime",steadyTime)
# print("decelStartTime",decelStartTime)
# print("startVelocity",startVelocity)
# print("t1",t1)
# print("t2",t2)
# print("t3",t3)
# print("T1",T1)
# print("T2",T2)
# print("T3",T3)
# print("F, f_1e, f_2e", F)


# # ###################################   MOVE 2   ###################################


# move2_startVelocity=velocity[-1]
# move2_startDistance=distance[-1]





# move2_accelStopTime=79687.500000
# move2_steadyTime=1032304.687500
# move2_decelStartTime=1111992.250000
# move2_totaltime=1205742


# move2_accelStopTime=round(move2_accelStopTime*clocksToSecMultiplier, rounding)
# move2_steadyTime=round(move2_steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
# move2_decelStartTime=round(move2_decelStartTime*clocksToSecMultiplier, rounding)
# move2_totaltime=round(move2_totaltime*clocksToSecMultiplier, rounding)

# move2_t1=move2_accelStopTime
# move2_t2=move2_decelStartTime
# move2_t3=move2_totaltime

# move2_T1=move2_accelStopTime
# move2_T2=move2_steadyTime

# move2_loop_iterator=totalTime_move2+time_interval

# move2_F=move2_f_2e=move2_f_1e=move2_startVelocity+acceleration*move2_T1
# move2_s_1e = move2_startDistance + move2_startVelocity*move2_T1 + 0.5*acceleration*move2_T1**2
# move2_s_2e = move2_s_1e + move2_F*move2_T2

# move2_timearray=[] 

# for k in np.arange(0,move2_loop_iterator,time_interval):
#     move2_timearray.append(round(k,rounding))


# for q in move2_timearray:
    
#     if (q>0 and q<=move2_t1):
#         # velocity.append("A")
#         vel = startVelocity + acceleration * q
#         velocity.append(vel)

#     elif (q>=move2_t1 and q<=move2_t2):
#         # velocity.append("B")
#         vel=move2_F
#         velocity.append(vel)

#     elif (q>=move2_t2 and q<=move2_t3):
#         # velocity.append("C")
#         vel = f_2e + deceleration * (q-move2_t2)
#         velocity.append(vel)


# for p in move2_timearray:

#     if (p>0 and p<=move2_t1):
#         # distance.append("A")
#         dist = move2_startDistance+move2_startVelocity*p + 0.5*acceleration*p**2
#         distance.append(dist)

#     elif (p>move2_t1 and p<=move2_t2):
#         # distance.append("B")
#         dist = move2_s_1e + move2_F * (p-move2_t1)
#         distance.append(dist)

#     elif (p>move2_t2 and p<=move2_t3):
#         # distance.append("C")
#         dist = move2_s_2e + move2_f_2e*(p-move2_t2) + 0.5*deceleration*(p-move2_t2)**2
#         distance.append(dist)


# for z in move2_timearray:

#     if (z>0 and z<=move2_t1):
#         # displacement.append("A")
#         displ = move2_startDistance - move2_startVelocity*z - 0.5*acceleration*z**2
#         displacement.append(displ)
#         last_dist=displacement[-1]

#     elif (z>move2_t1 and z<=move2_t2):
#         # displacement.append("B")
#         displ = last_dist - move2_F * (z-move2_t1)
#         displacement.append(displ)
#         last_dist_steady=displacement[-1]

#     elif (z>move2_t2 and z<=move2_t3):
#         # displacement.append("C")
#         displ = last_dist_steady - move2_f_2e*(z-move2_t2) - 0.5*deceleration*(z-move2_t2)**2
#         displacement.append(displ)



# # print("move2_accelStopTime",move2_accelStopTime)
# # print("move2_steadyTime",move2_steadyTime)
# # print("move2_decelStartTime",move2_decelStartTime)
# # print("move2_startVelocity",move2_startVelocity)
# # print("move2_t1",move2_t1)
# # print("move2_t2",move2_t2)
# # print("move2_t3",move2_t3)
# # print("move2_T1",move2_T1)
# # print("move2_T2",move2_T2)
# # print("move2_F",move2_F)
# # print("move2_s_1e,",move2_s_1e)
# # print("move2_s_2e",move2_s_2e)        


# print(len(time))
# print(len(velocity))
# print(len(distance))
# print(len(displacement))

# fig, axs = plt.subplots(3)
# # fig.suptitle('Motion Profiles')

# axs[0].plot(time, velocity)
# axs[0].set_title("Velocity Profile")
# axs[0].set_xlabel("time (s)")
# axs[0].set_ylabel("velocity (mm/s)")
# axs[0].axvline(x=accelStopTime, color='r', linestyle='dashed')
# axs[0].axvline(x=decelStartTime, color='g', linestyle='dashed')
# axs[0].axvline(x=move2_accelStopTime+totalTime_move1, color='b', linestyle='dashed')
# axs[0].axvline(x=move2_decelStartTime+totalTime_move1, color='c', linestyle='dashed')

# axs[1].plot(time, distance)
# axs[1].set_title("Distance Profile")
# axs[1].set_xlabel("time (s)")
# axs[1].set_ylabel("Distance (mm)")
# axs[1].axvline(x=accelStopTime, color='r', linestyle='dashed')
# axs[1].axvline(x=decelStartTime, color='g', linestyle='dashed')
# axs[1].axvline(x=move2_accelStopTime+totalTime_move1, color='b', linestyle='dashed')
# axs[1].axvline(x=move2_decelStartTime+totalTime_move1, color='c', linestyle='dashed')

# axs[2].plot(time, displacement)
# axs[2].set_title("Displacement Profile")
# axs[2].set_xlabel("time (s)")
# axs[2].set_ylabel("Displacement (mm)")
# axs[2].axvline(x=accelStopTime, color='r', linestyle='dashed')
# axs[2].axvline(x=decelStartTime, color='g', linestyle='dashed')
# axs[2].axvline(x=move2_accelStopTime+totalTime_move1, color='b', linestyle='dashed')
# axs[2].axvline(x=move2_decelStartTime+totalTime_move1, color='c', linestyle='dashed')




# fig.tight_layout(rect=[0, 0.001, 1, 1.03])
# plt.show()