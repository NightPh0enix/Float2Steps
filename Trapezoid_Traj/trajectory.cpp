
#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <string> 
#include "trajectory.h"


int decimal_digits(double sampling_interval) //number of digits in the decimal part of the sampling intervals i.e 0.1   0.01    0.001
{
    int rounding=0; 
    std::string interval_string = std::to_string(sampling_interval);
    for(int i=0;i<interval_string.length();i++)
    {
        rounding = rounding+1;

        if (interval_string[i]=='.')
        {
            rounding=0;
        }

        if (interval_string[i]=='1')
        {
            break;
        }
    }
    return rounding;
}

double convertTimeToSeconds(double time,Setup setup)
{
    time=time*setup.clocksToSecMultiplier;
    // time = roundf(time * pow(10,setup.rounding)) / pow(10,setup.rounding);
    return time;
}

double convertAccToSec(double accel,Setup setup)
{
    accel=accel*setup.clockRate*setup.clockRate;
    // accel = roundf(accel * pow(10,setup.rounding)) / pow(10,setup.rounding);
    return accel;
}

double convertVelToSec(double vel,Setup setup)
{
    vel=vel*setup.clockRate;
    //vel = roundf(vel * pow(10,setup.rounding)) / pow(10,setup.rounding);
    return vel;
}

std::vector<double> makeTimevector(std::vector<double> t_vector, double endtime, double signif_digits)
{
    // std::cout<<"endtime "<<endtime<<"\n";

    int loop_constant = ceil(endtime * (std::pow(10,decimal_digits(signif_digits))));


    for(int i=0;i<=loop_constant;i++)
    {
        float t = i/(std::pow(10,decimal_digits(signif_digits)));
        t_vector.push_back(t);
        //std::cout<<"t : "<<t <<"\n";

    }
    // std::cout<<"endtime"<<endtime;
    std::cout<<"last elem in time vec"<<t_vector.back();
    std::cout<<"len "<<t_vector.size();
    return t_vector;
}

std::vector<double> makeVelocityVector(std::vector<double> v_vector, std::vector<double> t_vector, MoveParameters m, Setup s)
{
    for (int t=0; t<=t_vector.size()-1;t++)
    {
        if (t_vector[t]>=0 and t_vector[t]<=m.t1)
        {
            double vel = m.startVelocity + s.acceleration * t_vector[t];
            v_vector.push_back(vel);
            std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }

        else if (t_vector[t]>=m.t1 and t_vector[t]<=m.t2)
        {
            double vel = m.F;
            v_vector.push_back(vel);
            std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }

        else if (t_vector[t]>=m.t2)
        {
            double vel = m.f_2e + s.deceleration * (t_vector[t]-m.t2);
            v_vector.push_back(vel);
            std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }
    }

    return v_vector;
}

std::vector<double> makeDistanceVector(std::vector<double> d_vector, std::vector<double> t_vector, MoveParameters m, Setup s)
{

    // std::cout<<"m.s_2e "<<m.s_2e<<"\n";
    // std::cout<<"m.s_1e "<<m.s_1e<<"\n";
    // std::cout<<"m.startDist "<<m.startDist<<"\n";
    // std::cout<<"m.startVelocity "<<m.startVelocity<<"\n";
    // std::cout<<"s.acceleration "<<s.acceleration<<"\n";
    std::cout<<"\n";
    for (int j=0; j<=t_vector.size()-1;j++)
    {
        // std::cout<<" j "<<j<<"\n";
        // std::cout<<" time "<<t_vector[j]<<"\n";
        //std::cout<<" C "<<t_vector[0]<<"\n";

        if (t_vector[j]>=0 and t_vector[j]<m.t1)
        {
            //std::cout<<"A";
            double dist = m.startDist + m.startVelocity * t_vector[j] + 0.5 * s.acceleration * t_vector[j] * t_vector[j];
            //std::cout << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
            d_vector.push_back(dist);
        }
        else if (t_vector[j]>=m.t1 and t_vector[j]<m.t2)
        {
            //std::cout<<"B";
            double dist = m.s_1e + m.F * (t_vector[j]-m.t1);
            //std::cout << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
 
            
            d_vector.push_back(dist);
        }
        else if (t_vector[j]>=m.t2)
        {
            //std::cout<<"C";
            double dist = m.s_2e + m.f_2e * (t_vector[j]-m.t2) + 0.5 * s.deceleration * (t_vector[j]-m.t2) * (t_vector[j]-m.t2);
            //std::cout << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
            d_vector.push_back(dist);

        }
    }
    return d_vector;
}


int main() 
{
    MoveParameters move1;
    MoveParameters move2;
    Setup setup;
    
    setup.acceleration = convertAccToSec(setup.acceleration,setup);
    setup.deceleration = setup.acceleration*(-1);

    move1.startDist = 0;
    move1.startVelocity = 0;

    move1.accelStopTime = 93750.000000;
    move1.steadyTime = 93750.000000;
    move1.decelStartTime = 187500.000000;
    move1.moveTime = 281250.000;

    move1.accelStopDistance= 2.500000;
    move1.steadyDistance = 5.000000;
    move1.decelStartDistance = 7.500000;

    move1.startVelocity = convertVelToSec(move1.startVelocity,setup);
    move1.accelStopTime = convertTimeToSeconds(move1.accelStopTime,setup);
    move1.steadyTime = convertTimeToSeconds(move1.steadyTime,setup);
    move1.decelStartTime = convertTimeToSeconds(move1.decelStartTime,setup);
    move1.moveTime = convertTimeToSeconds(move1.moveTime,setup);

    move1.t1 = move1.accelStopTime;
    move1.t2 = move1.decelStartTime;
    move1.t3 = move1.moveTime;
    move1.T1 = move1.accelStopTime;
    move1.T2 = move1.steadyTime;
    move1.T3 = move1.moveTime-move1.decelStartTime;

    move1.F = move1.f_1e = move1.f_2e = move1.startVelocity + setup.acceleration * move1.T1;

    move1.s_1e = move1.startDist + move1.startVelocity*move1.T1 + 0.5*setup.acceleration*move1.T1*move1.T1;
    move1.s_2e = move1.s_1e + move1.F*move1.T2;

    timevector = makeTimevector(timevector,move1.moveTime,setup.time_interval); // make time vector 
    velocityVector = makeVelocityVector(velocityVector,timevector,move1,setup);
    distanceVector = makeDistanceVector(distanceVector,timevector,move1,setup);

    // std::cout<<"last distamce "<<distanceVector.back();

    for (int c=0; c<=timevector.size()-1;c++)
    {
        std::cout<<"TIME : "<<timevector[c]<<"\n";
    }

    for (int c=0; c<=velocityVector.size()-1;c++)
    {
        std::cout<<"VEL : "<<velocityVector[c]<<"\n";
    }

    for (int c=0; c<=distanceVector.size()-1;c++)
    {
        std::cout<<"dis : "<<distanceVector[c]<<"\n";
    }
    //std::cout<<"time len : "<<timevector.size()<<"\n";

    // for (int c=0; c<=timevector.size()-1;c++)
    // {
    //     std::cout<<"dist : "<<distanceVector[c]<<"\n";
    // }
    // std::cout<<"move1.accelStopTime "<<move1.accelStopTime<<"\n";
    // std::cout<<"move1.steadyTime "<<move1.steadyTime<<"\n";
    // std::cout<<"move1.decelStartTime "<<move1.decelStartTime<<"\n";
    // std::cout<<"move1.moveTime "<<move1.moveTime<<"\n";
    // std::cout<<"move1.startVelocity :" <<move1.startVelocity<<"\n";
    // std::cout<<"setup.acceleration : " <<setup.acceleration<<"\n";
    // std::cout<<"move1.T1 : " <<move1.T1<<"\n";
    // std::cout<<"move1.T2 : " <<move1.T2<<"\n";
    // std::cout<<"move1.T3 : " <<move1.T3<<"\n";
    // std::cout<<"move1.t1 : " <<move1.t1<<"\n";
    //std::cout<<"move1.t2 : " <<move1.t2<<"\n";


// # ###################################   MOVE 2   ###################################

    move2.startDist = distanceVector.back() ;
    move2.startVelocity = velocityVector.back();

    move2.accelStopTime = 93750.000000;
    move2.steadyTime = 843750.000000;
    move2.decelStartTime = 937500.000000;
    move2.moveTime = 1031250.000;

    move2.accelStopDistance= 2.500000;
    move2.steadyDistance = 45.056252;
    move2.decelStartDistance = 47.556252;

    move2.startVelocity = convertVelToSec(move2.startVelocity,setup);
    move2.accelStopTime = convertTimeToSeconds(move2.accelStopTime,setup);
    move2.steadyTime = convertTimeToSeconds(move2.steadyTime,setup);
    move2.decelStartTime = convertTimeToSeconds(move2.decelStartTime,setup);
    move2.moveTime = convertTimeToSeconds(move2.moveTime,setup);

    move2.t1 = move2.accelStopTime;
    move2.t2 = move2.decelStartTime;
    move2.t3 = move2.moveTime;
    move2.T1 = move2.accelStopTime;
    move2.T2 = move2.steadyTime;
    move2.T3 = move2.moveTime-move2.decelStartTime;

    move2.F = move2.f_1e = move2.f_2e = move2.startVelocity + setup.acceleration * move2.T1;

    move2.s_1e = move2.startDist + move2.startVelocity*move2.T1 + 0.5*setup.acceleration*move2.T1*move2.T1;
    move2.s_2e = move2.s_1e + move2.F*move2.T2;



    return 0;



}

// ###################################    Parameters subject to change    ###################################

// ###################################    dist in mm     ################################### 

// ###################################    Times in clocks    ###################################  
 
// ###################################    Times in seconds    ###################################

// totalTime=round(totalTime*clocksToSecMultiplier, rounding)

// totalTime_move1=round(totalTime_move1*clocksToSecMultiplier, rounding)
// totalTime_move2=round(totalTime_move2*clocksToSecMultiplier, rounding)

// accelStopTime=round(accelStopTime*clocksToSecMultiplier, rounding)
// steadyTime=round(steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
// decelStartTime=round(decelStartTime*clocksToSecMultiplier, rounding)


// acceleration = round(acceleration*clockRate*clockRate,rounding)
// deceleration =acceleration*(-1)

// startVelocity=round(startVelocity*clockRate,rounding)

// t1=accelStopTime
// t2=decelStartTime
// t3=totalTime_move1

// T1=accelStopTime
// T2=steadyTime=decelStartTime-accelStopTime
// T3=totalTime_move1-decelStartTime


// ###################################    setting up time list    ################################### 
// loop_iterator=totalTime+time_interval
// # loop_iterator=totalTime_move1+time_interval
// time=[] 
// for i in np.arange(0,loop_iterator,time_interval):
//     # time.append(round(i,rounding)) # increments in time will be in 0.1 sec step i.e. 0.1 0.2 0.3
//     time.append(round(i,rounding))

// ###################################    velocity profile    ###################################

// f_2e=F=f_1e=startVelocity+acceleration*T1


// velocity=[]
// loopcounter=0
// for t in time:
//     loopcounter=loopcounter+1
//     if (t<=totalTime_move1):
//         # print("T value is ",t)
//         if (t>=0 and t<=t1):
//             # velocity.append("AA")
//             vel = startVelocity + acceleration * t
//             velocity.append(vel)
//             # print("A",t)
//         elif (t>=t1 and t<=t2):
//             # velocity.append("BB")
//             vel=F
//             velocity.append(vel)
//             # print("B",t)
//         elif (t>=t2 and t<=t3):
//             # velocity.append("CC")
//             vel = f_2e + deceleration * (t-t2)
//             velocity.append(vel)
//             # print("C",t)





// # ###################################    distance profile    ###################################

// s_1e = startDist + startVelocity*T1 + 0.5*acceleration*T1**2
// s_2e = s_1e + F*T2
// # print("t1 t2 t3",t1,t2,t3)


// distance=[]
// displacement=[]

// for t in time:
//     if (t<=totalTime_move1):
//         if (t>=0 and t<=t1):
//             dist = startDist+startVelocity*t + 0.5*acceleration*t**2
//             distance.append(dist)
//         elif (t>=t1 and t<=t2):
//             dist = s_1e+F*(t-t1)
//             distance.append(dist)
//         elif (t>=t2 and t<=t3):
//             dist = s_2e +f_2e*(t-t2) + 0.5*deceleration*(t-t2)**2
//             distance.append(dist)

// displacement=distance.copy()


// if displacement==distance:
//     print("NONONO")

// print("totalTime",totalTime)
// print("totalTime_move1",totalTime_move1)
// print("totalTime_move2",totalTime_move2)
// print("accelStopTime",accelStopTime)
// print("steadyTime",steadyTime)
// print("decelStartTime",decelStartTime)
// print("startVelocity",startVelocity)
// print("t1",t1)
// print("t2",t2)
// print("t3",t3)
// print("T1",T1)
// print("T2",T2)
// print("T3",T3)
// print("F, f_1e, f_2e", F)


// # ###################################   MOVE 2   ###################################


// move2_startVelocity=velocity[-1]
// move2_startDistance=distance[-1]





// move2_accelStopTime=79687.500000
// move2_steadyTime=1032304.687500
// move2_decelStartTime=1111992.250000
// move2_totaltime=1205742


// move2_accelStopTime=round(move2_accelStopTime*clocksToSecMultiplier, rounding)
// move2_steadyTime=round(move2_steadyTime*clocksToSecMultiplier, rounding) #for how long are we moving with const speed
// move2_decelStartTime=round(move2_decelStartTime*clocksToSecMultiplier, rounding)
// move2_totaltime=round(move2_totaltime*clocksToSecMultiplier, rounding)

// move2_t1=move2_accelStopTime
// move2_t2=move2_decelStartTime
// move2_t3=move2_totaltime

// move2_T1=move2_accelStopTime
// move2_T2=move2_steadyTime

// move2_loop_iterator=totalTime_move2+time_interval

// move2_F=move2_f_2e=move2_f_1e=move2_startVelocity+acceleration*move2_T1
// move2_s_1e = move2_startDistance + move2_startVelocity*move2_T1 + 0.5*acceleration*move2_T1**2
// move2_s_2e = move2_s_1e + move2_F*move2_T2

// move2_timearray=[] 

// for k in np.arange(0,move2_loop_iterator,time_interval):
//     move2_timearray.append(round(k,rounding))


// for q in move2_timearray:
    
//     if (q>0 and q<=move2_t1):
//         # velocity.append("A")
//         vel = startVelocity + acceleration * q
//         velocity.append(vel)

//     elif (q>=move2_t1 and q<=move2_t2):
//         # velocity.append("B")
//         vel=move2_F
//         velocity.append(vel)

//     elif (q>=move2_t2 and q<=move2_t3):
//         # velocity.append("C")
//         vel = f_2e + deceleration * (q-move2_t2)
//         velocity.append(vel)


// for p in move2_timearray:

//     if (p>0 and p<=move2_t1):
//         # distance.append("A")
//         dist = move2_startDistance+move2_startVelocity*p + 0.5*acceleration*p**2
//         distance.append(dist)

//     elif (p>move2_t1 and p<=move2_t2):
//         # distance.append("B")
//         dist = move2_s_1e + move2_F * (p-move2_t1)
//         distance.append(dist)

//     elif (p>move2_t2 and p<=move2_t3):
//         # distance.append("C")
//         dist = move2_s_2e + move2_f_2e*(p-move2_t2) + 0.5*deceleration*(p-move2_t2)**2
//         distance.append(dist)


// for z in move2_timearray:

//     if (z>0 and z<=move2_t1):
//         # displacement.append("A")
//         displ = move2_startDistance - move2_startVelocity*z - 0.5*acceleration*z**2
//         displacement.append(displ)
//         last_dist=displacement[-1]

//     elif (z>move2_t1 and z<=move2_t2):
//         # displacement.append("B")
//         displ = last_dist - move2_F * (z-move2_t1)
//         displacement.append(displ)
//         last_dist_steady=displacement[-1]

//     elif (z>move2_t2 and z<=move2_t3):
//         # displacement.append("C")
//         displ = last_dist_steady - move2_f_2e*(z-move2_t2) - 0.5*deceleration*(z-move2_t2)**2
//         displacement.append(displ)