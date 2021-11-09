
#include <iostream>
#include <vector>
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
    return time;
}

double convertAccToSec(double accel,Setup setup)
{
    accel=accel*setup.clockRate*setup.clockRate;
    return accel;
}

double convertVelToSec(double vel,Setup setup)
{
    vel=vel*setup.clockRate;
    return vel;
}

std::vector<double> makeTimevector(std::vector<double> t_vector, double endtime, double signif_digits)
{
    int loop_constant = ceil(endtime * (std::pow(10,decimal_digits(signif_digits))));

    for(int i=0;i<=loop_constant;i++)
    {
        float t = i/(std::pow(10,decimal_digits(signif_digits)));
        t_vector.push_back(t);
    }
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
            //std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }

        else if (t_vector[t]>=m.t1 and t_vector[t]<=m.t2)
        {
            double vel = m.F;
            v_vector.push_back(vel);
            //std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }

        else if (t_vector[t]>=m.t2)
        {
            double vel = m.f_2e + s.deceleration * (t_vector[t]-m.t2);
            v_vector.push_back(vel);
            //std::cout << "t " << t << " time "<< t_vector[t] << " vel " << vel << "\n";
        }
    }

    return v_vector;
}

std::vector<double> makeDistanceVector(std::vector<double> d_vector, std::vector<double> t_vector, MoveParameters m, Setup s)
{
    for (int j=0; j<=t_vector.size()-1;j++)
    {
        if (t_vector[j]>=0 and t_vector[j]<m.t1)
        {

            double dist = m.startDist + m.startVelocity * t_vector[j] + 0.5 * s.acceleration * t_vector[j] * t_vector[j];
            //std::cout<< "A " << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
            d_vector.push_back(dist);
        }
        else if (t_vector[j]>=m.t1 and t_vector[j]<m.t2)
        {
            double dist = m.s_1e + m.F * (t_vector[j]-m.t1);
            //std::cout<< "B " << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
            d_vector.push_back(dist);
        }
        else if (t_vector[j]>=m.t2)
        {
            double dist = m.s_2e + m.f_2e * (t_vector[j]-m.t2) + 0.5 * s.deceleration * (t_vector[j]-m.t2) * (t_vector[j]-m.t2);
            //std::cout<< "C " << "j " << j << " time "<< t_vector[j] << " dist " << dist << "\n";
            d_vector.push_back(dist);
        }
    }
    return d_vector;
}

void printDebug(MoveParameters m, Setup s)
{

    std::cout<<"time len : "<<timevector.size()<<"\n";
    std::cout<<"accelStopTime : "<<m.accelStopTime<<"\n";
    std::cout<<"steadyTime : "<<m.steadyTime<<"\n";
    std::cout<<"decelStartTime : "<<m.decelStartTime<<"\n";
    std::cout<<"moveTime : "<<m.moveTime<<"\n";
    std::cout<<"startVelocity : "<<m.startVelocity<<"\n";
    std::cout<<"startDist : "<<m.startDist<<"\n";
    std::cout<<"setup.acceleration : "<<s.acceleration<<"\n";
    std::cout<<"T1 : " <<m.T1<<"\n";
    std::cout<<"T2 : " <<m.T2<<"\n";
    std::cout<<"T3 : " <<m.T3<<"\n";
    std::cout<<"t1 : " <<m.t1<<"\n";
    std::cout<<"t2 : " <<m.t2<<"\n";
    std::cout<<"F = f_1e = f_2e " <<m.F <<"\n";
    std::cout<<"s_1e : "<<m.s_1e<<"\n";
    std::cout<<"s_2e : "<<m.s_2e<<"\n";
   

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
        std::cout<<"DIS : "<<distanceVector[c]<<"\n";
    }

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
    move1.steadyTime = 140625.000000; 
    move1.decelStartTime = 234375.000000; 
    move1.moveTime = 234375.0;

    // move1.accelStopDistance= 2.500000;
    // move1.steadyDistance = 5.000000;
    // move1.decelStartDistance = 10;

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

    printDebug(move1,setup);

// # ###################################   MOVE 2   ###################################

    move2.startDist = distanceVector.back() ;
    move2.startVelocity = velocityVector.back();
    
    move2.accelStopTime = 0.000000;
    move2.steadyTime = 234375.000000;
    move2.decelStartTime = 234375.000000;
    move2.moveTime = 328125;

    // move2.accelStopDistance= 0.000000;
    // move2.steadyDistance = 12.500000;
    // move2.decelStartDistance = 12.500000;

    //move2.startVelocity = convertVelToSec(move2.startVelocity,setup);
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

    timevector.clear(); // IT IS IMPORTANT TO CLEAR THE TIME VECTOR BEFORE PROCEEDING WITH THE NEXT MOVE!!!!!!!

    timevector = makeTimevector(timevector,move2.moveTime,setup.time_interval); // make time vector for the new move so we will have to clear the previous time array
    velocityVector = makeVelocityVector(velocityVector,timevector,move2,setup);
    distanceVector = makeDistanceVector(distanceVector,timevector,move2,setup);

    printDebug(move1,setup);

    return 0;
}

