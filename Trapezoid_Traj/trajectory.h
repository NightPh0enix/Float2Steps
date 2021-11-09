//Constant params





std::vector<double> timevector;
std::vector<double> velocityVector;
std::vector<double> distanceVector;

class MoveParameters
{
    
    public:
        double startVelocity;
        double startDist;
        double accelStopDistance;
        double steadyDistance;
        double decelStartDistance;
        double accelStopTime;
        double steadyTime;
        double decelStartTime;
        double moveTime;

        double F,f_1e,f_2e;

        double s_1e;
        double s_2e;



        // double totalTime = totalTime_move1 + totalTime_move2; this will be move1.movetime +  move2.movetime

        double t1 = accelStopTime;
        double t2 = decelStartTime;
        double t3 = moveTime;

        double T1 = accelStopTime;
        double T2 = steadyTime=decelStartTime-accelStopTime;
        double T3 = moveTime-decelStartTime;



};

class Setup 
{
    public:
        friend class MoveParameters;
        //int rounding;
        double startTime=0.0;
        double clockRate=937500.0;
        double clocksToSecMultiplier=1/clockRate;
        double acceleration = 5.6889e-10; //mm/cks^2
        double deceleration;
        double time_interval=0.01;
};

