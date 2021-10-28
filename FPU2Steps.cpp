
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream> 

uint32_t float2steps(float floatDistance, int StepsPerMm)
{
  float steps = floatDistance*StepsPerMm;
  
  std::cout <<int(steps)<< '\n';
  return steps;
}

void writeResult(std::vector<std::string> &result)  //for pass by reference look at the commented code below
{
  int EOLcount=0;
  std::ofstream merifile;
  merifile.open ("F2Step_result.txt"); // open file to write data
  for (size_t i = 0; i < result.size()-4; ++i) // -4 since the last 4 entries repeat 
  {
    EOLcount++;                                 // count which puts data in X Y Z E places 
    float my_res=std::stof(result[i]);               // string to float
    int steps2write = float2steps(my_res,80);  // fpu to mm
    if (EOLcount==1)
    {
      merifile <<"X"<<steps2write<<" ";
    }
    if (EOLcount==2)
    {
      merifile <<"Y"<<steps2write<<" ";
    }
    if (EOLcount==3)
    {
      merifile <<"Z"<<steps2write<<" ";
    }
    if (EOLcount==4)
    {
      merifile <<"E"<<steps2write<<'\n';
      EOLcount=0; 
    }  
  }
  merifile.close();
}

void readData()
{
  std::vector<std::string> result;
  std::ifstream myfile;
  std::string myline;
  myfile.open("desired.txt");
  if ( myfile.is_open() ) 
  {
    while (myfile) 
    { 
      std::getline (myfile, myline); // get a single row from the data
      std::stringstream ss(myline); // create a stream from the "extracted" line 
      while(ss.good())  // while(ss) ie. if the stream is true
      {
        std::string substr;
        std::getline(ss, substr, ','); // get the elements in the row seperated by ','
        result.push_back( substr ); //pushing each element in the vector 
      }
    }
  }
  else 
  {
    std::cout << "Couldn't open file\n";
  }
  
  writeResult(result); //passing the vector to be written in the desired format
}

int main ()
{ 
  readData();
  return 0;    
}


//#######################################################################################
//PASSING BY REFERENCE

// void function2(std::vector<std::string> &number)
// {
//   for (size_t i = 0; i < number.size(); ++i)
//   {
//     std::cout<<number[i]<<std::endl;
//   } 
  
// }
// void function1()
// {
//     std::vector<std::string> colour {"Blue", "Red", "Orange"};
//     function2(colour);
// }

// int main ()
// { 

//   function1();
//   return 0;    
// }
//#######################################################################################