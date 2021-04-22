#include<iostream>
#include<cstdlib> // EXIT_SUCCESS
//#include<iomanip> // for setprecision() and fixed
//#include<cstddef> // for operator sizeof() and the data type size_t
//#include<climits> // for INT_MAX, SHRT_MIN, etc.
//#include<ctime> // clock_t begin = clock() , CLOCKS_PER_SEC
#include<string>
#include <fstream>
//#define _USE_MATH_DEFINES // M_PI, M_E, and other math constants
//#include "math.h" // Also needed for the math constants
#include<cmath> // trig , hyperbolic functions
using namespace std;

const double c = 2.99792458e8; //m/s
const double mach1 = 343.0; //m/s
const double pi = 3.14159265359;

void msg(int n){
    switch(n) {
    case 0:
        cout << "C++ Program Start\n";
        break;
    case 1:
        cout << "Exiting program\n";
        break;
    default:
        cout << "Default msg... Choose msg 0 or 1 instead\n";
        break;
    }
}

void save(std::string fname, std::string line) {
    // std:cout << str1 << std::endl;  //test string input
    ofstream myfile (fname, ios::out);
    myfile << line;
    myfile.close();
}

void append(std::string fname, std::string line) {
    // std:cout << str1 << std::endl;  //test string input
    ofstream myfile (fname, ios::out | ios::app );
    myfile << line;
    if (line.compare("") == 0)
        myfile.close();
}

std::string math2D(std::string fname, double x) {
    // takes function name, and if it exists, evaluates the function at x
    std::string line = "";
    if (!fname.compare("sinh")) {
        //double y = sinh(x);
        line = std::to_string(x) + "    " + std::to_string( sinh(x) );
    }
    else if (!fname.compare("cosh")) {
        //double y = sinh(x);
        line = std::to_string(x) + "    " + std::to_string( cosh(x) );
    }
    else {
        line = "Unavailable function, choose again.\n";
    }
    return(line);
}

void func_a_to_b(std::string file, std::string func, double a, double b, int n) {
    //saves data for function from a to b with n points , to file fname
    double x;
    double inc = (b-a)/static_cast<double>(n);
    save(file,"#x \t\t y\n");  // put # first so python recognizes as comment
    for (int i=0; i < n; i++) {
        x = a + inc*i;
        append(file,math2D(func, x)+"\n");
    }
    append(file,""); // to close file
}

void run_python(string file="", string args="") {
    string command = "python " + file + " " + args;
    //cout << command;
    //cout << command.c_str();
    //cout << command.compare(command.c_str());
    //system(command.c_str());
    system(command.c_str());
}

int main()
{
    msg(0);

    save("data.txt","3.1415926    0.21415\n");
    append("data.txt","3.1415926    0.21415\n");
    append("data.txt",math2D("sinh",0.5));

    func_a_to_b("sinh.txt","sinh",-2*pi,2*pi,200);
    func_a_to_b("cosh.txt","cosh",-2*pi,2*pi,200);

    string args1 = "-f cosh.txt -lw 1 -xl x -yl y -sh True -t coshx -l cosh(x) -s cosh.png";
    string args2 = "-f sinh.txt -c r -lw 1 -xl x -yl y -t \"sinh(x)\" -l sinh(x) -sh True -s sinh.png";
    string args3 = "-f sinh.txt cosh.txt -s myplot.png -c b r -lw 1 1 -xl x -yl y -t \"sinh & cosh\" -l sinh(x) cosh(x) -dpi=1600 -sh True -sf False";

    run_python("plot2D.py",args1);
    run_python("plot2D.py",args2);
    run_python("plot2D.py",args3);

    msg(1);
    return(EXIT_SUCCESS);
}