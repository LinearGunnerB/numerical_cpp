#include <string>
#include <fstream>
#include <cmath>

void save(std::string fname, std::string line) {
    ofstream myfile (fname, ios::out);
    myfile << line;
    myfile.close();
}

void append(std::string fname, std::string line) {
    ofstream myfile (fname, ios::out | ios::app );
    myfile << line;
    if (line.compare("") == 0)
        myfile.close();
}

std::string math2D(std::string fname, double x) {
    // takes function name, and if it exists, evaluates the function at x
    std::string line = "";
    if (!fname.compare("sinh")) {
        line = std::to_string(x) + "    " + std::to_string( sinh(x) );
    }
    else if (!fname.compare("cosh")) {
        line = std::to_string(x) + "    " + std::to_string( cosh(x) );
    }
    // Add your functions here
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
    system(command.c_str());
}