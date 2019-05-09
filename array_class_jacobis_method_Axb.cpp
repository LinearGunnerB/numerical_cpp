#include<iostream>
#include<cmath>
#include<cstdlib> // atoi() EXIT_SUCCESS, EXIT_FAILURE, and exit()
#include<cstddef> // ptrdiff_t
#include<cstring> // strcpy
//#include<array>   // C++11
#include<vector> // C++03
#include<algorithm> // std::sort
#include <string>
#include <sstream>
using namespace std;

// Array class
class Array {
public:
    int rows, cols;
    double **array;
    Array(int m, int n) : rows(m) , cols(n)  // Constructor
    { 
        array = new double *[m];
        for (int i = 0; i < m; i++)
        {
            array[i] = new double[n];  //  array[count] = new int[cols];
        }
    }
    ~Array()  // Destructor
    {
        // Deallocate
        for (int i = 0; i < rows; i++ )
            delete[] array[i];
        delete[] array; //this needs to be done last
    }
    Array(const Array& m)  // Copy Constructor
    {
        rows = m.rows;
        cols = m.cols;
        array = new double*[rows];
    
        for (int i=0; i<rows; i++)
        {
            array[i] = new double[cols];
        }
    
        for (int i=0; i<rows; i++)
        {
            for (int j=0; j<cols; j++)
            {
                array[i][j] = m.array[i][j];
            }
        }
    }

    void resize(int m, int n)
    {
        Array temp(m,n);
        for (int i = 0; i < rows; i++)
        {
            for (int j=0; j < cols ; j++)
            {
                temp.array[j][i] = array[i][j];
            }
        }
        
        // print temp.array
        //cout << "temp = ";
        //temp.print_array();


        // deallocate array
        for (int i = 0; i < rows; i++ )
            delete[] array[i];
        delete[] array; //this needs to be done last

        //reallocate array, with new sizes
        array = new double *[m];
        for (int i = 0; i < m; i++)
        {
            array[i] = new double[n];  //  array[count] = new int[cols];
        }

        // reassign rows and cols
        rows = m;
        cols = n;

        // reassign values of temp to array
        for (int i = 0; i < rows; i++)
        {
            for (int j=0; j < cols ; j++)
            {
                array[i][j] = temp.array[i][j];
            }
        }
        //cout << "array = ";
        //print_array();

    } // end of resize


    double** operator + (Array& B) 
    {
        if ( rows !=  B.get_rows() || cols != B.get_cols() ) {
            // exit();
        }
        else {
            for (int i = 0; i < rows; i++) {
                for (int j=0; j < cols; j++) {
                    array[i][j] += B.array[i][j];
                    //cout << "array" << array[i][j] << endl;
                }
            }
        }
        return(this -> array);
    }

    //overload operator= (copy assignment operator) - does not create new matrix
    /*  
        http://www.cplusplus.com/forum/beginner/200934/   
        http://www.cplusplus.com/doc/tutorial/classes2/#copy_constructor
    */
    Array operator=(const Array& m)
    {
        if ( rows !=  m.get_rows() || cols != m.get_cols() ) {
            cout << "\nError: rhs rows & cols different than lhs array\n";
        }
        else {
            for (int i = 0; i < rows; i++) {
                for (int j=0; j < cols; j++) {
                    this->array[i][j] = m.array[i][j];
                }
            }
        } // end of else
        return(*this);
    }

    double** operator = (std::initializer_list<double> s ) 
    {
        double temp[s.size()];
        std::copy(s.begin(), s.end(), temp);

        //for (int i=0; i< s.size(); i++) { cout << temp[i] << " ";}

        int k = 0;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                array[i][j] = temp[k];
                //cout << "array ij = " << array[i][j] << endl;
                k += 1;
            }
        }

        return(this ->array);
    }
    

    void fill(double a, double b) 
    { 
        double inc = (b-a)/double(rows*cols);
        for (int i=0; i < rows; i++)  {
            for (int j=0; j< cols; j++) {
                array[i][j] = a + inc * double(i*cols + j);
            }

        }
    }

    double max()
    {
        //cout << "\nin max!" << endl;
        //cout << "array[0][0] = " << array[0][0] << endl;
        double temp = array[0][0];
        //cout << rows << " , " << cols << endl;

        for (int i=0; i < (rows); i++)
        {
            for (int j=0; j < (cols); j++)
            {
                //cout << "temp = " << temp << " , " << "array[i][j] = " << array[i][j] << endl;
                if ( (temp < array[i][j]) )
                {
                    //cout << "Duck Duck Goose:" << temp << endl;
                    temp = array[i][j];
                }
                else
                {
                    //temp = array[i][j];
                    //cout << "Red Riding Hood: " << temp << endl;
                    continue;
                }
            }
        }
        return(temp);
    }

    int get_rows() const { return this->rows; }
    int get_cols() const { return this->cols; }

    double element(int row, int col)
    {
        //return *(pt_a + row*cols + col);
        return (array[row][col]);
    }

    void print_array() 
    {
        cout << "\n";
        for (int i = 0; i < rows; i++ ) {
            for (int j = 0; j < cols - 1; j++) {
                cout << element(i, j) << " ";
            }
    
            cout << element(i, cols-1) << "\n";

        }
    }
};


Array Jacobi(Array& A, Array& x0, Array& b, int N=10, double tol=1e-4)
{
    // Jacobi's Method
    cout << "\nJacobi's Method:" << endl;
    
    // Check if x0, b, and A are correct shape
    if ( (x0.get_rows() == 1) && (x0.get_cols() > 1)  ) 
    {
        // transpose
        //cout << "in transpose: " << endl;
        //x0.print_array();
        x0.resize( x0.get_cols() , x0.get_rows() );
        //x0.print_array();
        cout << endl;
    } 

    Array x(x0.get_rows(), x0.get_cols() );
    x.fill(0,0);
    //cout << "x rows =  " << x.get_rows() << ", x cols = "<< x.get_cols() << ", x = ";

    int colOne = x0.get_cols() - 1;
    double sum = 0;
    double xmax = 0, x0max = 0, error = 0;

    int k = 0;
    while ( k < N ) 
    {
        for ( int i = 0; i < A.get_rows(); i++ )
        {
            sum = 0;
            for (int j = 0; j < A.get_cols(); j++ )
            {
                if ( i != j )
                    sum -= A.array[i][j] * (x0.array[j][colOne]); // minus is important
                else
                {
                    //continue;
                }
            }
            x.array[i][colOne] = (sum + (b.array[i][colOne])) / (A.array[i][i]);
        }

        
        cout << "x(" << k << ") = ";
        x.print_array();

        /* //Debug
        xmax  = x.max();
        cout << "\nx.max  = " << xmax << endl;

        x0max = x0.max();
        cout << "x0.max = " << x0max << endl;
        */

        xmax  = x.max();
        x0max = x0.max();
        error = xmax - x0max;
        if ( abs(error) < tol )
        {
            cout << "\nSuccess! Reached tol in "<< k <<" iterations, run complete! " << endl;
            break;
        }
        else
        {         
            cout << "\nerror = " << abs( error ) << ", tol = " << tol << endl;
            x0 = x;
        }
        ++k;
    }

    if ( k >= N )
    {
        cout << "\nReached max iterations: " << k << ", error is: " << abs(error) << endl;
    }
    return(x);
}


Array GuassSidel(Array& A, Array& x0, Array& b, int N=10, double tol=1e-4)
{
    // Guass-Sidel Method
    // Finish Later
    cout << "\nGuass-Sidel Method:" << endl;
    
    // Check if x0, b, and A are correct shape, if not then transpose 
    if ( (x0.get_rows() == 1) && (x0.get_cols() > 1)  ) 
    {
        // transpose
        x0.resize( x0.get_cols() , x0.get_rows() );
        cout << endl;
    } 

    Array x(x0.get_rows(), x0.get_cols() );
    x.fill(0,0);

    int colOne = x0.get_cols() - 1;
    double sum = 0;
    double xmax = 0, x0max = 0, error = 0;

    int k = 0;
    while ( k < N ) 
    {
        for ( int i = 0; i < A.get_rows(); i++ )
        {
            sum = 0;
            for (int j = 0; j < A.get_cols(); j++ )
            {
                if ( i != j )
                    sum -= A.array[i][j] * (x0.array[j][colOne]); // minus is important
                else
                {
                    //continue;
                }
            }
            x.array[i][colOne] = (sum + (b.array[i][colOne])) / (A.array[i][i]);
        }

        cout << "x(" << k << ") = ";
        x.print_array();

        xmax  = x.max();
        x0max = x0.max();
        error = xmax - x0max;
        if ( abs(error) < tol )
        {
            cout << "\nSuccess! Reached tol in "<< k <<" iterations, run complete! " << endl;
            break;
        }
        else
        {         
            cout << "\nerror = " << abs( error ) << ", tol = " << tol << endl;
            x0 = x;
        }
        ++k;
    }

    if ( k >= N )
    {
        cout << "\nReached max iterations: " << k << ", error is: " << abs(error) << endl;
    }
    return(x);
}




int main(int argc, char * argv[])
{
    // Reference: from https://stackoverflow.com/questions/2076624/c-matrix-class
    cout << "Starting program!" << endl;
    Array A(4,4); 
    cout << A.get_rows() << endl;
    cout << A.get_cols() << endl;
    
    A = { 10, -1, 2, 0, -1, 11, -1, 3, 2, -1, 10, -1, 0, 3, -1, 8 };
    cout << "\nA = ";
    A.print_array();

    Array x0(1,4);  
    x0.fill(0,0);
    cout << "\nx0 = ";
    x0.print_array();

    Array b(4,1);  
    b = { 6, 25, -11, 15};
    cout << "b =  ";
    b.print_array();

    Array x(4,1);

    // Calculate x with Jacobi's Method
    if (argc ==  3)
    {
        x = Jacobi(A,x0,b,atoi(argv[1]),atof(argv[2]));
    }
    else 
    {
        x = Jacobi(A,x0,b);
    }

    x.print_array();


    return(EXIT_SUCCESS);
}