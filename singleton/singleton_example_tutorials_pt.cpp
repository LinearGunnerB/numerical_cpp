// [1] https://www.tutorialspoint.com/Explain-Cplusplus-Singleton-design-pattern

#include <iostream>

using namespace std;

// Can Static function can only access static variables?
static const double PI_76th = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
const double nonstatic_PI_76th = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
static double CalcPi(size_t NumTaylor)
{
  double arctan1 = 0.0;
  for (size_t k = 0; k < NumTaylor; k++)
  {
    if ( k % 2 )
      arctan1 += (-1.0) / ( 2.0*k + 1.0 );
    else
      arctan1 += (1.0) / ( 2.0*k + 1.0 );
  }

  arctan1 *= 4;
  cout << "arctan1 is : " << arctan1 << endl;
  //cout << "PI_76th is : " << PI_76th << endl;
  cout << "PI_76th is : " << nonstatic_PI_76th << endl;   // For free static functions they seem to access nonstatic PI variable, 
                                                          // but probably not true for static member functions of classes
  cout << "residual ( PI - arctan1 ) : " << (PI_76th - arctan1) << endl;
  return arctan1;
}


class A{
	int a;//non static
public:
  static int b;
	//Static function
	static int GetValue(){
		// a = 10;  // compile error
    b = 10;
		//return a;		// compile error
    return b;
	}
};
//initialize static variable
int A::b = 50;


class Singleton {
   static Singleton *instance;
   int data;
 
   // Private constructor so that no objects can be created.
   Singleton() {
      data = 0;
   }

   public:
   static Singleton *getInstance() {
      if (!instance)
      instance = new Singleton;
      return instance;
   }

   int getData() {
      return this -> data;
   }

   void setData(int data) {
      this -> data = data;
   }
};

//Initialize pointer to zero so that it can be initialized in first call to getInstance
Singleton *Singleton::instance = 0;


int main(int argc, char *argv[]){
   Singleton *s = s->getInstance();
   cout << s->getData() << endl;
   s->setData(100);
   cout << s->getData() << endl;

  std::cout.precision(76);
  if (argv[1] > 0)
    CalcPi(atoi(argv[1]));
  else
    CalcPi(10);

  cout << A::GetValue();

   return 0;
}
