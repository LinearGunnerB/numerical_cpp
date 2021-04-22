#include <iostream>
#include <cstdlib> // EXIT_SUCCESS
#include <cassert> // assert
#include <thread>
#include <chrono>
#include <vector>

using namespace std;

void execute(double (*myfunc) (double, double), double a, double b  )
{
  std::cout << myfunc(a, b) << "\n";
}

double add(double a, double b)
{
  return (a+b);
}

void foo(int a)
{
    std::cout << "foo: " << a << "\n";
}

void fuz(int a, int b)
{
    std::cout << "fuz: (" << a << ", " << b << ")\n";
}

class Bar
{
public:
    void foo(int a)
    {
        std::cout << a << '\n';
    }
};

class BarFoo
{
public:
    void operator()(int a)
    {
        std::cout << a << '\n';
    }
};


template <typename T>
void sum_vector(std::vector<T> *v, T *sum)
{
  for (typename std::vector<T>::iterator it = v->begin(); it != v->end(); it++)
  {
    cout << "it : " << *it << ", ";
    *sum += *it;
    cout << "\n";
  }
}

template <typename T>
void add_two_elements(T &a, T &b, T &c)
{
   c = a + b;
}

template <typename T=double>
void addRef(T& a, T& b, T& c)
{
  c = a + b;
}

template <typename T>
void add_with_threads(size_t numThreads, vector<T>& a, vector<T> &b, vector<T> &c)
{
  static const size_t MAX_THREADS_ALLOWED = 100;
  assert(numThreads <= MAX_THREADS_ALLOWED);

  std::thread myThread[numThreads];
  size_t tid = 0;
  size_t vid = 0;
  
  if ( a.size() > numThreads )
  {  
    //size_t inc = a.size() / numThreads;
    //size_t rem = vid % numThreads;
    for ( size_t vid = 0; vid < a.size(); vid++ )
    {
      for (size_t tid = 0; tid < numThreads; tid++ ) 
      {
        myThread[tid] = std::thread(add_two_elements<double>, std::ref(a[vid]), std::ref(b[vid]), std::ref(c[vid])); // function template needs explicit type
      }
      
      for (size_t tid = 0; tid < numThreads; tid++ ) 
      {
        myThread[tid].join();
      }
      
    }
    
  } // end of if
  else 
  {
    for (size_t tid = 0; tid < numThreads; tid++ ) 
    {
      if ( vid < a.size() ) {
        myThread[tid] = std::thread(add_two_elements<double>, std::ref(a[vid]), std::ref(b[vid]), std::ref(c[vid])); // function template needs explicit type
        vid += 1;
      }
    }
    for (size_t i=0; i < numThreads; i++)
    {
      myThread[i].join();
    }
  } // end of else
  
  //std::thread myThread(add_two_elements<double>, std::ref(a[i]), std::ref(b[i]), std::ref(c[i])); // WORKS for single add
  //myThread.join();  // WORKS for single add
}

template <typename T>
void add_chunk(size_t idxStart, size_t idxEnd, vector<T>& a, vector<T> &b, vector<T> &c)
{
  if ( idxEnd > a.size() ) { idxEnd = a.size() - 1; }
    for ( size_t i = idxStart;  i < idxEnd; i++ )
    {
      c[i] = a[i] + b[i];
    }
}

template <typename T>
void add_chunk_with_few_threads(size_t numThreads, size_t chunkSize, vector<T>& a, vector<T> &b, vector<T> &c)
{
  static const size_t MAX_THREADS_ALLOWED = 32;
  assert(numThreads <= MAX_THREADS_ALLOWED && a.size() == b.size() && b.size() == c.size() );
#if DEBUG
  cout << "Past Assert" << endl;
#endif

  std::thread myThread[numThreads];
  size_t tid = 0; // need tid at highest scope
  size_t vid = 0;
  size_t joinIdxEnd = 0;
  //size_t chunkIdx = chunkSize - 1; // 250  becoms 249 because of zero indexing
  
  //if ( chunkSize * numThreads <= a.size() )
  //{  
    size_t numLoops = a.size() / (chunkSize * numThreads) + 1; // + 1 for remainder if any
#if DEBUG
    cout << "numLoops : " << numLoops << endl;
#endif
    //size_t rem = vid % numThreads;
    //for ( size_t vid = 0; vid < a.size(); vid+= chunkSize * numThreads )
    for ( size_t n = 0; n < numLoops; n++)
    {
      for (tid = 0; tid < numThreads; ++tid ) 
      {
#if DEBUG
        cout << "numLoop [  " << n << " ] , thread [ " << tid << " ]" << endl;
#endif
        
        if ( vid >= a.size() ) // need this so when we have exactly enough threads to finish in numLoop 0, that we do not try to join extra non-existent threads
        {
          if ( tid > 0 ) { tid -=1; }   // when we use extra threads in 2nd numLoop we need to -1 from tid in join
          break;
        }
#if DEBUG
        cout << "vid : " << vid << endl;
        cout << "vid + chunkSize: " << vid + chunkSize << endl;
#endif
        myThread[tid] = std::thread(add_chunk<double>, vid, vid + chunkSize, std::ref(a), std::ref(b), std::ref(c)); // function template needs explicit type
        vid += chunkSize;
#if DEBUG
        cout << "numLoop [  " << n << " ] , thread [ " << tid << " ] : " << "bottom of loop" << endl;
#endif
      }
#if DEBUG
      cout << "+++++\n";
      cout << "tid : " << tid << endl;
      cout << "joinIdxEnd : " << joinIdxEnd << endl;
#endif
      if ( tid == numThreads )
      { joinIdxEnd = numThreads;  }
      else if ( tid == 0 )
      { joinIdxEnd = 1;  }
      else
      { joinIdxEnd = tid ; }
#if DEBUG  
      cout << "---------\n";
      cout << "tid : " << tid << endl;
      cout << "joinIdxEnd : " << joinIdxEnd << endl;
#endif
      for ( size_t i = 0; i < joinIdxEnd; i++ )
      {
#if DEBUG        
        cout << "start join() : " << i << " , tid : " << tid << " , joinIdxEnd : " << joinIdxEnd << endl;
#endif
        myThread[i].join();
#if DEBUG
        cout << "end join() : " << i << endl;
#endif
      }
      if ( vid >= a.size() ) // need this so when we have exactly enough threads to finish in numLoop 0, that we do not try to join extra non-existent threads
      {
        break;
      }
      //cout << "tid : " << tid << endl;       
    }
    cout << "Completed : add_chunk_with_few_threads\n";
  //} // end of if
  /*
  else 
  {
    size_t tid;
    for ( tid = 0; tid < numThreads; tid++ ) 
    {
        myThread[tid] = std::thread(add_chunk<double>, vid, vid + chunkIdx, std::ref(a), std::ref(b), std::ref(c) ); // function template needs explicit type
        vid += chunkSize;
        if ( vid >= a.size() )
        {
          break;
        }
    }
    for (size_t i=0; i < tid; i++)  // tid instead of numThreads in case we break before tid reaches numThreads
    {
      myThread[i].join();
    }
  } // end of else
  */
}

template <typename T>
void print_vector(vector<T>& vector)
{
  for (auto &vi: vector)
  {
    cout << vi << ", ";
  }
  cout  << "\n";
}

// overload + for vectors
template <typename T>
std::vector<T> operator+(std::vector<T>& a, std::vector<T> &b)
{
  assert(a.size() == b.size());
  std::vector<T> c(a.size() );
  for ( size_t i = 0; i < a.size(); i++)
  {
    c[i] = a[i] + b[i];
  }
  return c;
}


int main()
{
  {
  // Function pointer review
  execute(add, 5, 4);

  // Array of pointer review
  // std::string strArr[5] = {"Boys", "girls", "Cats", "Dogs", "12345"};
  std::string strArr[5];
  strArr[0] = "Line 1";
  strArr[1] = "Roses are Red";
  strArr[2] = "My name is Brian";
  std::cout << strArr[0] << "\n" << strArr[1] << "\n" << strArr[2] << "\n";

  char *Arr[5];
  char greeting[] = "hello";
  char name[]     = "Brian";
  Arr[0] = &greeting[0];
  Arr[1] = name;
  cout << Arr[0] << ", " << Arr[1] << "\n";


  /////////////////////////////////////////////////////////////
  //    Free Function Example
  /////////////////////////////////////////////////////////////
  // Create and execute the thread
  std::thread thread(foo, 10); // foo is the function to execute, 10 is the
                               // argument to pass to it
  // Keep going; the thread is executed separately
  // Wait for the thread to finish; we stay here until it is done
  thread.join();

  std::thread thread2(fuz, 10, 20);
  thread2.join();

  /////////////////////////////////////////////////////////////
  //    Member Function Example
  /////////////////////////////////////////////////////////////
  Bar bar;
  // Create and execute the thread
  std::thread thread3(&Bar::foo, &bar, 10); // Pass 10 to member function
  // The member function will be executed in a separate thread
  // Wait for the thread to finish, this is a blocking operation
  thread3.join();

  /////////////////////////////////////////////////////////////
  //    Functor Function Example
  /////////////////////////////////////////////////////////////
  BarFoo barfoo;  
  // Create and execute the thread
  //std::thread thread4(&BarFoo::operator(), &barfoo, 10); // Pass 10 to functor object
  std::thread thread4(barfoo, 10);
  // The functor object will be executed in a separate thread
  // Wait for the thread to finish, this is a blocking operation
  thread4.join();

  /////////////////////////////////////////////////////////////
  //    Lambda Expression Example
  /////////////////////////////////////////////////////////////
  auto lambda = [](int a) { std::cout << "lambda : " << a << "\n"; };
  std::thread thread_L(lambda, 4);
  thread_L.join();
  }

  ////////////////////////////////////////////////////////////
  //  Add vectors with multiple threads by Brian :) 
  ////////////////////////////////////////////////////////////
  {
  vector<double> x = {1 , 2, 3, 4, 5, 6, 7, 8, 9, 10};
  vector<double> y = {1 , 2, 3, 4, 5, 6, 7, 8, 9, 10};
  vector<double> z(10,0.0); // z(length, initial value of each element)

  print_vector(x);
  print_vector(y);
  print_vector(z);

  addRef(x[0], y[0], z[0]);
  print_vector(z);
  z[0] = 0;

  std::thread testThread(foo, x[0]);
  testThread.join();

  double sum = 0;
  std::thread myThread(sum_vector<double>, &x, &sum);
  myThread.join();
  cout << "sum ( x ) : " << sum << endl;
  
  add_with_threads(10, x, y, z);
  print_vector(z);

  x.push_back(11.0);
  x.push_back(12.0);
  y.push_back(11.0);
  y.push_back(12.0);
  z.push_back(0.0);
  z.push_back(0.0);
  add_with_threads(12, x, y, z);
  print_vector(z);
  }

  /////////////////////////////////////////////////////////////////
  // Large vector addition with multiple threads
  /////////////////////////////////////////////////////////////////
  size_t N = 100'000'000;
  size_t usPerSec = 1000000;
  vector<double> X( N, 3.14);
  vector<double> Y( N, 2.5);
  vector<double> Z( N, 0.0);
  //print_vector(X);
  //print_vector(Y);

  using namespace std::chrono;
  ///////////////////////////////////////////////////
  //   Single Thread Example
  ///////////////////////////////////////////////////
  {
    cout << "\n\n+++++++   Single Thread Vector Addition  +++++++++\n";
    cout << "Vector Size N : " << N << endl;
    auto start = high_resolution_clock::now();
    Z = X + Y;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "seconds : " <<  (double)duration.count() / (double)usPerSec << endl; 
    //print_vector(Z);
    cout << "Z[200] : " << Z[200] << endl;
    std::fill(Z.begin(), Z.end(), 0.0);
    cout << "Z[200] : " << Z[200] << endl;
  }

  ///////////////////////////////////////////////////
  //   Multiple Thread Example
  ///////////////////////////////////////////////////
  // Large vector addition with Chunks and fewer threads
  
  cout << "\n\n\n\n\n\n\n\n+++++++   add_chunk_with_few_threads result  +++++++\n";
  cout << "aka multi thread example :\n" ;
  {
    auto start = high_resolution_clock::now();
    size_t numThreads = 12;
    size_t myChunk = N / numThreads;
    
    cout << "N : " << N << endl;
    cout << "numThreads : " << numThreads << endl;
    cout << "myChunk : " << myChunk << " , (double)myChunk : " << (double)(N) / (double)numThreads << endl;

    add_chunk_with_few_threads(numThreads, myChunk, X, Y, Z);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "seconds : " << (double)duration.count() / (double)usPerSec << endl; 
    //print_vector(Z);
    //cout << "Z[200] : " << Z[200] << endl;
  }

  ///////////////////////////////////////////////////
  //   Multiple Thread Example
  ///////////////////////////////////////////////////
  //  Is slower than single thread because too many threads are created with small chunkSize
  /*
  std::fill(Z.begin(), Z.end(), 0.0);
  // Large vector addition with more threads
  cout << "\n\n\n\n\n\n\n\n add_with_threads result : \n";
  {
    auto start = high_resolution_clock::now();
    add_with_threads(8, X, Y, Z2);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "micro seconds :" << duration.count() << endl; 
    //print_vector(Z2);
  }
  */
  

  return (EXIT_SUCCESS);
}