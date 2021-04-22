// https://www.tutorialspoint.com/cplusplus/cpp_multithreading.htm
// compile with > g++ multithreadingFile.cpp -lprthread
#include <iostream>
#include <cstdlib> // intptr_t
#include <cstdint>
#include <string>
#include <pthread.h>
#include<windows.h>  // Sleep(milliseconds)  Win32 API
#include <unistd.h>  // sleep(), usleep()  windows and linux maybe?


using namespace std;
#define NUM_THREADS 5
const size_t us = 1;
const size_t ms = us * 1000;
const size_t  s = ms * 1000;

//////////////////////////////////////////////////////////////////
//   Creating Threads that print "Hello World"
//////////////////////////////////////////////////////////////////
void *PrintHello(void *threadid) {
   long tid;
   tid = (intptr_t)threadid;
   //tid = *( (long*)(&threadid) );
   cout << "Hello World! Thread ID, " << tid << endl;
   pthread_exit(NULL);
   return (nullptr);  // removes warning 
}

//////////////////////////////////////////////////////////////////
//   Passing Arguments to Threads
//////////////////////////////////////////////////////////////////

struct thread_data {
   int thread_id;
   //char *message;
   string message;
};


void *PrintData(void *threadarg) 
{
  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg;

  cout << "Thread ID : " << my_data->thread_id ;
  cout << " Message : " << my_data->message << endl;

  //usleep(100000);

  pthread_exit(NULL);
  return (nullptr);
}

void *wait(void *t) {
   int i;
   long tid;

   //tid = (long)t;  // error
   tid = (intptr_t)t;

   sleep(3);
   cout << "Sleeping in thread " << endl;
   cout << "Thread with id : " << tid << "  ...exiting " << endl;
   pthread_exit(NULL);
   return (nullptr);
}

int main () {
  /////////////////////////////////////////////////////////////////
  // Thread Example Hello World
  /////////////////////////////////////////////////////////////////
  {
   pthread_t threads[NUM_THREADS];
   int rc;
   int i;
   for( i = 0; i < NUM_THREADS; i++ ) {
      cout << "main() : creating thread, " << i << endl;
      //rc = pthread_create(&threads[i], NULL, PrintHello, (void*)i);
      rc = pthread_create(&threads[i], NULL, PrintHello, &i);
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
   }
   //pthread_exit(NULL);  // had to comment for next block of code to execute
   
  }
  //////////////////////////////////////////////////////////////////
  //   Passing Arguments to Threads
  //////////////////////////////////////////////////////////////////
  {
    //Sleep(500);
    usleep(ms);
    cout << "\n\n";

    pthread_t threads[NUM_THREADS];
    struct thread_data td[NUM_THREADS];
    int rc;
    int i;

    for( i = 0; i < NUM_THREADS; i++ ) 
    {
      cout <<"main() : creating thread, " << i << endl;
      td[i].thread_id = i;
      td[i].message = "This is message";
      rc = pthread_create(&threads[i], NULL, PrintData, (void *)&td[i]);
      
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
    }
    usleep(ms);
    //pthread_exit(NULL);
  }

  //////////////////////////////////////////////////////////////////
  //   Joining and Detaching Threads
  //////////////////////////////////////////////////////////////////
  /*    pthread_join (threadid, status)
        pthread_detach (threadid) 
        blocks the calling thread until the specified 'threadid' thread terminates
        one of its attributes defines whether it is joinable or detached
  */
  {
    usleep(ms);
    cout << "\n\n";

    int rc;
    int i;
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    void *status;

    // Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for( i = 0; i < NUM_THREADS; i++ ) {
      cout << "main() : creating thread, " << i << endl;
      //rc = pthread_create(&threads[i], &attr, wait, (void *)i );
      rc = pthread_create(&threads[i], &attr, wait, &i ); // removes warning on i
      if (rc) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
    }

    // free attribute and wait for the other threads
    pthread_attr_destroy(&attr);
    for( i = 0; i < NUM_THREADS; i++ ) {
      rc = pthread_join(threads[i], &status);
      if (rc) {
         cout << "Error:unable to join," << rc << endl;
         exit(-1);
      }
      cout << "Main: completed thread id :" << i ;
      cout << "  exiting with status :" << status << endl;
    }

    cout << "Main: program exiting." << endl;
    pthread_exit(NULL);
  }

  return(0);
}