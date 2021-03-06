#include <pthread.h>
#include <unistd.h>
#include <stdio.h>

// Example 1 Stuff
void wait(void)
{
    sleep(2);

    printf("Done.\n");
}

// Example 2 Stuff
pthread_mutex_t lock;
int j;

void do_process()
{
    pthread_mutex_lock(&lock);
    int i = 0;

    j++;

    while(i < 5)
    {
        printf("%d", j);
        sleep(1);

        i++;
    }

    printf("...Done\n");

    pthread_mutex_unlock(&lock);
}

int main(void)
{
  // Example 1
  {
    pthread_t thread;
    int err;

    err = pthread_create(&thread, NULL, wait, NULL);

    if (err)
    {
        printf("An error occured: %d", err);
        return 1;
    }

    printf("Waiting for the thread to end...\n");

    pthread_join(thread, NULL);

    printf("Thread ended.\n");    
  }

  // Example 2
  {
    int err;
    pthread_t t1, t2;

    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("Mutex initialization failed.\n");
        return 1;
    }

    j = 0;

    pthread_create(&t1, NULL, do_process, NULL);
    pthread_create(&t2, NULL, do_process, NULL);

    pthread_join(t1, NULL);
    pthread_join(t2, NULL);
  }

    return 0;
}
