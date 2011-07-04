#ifndef _THREADPOOL_H_
#define _THREADPOOL_H_

#include <pthread.h>
#include <semaphore.h>
#include "ringbuffer.h"

extern "C" void *_threadBody(void *self);

class Task {
  friend void *_threadBody(void*);
  
protected:
  bool done;
  
public:
  Task() : done(false) {};
  virtual ~Task() {};

  virtual void run(unsigned thread)=0;
};

class Threadpool {
  friend void *_threadBody(void*);

protected:
  unsigned size;
  pthread_t *threads;
  RingBuffer<Task *> tasks;
  sem_t sem;
  unsigned taskCount;

  struct threadarg {
    Threadpool *instance;
    unsigned id;
  };
  threadarg *args;

  int freeThread();

public:
  Threadpool(unsigned size);
  ~Threadpool();
  
  void addTask(Task *t);
  void join();
};

#endif
