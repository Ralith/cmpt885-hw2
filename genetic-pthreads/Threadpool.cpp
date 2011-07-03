#include "Threadpool.h"

extern "C" void *_threadBody(void *self) {
  Threadpool *pool = static_cast<Threadpool*>(self);

  Task *t;
  while((t = pool->tasks.get())) {
    t->run();
    sem_post(&pool->sem);
  }

  return 0;
}

Threadpool::Threadpool(unsigned size_) :
  size(size_), threads(new pthread_t[size_]),
  tasks(size_), taskCount(0) {
  sem_init(&sem, 0, 0);
  for(unsigned i = 0; i < size; ++i) {
    pthread_create(&threads[i], NULL, _threadBody, this);
  }
}
Threadpool::~Threadpool() {
  for(unsigned i = 0; i < size; ++size) {
    tasks.put(0);
  }
  for(unsigned i = 0; i < size; ++size) {
    pthread_join(threads[i], 0);
  }
  delete[] threads;
  sem_destroy(&sem);
}

void Threadpool::addTask(Task *t) {
  tasks.put(t);
  ++taskCount;
}

void Threadpool::join() {
  for(unsigned i = 0; i < taskCount; ++i) sem_wait(&sem);
  taskCount = 0;
}
