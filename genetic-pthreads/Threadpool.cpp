#include "Threadpool.h"

extern "C" void *_threadBody(void *self) {
  Threadpool::threadarg *a = static_cast<Threadpool::threadarg*>(self);

  Task *t;
  while((t = a->instance->tasks.get())) {
    t->run(a->id);
    sem_post(&a->instance->sem);
  }

  return 0;
}

Threadpool::Threadpool(unsigned size_) :
  size(size_), threads(new pthread_t[size_]),
  tasks(size_), taskCount(0), args(new threadarg[size_]) {
  sem_init(&sem, 0, 0);
  for(unsigned i = 0; i < size; ++i) {
    args[i].instance = this;
    args[i].id = i;
    pthread_create(&threads[i], NULL, _threadBody, &args[i]);
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
  delete[] args;
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
