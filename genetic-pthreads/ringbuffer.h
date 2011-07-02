#ifndef _RINGBUFFER_H_
#define _RINGBUFFER_H_

#include <pthread.h>
#include <vector>
#include <string>

template<class T>
class RingBuffer {
private:
  T *buffer;
  unsigned tail, head, count, _size;
  pthread_cond_t cond;
  pthread_mutex_t mutex;

public:
  RingBuffer(unsigned size);
  ~RingBuffer();

  T get();
  void put(T value);
};

template<class T>
RingBuffer<T>::RingBuffer(unsigned size) :
  buffer(new T[size]),
  tail(0),
  head(0),
  count(0),
  _size(size) {
  pthread_mutex_init(&mutex, 0);
  pthread_cond_init(&cond, 0);
}

template<class T>
RingBuffer<T>::~RingBuffer() {
  delete[] buffer;
  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond);
}

template<class T>
T RingBuffer<T>::get() {
  T result;

  pthread_mutex_lock(&mutex);
  while(count == 0) {
    // Empty buffer, waiting for data.
    pthread_cond_wait(&cond, &mutex);
  }
  result = buffer[tail];
  --count;
  tail = (tail + 1) % _size;
  pthread_cond_signal(&cond);
  pthread_mutex_unlock(&mutex);
  
  return result;
}

template<class T>
void RingBuffer<T>::put(T value) {
  pthread_mutex_lock(&mutex);
  while(count == _size) {
    // Full buffer, waiting for space.
    pthread_cond_wait(&cond, &mutex);
  }
  buffer[head] = value;
  head = (head + 1) % _size;
  ++count;
  pthread_cond_signal(&cond);
  pthread_mutex_unlock(&mutex);
}

#endif
