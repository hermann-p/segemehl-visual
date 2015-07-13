/* 
 * A thread-safe FIFO queue for a one creator - one consumer scenario
 * Avoids use of lock/mutex structures to speed up access.
 * Pushes non-blocking, pops blocking.
 * When input is finished, ::done() must be called, else the last
 * pop will wait infinitely for more elements to appear.
 */


#ifndef LOCKFREEQUEUE
#define LOCKFREEQUEUE

#include <atomic>
#include <thread>

template <class T>
class LockFreeQueue {
private:
  struct Node {
    Node ( T val ) : value(val), next(nullptr) {}
    T value;
    Node* next;
  };

  Node* head;                      // marks first element
  std::atomic<Node*> divider;      // marks first not-consumed element
  std::atomic<Node*> tail;         // marks last element for appending nodes
  std::atomic<bool> done_;         // has appending finished?

public:
  LockFreeQueue ();
  ~LockFreeQueue ();
  void push ( const T& data );     // push one data element to queue, nonblocking
  bool pop ( T& dest );            // delete first element from queue and store it in dest
  bool isDone();                   // check if all elenemts have been pushed
  void done ( bool flag = true );  // set or reset done flag
};




template <class T>
LockFreeQueue<T>::LockFreeQueue()  {
  head = divider = tail = new Node( T() ); // create dummy
  done_ = false;
}

template <class T>
LockFreeQueue<T>::~LockFreeQueue () { // free remaining pointers
  while (head != nullptr) {
    Node* tmp = head->next;
    delete head;
    head = tmp;
  }
}

template <class T>
void LockFreeQueue<T>::push ( const T& data ) {
  Node* pTail = tail.load();      // get atomic state
  pTail->next = new Node(data);   // construct new node
  tail = pTail->next;             // update atomic state
  while (head != divider) {       // clean up data that has already been consumed...
    Node* tmp = head;
    head = head->next;
    delete tmp;
  }
}

template <class T>
bool LockFreeQueue<T>::pop ( T& dest ) {
  while (!done_ && divider == tail)
    std::this_thread::sleep_for(std::chrono::milliseconds(1));  // no elements, but elements expected: wait for element
  if (divider != tail) {           // if unconsumed elements exist...
    Node* div = divider.load();    // load atomic state
    dest = div->next->value;
    divider = div->next;           // advance atomic state
    return true; // success
  }
  return false;  // done was set and no more elements found - return false to end consumer loop
}

template <class T>
bool LockFreeQueue<T>::isDone() { return done_; }

template <class T>
void LockFreeQueue<T>::done ( bool flag ) { done_ = flag; }


#endif // LOCKFREEQUEUE
