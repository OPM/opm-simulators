#ifndef OPM_THREADHANDLE_HPP
#define OPM_THREADHANDLE_HPP

#include <cassert>
#include <dune/common/exceptions.hh>

#include <thread>
#include <mutex>
#include <queue>

namespace Opm
{

  class ThreadHandle
  {
  public:
    class ObjectInterface
    {
    protected:
      ObjectInterface() {}
    public:
      virtual ~ObjectInterface() {}
      virtual void run() = 0;
      virtual bool isEndMarker () const { return false; }
    };

    template <class Object>
    class ObjectWrapper : public ObjectInterface
    {
      Object obj_;
    public:
      ObjectWrapper( Object&& obj ) : obj_( std::move( obj ) ) {}
      void run() { obj_.run(); }
    };

  protected:
    class EndObject : public ObjectInterface
    {
    public:
      void run () { }
      bool isEndMarker () const { return true; }
    };

    ////////////////////////////////////////////
    // class ThreadHandleQueue
    ////////////////////////////////////////////
    class ThreadHandleQueue
    {
    protected:
      std::queue< std::unique_ptr< ObjectInterface > > objQueue_;
      std::mutex  mutex_;

      // no copying
      ThreadHandleQueue( const ThreadHandleQueue& ) = delete;

    public:
      //! constructor creating object that is executed by thread
      ThreadHandleQueue()
        : objQueue_(), mutex_()
      {
      }

      //! insert object into threads queue
      void push_back( std::unique_ptr< ObjectInterface >&& obj )
      {
        // lock mutex to make sure objPtr is not used
        mutex_.lock();
        objQueue_.emplace( std::move(obj) );
        mutex_.unlock();
      }

      //! do the work until the queue received an end object
      void run()
      {
        // wait until objects have been pushed to the queue
        while( objQueue_.empty() )
        {
          // sleep one second
          std::this_thread::sleep_for( std::chrono::seconds(1) );
        }

        {
            // lock mutex for access to objQueue_
            mutex_.lock();

            // get next object from queue
            std::unique_ptr< ObjectInterface > obj( objQueue_.front().release() );
            // remove object from queue
            objQueue_.pop();

            // unlock mutex for access to objQueue_
            mutex_.unlock();

            // if object is end marker terminate thread
            if( obj->isEndMarker() ){
                if( ! objQueue_.empty() ) {
                    OPM_THROW(std::logic_error,"ThreadHandleQueue: not all queued objects were executed");
                }
                return;
            }

            // execute object action
            obj->run();
        }

        // keep thread running
        run();
      }
    }; // end ThreadHandleQueue

    ////////////////////////////////////////////////////
    //  end ThreadHandleQueue
    ////////////////////////////////////////////////////

    // start the thread by calling method run
    static void startThread( ThreadHandleQueue* obj )
    {
       obj->run();
    }

    ThreadHandleQueue threadObjectQueue_;
    std::thread thread_;

  private:
    // prohibit copying
    ThreadHandle( const ThreadHandle& ) = delete;

  public:
    //! default constructor
    ThreadHandle()
      : threadObjectQueue_(),
        thread_( startThread, &threadObjectQueue_ )
    {
      // detach thread into nirvana
      thread_.detach();
    } // end constructor

    //! dispatch object to queue of separate thread
    template <class Object>
    void dispatch( Object&& obj )
    {
      typedef ObjectWrapper< Object >  ObjectPointer;
      ObjectInterface* objPtr = new ObjectPointer( std::move(obj) );

      // add object to queue of objects
      threadObjectQueue_.push_back( std::unique_ptr< ObjectInterface > (objPtr) );
    }

    //! destructor terminating the thread
    ~ThreadHandle()
    {
      // dispatch end object which will terminate the thread
      threadObjectQueue_.push_back( std::unique_ptr< ObjectInterface > (new EndObject()) ) ;
    }
  };

} // end namespace Opm
#endif
