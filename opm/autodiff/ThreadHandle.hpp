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

  protected:
    class EndObject : public ObjectInterface
    {
    public:
      void run () { }
      bool isEndMarker () const { return true; }
    };

    ////////////////////////////////////////////
    // class ThreadHandleObject
    ////////////////////////////////////////////
    class ThreadHandleObject
    {
    protected:
      std::queue< ObjectInterface* > objPtr_;
      std::mutex  mutex_;

      // no copying
      ThreadHandleObject( const ThreadHandleObject& ) = delete;

    public:
      // constructor creating object that is executed by thread
      ThreadHandleObject()
        : objPtr_(), mutex_()
      {
      }

      //! insert object into threads queue
      void push_back( ObjectInterface* obj )
      {
        // lock mutex to make sure objPtr is not used
        mutex_.lock();
        objPtr_.emplace( obj );
        mutex_.unlock();
      }

      // do the work
      void run()
      {
        // wait until objects have been pushed to the queue
        while( objPtr_.empty() )
        {
          // sleep one second
          std::this_thread::sleep_for( std::chrono::seconds(1) );
        }

        {
            // lock mutex for access to objPtr_
            mutex_.lock();

            // get next object from queue
            std::unique_ptr< ObjectInterface > obj( objPtr_.front() );
            objPtr_.pop();

            // unlock mutex for access to objPtr_
            mutex_.unlock();

            // if object is end marker terminate thread
            if( obj->isEndMarker() ){
                return;
            }

            // execute object action
            obj->run();
        }

        // keep thread running
        run();
      }
    }; // end ThreadHandleObject

    ////////////////////////////////////////////////////
    //  end ThreadHandleObject
    ////////////////////////////////////////////////////

    static void startThread( ThreadHandleObject* obj )
    {
       obj->run();
    }

    ThreadHandleObject threadObject_;
    std::thread thread_;

  private:
    // prohibit copying
    ThreadHandle( const ThreadHandle& ) = delete;

  public:
    // default constructor
    ThreadHandle()
      : threadObject_(),
        thread_( startThread, &threadObject_ )
    {
      // detach thread into nirvana
      thread_.detach();
    } // end constructor

    //! dispatch object to separate thread
    void dispatch( ObjectInterface* obj )
    {
      // add object to queue of objects
      threadObject_.push_back( obj ) ;
    }

    ~ThreadHandle()
    {
      // dispatch end object which will terminate the thread
      threadObject_.push_back( new EndObject() ) ;
    }
  };

} // end namespace Opm
#endif
