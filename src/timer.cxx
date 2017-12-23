//--------------------------------------------------------------------
//
//  Timer
//
//--------------------------------------------------------------------

#include <timer.h>

//--------------------------------------------------------------------
// Starts the timer.
void Timer::start (){
    //startTime = times( NULL );
    startTime = clock();
}

//--------------------------------------------------------------------
// Stops the timer.
void Timer::stop (){
    //stopTime = times( NULL );
    stopTime = clock();
}

//--------------------------------------------------------------------
// Computes the length of the time interval from startTime to
// stopTime (in seconds).
double Timer::get_time () const {
    return  double ( stopTime - startTime ) / (CLOCKS_PER_SEC) ;
}
