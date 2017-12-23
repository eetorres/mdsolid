//--------------------------------------------------------------------
//
//  Class declaration for the Timer
//--------------------------------------------------------------------

#ifndef _TIMER_H_
#define _TIMER_H_

// Insert a declaration for SystemTime here.
#include <iostream>
#include <ctime>
#include <time.h>

typedef clock_t SystemTime;

class Timer{

public:
    // Start and stop the timer
    void start ();
    void stop ();
    // Compute the elapsed time (in seconds)
    double get_time () const;

private:
    SystemTime startTime,   // Time that the timer was started
    stopTime;    // Time that the timer was stopped
};

#endif

