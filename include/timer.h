#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace Tailor
{
    typedef std::chrono::time_point<std::chrono::system_clock> time_point;
    typedef std::chrono::duration<double> duration;

    class Timer
    {
        time_point start_, end_;
        double elapsed_seconds_;

        public:

        Timer();

        time_point start_stop() const;
        void start();
        void stop();
        double elapsed() const;
    };
}

#endif
