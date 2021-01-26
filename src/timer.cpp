#include "timer.h"

namespace Tailor
{
    Timer::Timer(): elapsed_seconds_(0.)
    {
        start();
    }

    time_point Timer::start_stop() const
    {
        return std::chrono::system_clock::now();
    }

    void Timer::start()
    {
        start_ = start_stop();
    }

    void Timer::stop()
    {
        end_ = start_stop();
        elapsed_seconds_ = elapsed();
    }

    double Timer::elapsed() const
    {
        duration d = end_ - start_;
        return d.count(); 
    }
}
