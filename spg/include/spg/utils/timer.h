#pragma once

#include <chrono>

namespace spg
{
struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    inline void start()
    {
        m_running = true;
        m_start = Clock::now();
    }

    inline void stop()
    {
        m_stop = Clock::now();
        m_running = false;
    }

    inline float getSeconds() const { return getMilliseconds() / 1000.0f; }
    inline float getMilliseconds() const { return static_cast<float>(getMicroseconds()) / 1000.0f; }
    inline long long getMicroseconds() const
    {
        return m_running ? std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - m_start).count()
                         : std::chrono::duration_cast<std::chrono::microseconds>(m_stop - m_start).count();
    }

private:
    bool m_running{false};
    Clock::time_point m_start;
    Clock::time_point m_stop;
};
}  // namespace spg