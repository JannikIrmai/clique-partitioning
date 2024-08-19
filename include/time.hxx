#include <chrono>

using Time = std::chrono::steady_clock;
using ms = std::chrono::milliseconds;
using float_sec = std::chrono::duration<float>;
using float_time_point = std::chrono::time_point<Time, float_sec>;