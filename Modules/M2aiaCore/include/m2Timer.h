/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/
#pragma once

#include <chrono>
#include <iostream>
#include <mbilog.h>
#include <string>
#include <functional>
#include <mitkCommon.h>

namespace m2
{
  class Timer
  {
  public:
    using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
    using Duration = std::chrono::duration<float>;
    TimePoint start;
    Duration duration;
    std::string m_what;
    std::function<bool(Duration)> printIf = [](Duration){return true;};
    int m_N = 1;

    ~Timer()
    {
      auto now = std::chrono::high_resolution_clock::now();
      duration = now - start;
      if(printIf(duration))
        MITK_INFO << m_what << ": "<< duration.count()<< " Average(" << m_N <<" runs):" << duration.count()/double(m_N) << "s";
    }

    explicit Timer(std::string what) : m_what(what)
    {
      start = std::chrono::high_resolution_clock::now();
    }

    template<class F>
    explicit Timer(std::string what, int n, F && lambda) : Timer(what)
    {
      m_N = n;
      for(int i = 0; i < m_N; ++i)
        lambda();
    }
  };

} // namespace mitk