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

#include <M2aiaCoreExports.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include <mitkExceptionMacro.h>
#include <numeric>
#include <vector>
#include <signal/m2SignalCommon.h>

namespace m2
{
  namespace Signal
  {
    /*!
     * TotalIonCurrent: This function approximates the area under the curve of an given spectrum.
     *
     * \param mIt0 GetXAxis start position
     * \param mItEnd GetXAxis end position
     * \param iIt0 Intensity array start (ImzMLImageSource of intensities is the same size as the source of masses)
     * \return tic Value of the TotalIonCount
     */
    template <class MzItFirst, class MzItLast, class IntItFirst>
    double TotalIonCurrent(MzItFirst mIt0, MzItLast mItEnd, IntItFirst iIt0) noexcept
    {
      double TIC = 0;
      auto mIt1 = std::next(mIt0);
      auto iIt1 = std::next(iIt0);
      for (; mIt1 != mItEnd; ++mIt0, ++mIt1, ++iIt0, ++iIt1)
      {
        TIC += ((*iIt0) + (*iIt1)) * 0.5 * ((*mIt1) - (*mIt0));
      }
      return TIC;
    }

    template <class IntItFirst, class IntItLast>
    double RootMeanSquare(IntItFirst iIt0, IntItLast iItEnd) noexcept
    {
      double N = std::distance(iIt0,iItEnd);
      double sum = 0;
      for (; iIt0 != iItEnd; ++iIt0)
        sum += *iIt0**iIt0;
      return std::sqrt(sum/N);
    }

    template <class ContainerType>
    double Median(ContainerType &ints) noexcept
    {
      return m2::Signal::Median(std::begin(ints), std::end(ints));
    }

    template <class ItFirst, class ItLast>
    double Median(ItFirst first, ItLast last) noexcept
    {
      const size_t n = std::distance(first, last);
      if ((n % 2) == 0)
      {
        std::nth_element(first, first + n / 2, last);
        std::nth_element(first, first + n / 2 + 1, last);
        return 0.5 * (*(first + n / 2) + *(first + (n / 2 + 1)));
      }
      else
      {
        std::nth_element(first, first + (n / 2), last);
        return *(first + n / 2);
      }
    }

  }; // namespace Signal
} // namespace m2
