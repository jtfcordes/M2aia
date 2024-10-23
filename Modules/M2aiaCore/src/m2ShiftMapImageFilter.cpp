/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2ImzMLSpectrumImage.h>
#include <m2Process.hpp>
#include <m2ShiftMapImageFilter.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>

void m2::ShiftMapImageFilter::GenerateData(){

      m_Outputs[0] = static_cast<mitk::Image *>(m_Input.GetPointer())->Clone();
      m_Outputs[1] = mitk::Image::New();
      auto geom = m_Input->GetGeometry()->Clone();
      m_Outputs[1]->Initialize(mitk::MakeScalarPixelType<m2::ShiftImageType>(), *geom);
  
      auto centerProp = m_Input->GetProperty("m2aia.xs.selection.center");
      auto toleranceProp = m_Input->GetProperty("m2aia.xs.selection.tolerance");

      if (!centerProp || !toleranceProp)
        MITK_WARN << "Required properties (center, tolerance) not set.";

      auto centerFloatProp = dynamic_cast<mitk::DoubleProperty *>(centerProp.GetPointer());
      auto toleranceFloatProp = dynamic_cast<mitk::DoubleProperty *>(toleranceProp.GetPointer());
      double center = centerFloatProp->GetValue();
      double tolerance = toleranceFloatProp->GetValue();

      auto xAxis = m_Input->GetXAxis();
      auto lb = std::lower_bound(std::begin(xAxis), std::end(xAxis), center - tolerance);
      auto ub = std::lower_bound(std::begin(xAxis), std::end(xAxis), center + tolerance);

      auto ptrDiffLb = std::distance(std::begin(xAxis), lb);
      auto ptrDiffUb = std::distance(std::begin(xAxis), ub);

      auto meanSpectrum = m_Input->GetMeanSpectrum();
      auto meanMaxElement =
        std::max_element(std::begin(meanSpectrum) + ptrDiffLb, std::begin(meanSpectrum) + ptrDiffUb);
      auto meanMaxElementOffset = std::distance(std::begin(meanSpectrum), meanMaxElement);
      auto meanMzVal = xAxis[std::distance(std::begin(meanSpectrum), meanMaxElement)];


      mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(GetOutput(0));
      mitk::ImagePixelWriteAccessor<m2::ShiftImageType, 3> acc2(GetOutput(1));

      auto * d = GetOutput(1)->GetDimensions();
      auto N = std::accumulate(d, d+3, 1, std::multiplies<>());
      std::transform(acc2.GetData(),acc2.GetData()+N,acc2.GetData(),[](auto){return 0;});

      const auto &spectra = m_Input->GetSpectra();


      m2::Process::Map(m_Input->GetNumberOfValidPixels(),
                       20,
                       [&](auto /*thread id*/, auto a, auto b)
                       {
                         std::vector<double> ys, xs;
                         for (unsigned int i = a; i < b; ++i)
                         {
                          m_Input->GetSpectrum(i, xs, ys);
                          
                          auto lb = center < xs.front() ? xs.begin() : std::lower_bound(std::begin(xs), std::end(xs), center - tolerance);
                          auto ub = center > xs.back()? xs.end() : std::lower_bound(std::begin(xs), std::end(xs), center + tolerance);
                          auto ptrDiffLb = std::distance(std::begin(xs), lb);
                          auto ptrDiffUb = std::distance(std::begin(xs), ub);
                          auto maxElement = std::max_element(std::begin(ys) + ptrDiffLb, std::begin(ys) + ptrDiffUb);
                          auto maxElementOffset = std::distance(std::begin(ys), maxElement);
                         
                          if(*maxElement > *meanMaxElement*0.3){
                            auto mzVal = xAxis[std::distance(std::begin(ys), maxElement)];
                            acc.SetPixelByIndex(spectra[i].index, std::abs(mzVal - meanMzVal));
                            acc2.SetPixelByIndex(spectra[i].index, maxElementOffset - meanMaxElementOffset);
                          }else{
                            acc.SetPixelByIndex(spectra[i].index, 0);
                            acc2.SetPixelByIndex(spectra[i].index, 0);
                          }
                         }
                       });
         

 
    
}
