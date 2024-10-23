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

#include <itkObject.h>
#include <M2aiaCoreExports.h>
#include <m2ImzMLSpectrumImage.h>
  
namespace m2
{

class M2AIACORE_EXPORT ShiftMapImageFilter : public itk::Object
{

private:

std::map<int, mitk::Image::Pointer> m_Outputs;
m2::ImzMLSpectrumImage::Pointer m_Input;

public:
mitkClassMacroItkParent(m2::ShiftMapImageFilter, itk::Object);
itkFactorylessNewMacro(m2::ShiftMapImageFilter);
void SetInput(m2::ImzMLSpectrumImage * image)
{
    m_Input = image;
}

mitk::Image::Pointer GetOutput(int idx)
{
    return m_Outputs[idx];
}

void GenerateData();

};

} 
