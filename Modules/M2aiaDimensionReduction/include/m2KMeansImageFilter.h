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
#include <M2aiaDimensionReductionExports.h>
#include <m2IntervalVector.h>
#include <mitkLabelSetImage.h>
#include <m2ImzMLSpectrumImage.h>


// Eigen
// #include "Eigen/Dense"

namespace m2 {


class M2AIADIMENSIONREDUCTION_EXPORT KMeansImageFilter : public itk::Object
{

private:

    std::map<int, mitk::Image::Pointer> m_Outputs;
    m2::ImzMLSpectrumImage::Pointer m_Input;

public:
    mitkClassMacroItkParent(KMeansImageFilter, itk::Object);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);
    typedef mitk::LabelSetImage OutputType;

    itkSetMacro(NumberOfClusters, unsigned int);
    void GenerateData();
    
    void SetInput(m2::ImzMLSpectrumImage * image)
    {
        m_Input = image;
    }
    
    void SetIntervals(std::vector<m2::Interval> intervals){
        m_Intervals = intervals;
    }
    

    mitk::LabelSetImage::Pointer GetOutput(int idx)
    {
        if(m_Outputs.find(idx) == m_Outputs.end())
        {
            m_Outputs[idx] = mitk::LabelSetImage::New();
        }
        return dynamic_cast<mitk::LabelSetImage*>(m_Outputs[idx].GetPointer());
         
    }
 
private:

    void DoKMeans(const Eigen::MatrixXd& data, int k, std::vector<int>& clusterAssignments);

    std::vector<m2::Interval> m_Intervals;
    unsigned int m_NumberOfClusters = 0;
};

} // end namespace m2

