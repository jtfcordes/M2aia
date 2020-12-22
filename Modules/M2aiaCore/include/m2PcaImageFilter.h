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

#include <MitkM2aiaCoreExports.h>
//#include <Eigen/Dense>
#include <itkBarrier.h>
#include <m2MassSpecVisualizationFilter.h>
#include <m2ImzMLMassSpecImage.h>
#include <mitkImage.h>
#include <mitkImageToImageFilter.h>
#include <vector>

namespace m2
{
  class MITKM2AIACORE_EXPORT PcaImageFilter : public m2::MassSpecVisualizationFilter
  {
  public:
    mitkClassMacro(PcaImageFilter, MassSpecVisualizationFilter);
    itkFactorylessNewMacro(Self);
    itkCloneMacro(Self);
    vnl_matrix<double> PcCalculation();
    void initMatrix();

  protected:
    //using MatrixXd = Eigen::Matrix<PixelType, Eigen::Dynamic, Eigen::Dynamic>;
    //MatrixXd m_DataMatrix;
    PcaImageFilter() = default;
    void GenerateData() override;

    /*This method is used to equalize the directions of the principal components. This allows the visualization of
     * m2PcaImageFilter, m2PcaVnlSvdImageFilter and m2PcaEigenImageFilter to be the same for all*/
    //void EqualizePcAxesEigen(MatrixXd *);
    //    ~m2PcaImageFilter() override;
  private:

  };
} // namespace m2
