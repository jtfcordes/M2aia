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
#include <m2SpectrumImage.h>
#include <m2ElxRegistrationHelper.h>
#include <memory>

namespace m2
{
  /**
   * MS3DComposeDataObject represents a Composition of multiple mitk::MSDataObjects
   *
   **/
  class M2AIACORE_EXPORT SpectrumImageStack : public SpectrumImage
  {
  public:
    mitkClassMacro(SpectrumImageStack, SpectrumImage);
    mitkNewMacro2Param(Self, unsigned int, double);

  protected:
    
    /// @brief Deleted constructor without parameters
    SpectrumImageStack() = delete;

    /// @brief Constructor with one parameter
    /// @param stackSize Initial size of the stack
    /// @param spacingZ Image spacing along z axis
    explicit SpectrumImageStack(unsigned int stackSize, double spacingZ);

    /// @brief Type for storing slices of the stack. The m2::ElxRegistrationHelper is used
    /// to keep track of the referenced images
    using SpectrumImageStackVectorType = std::vector<std::shared_ptr<m2::ElxRegistrationHelper>>;

    /// @brief The vector holding images and transformations (elastix based image registration objects)
    SpectrumImageStackVectorType m_SliceTransformers;

    /// @brief Image data were transformed and inserted in the stack object at the given slice index
    /// @param warped Co-registered image
    /// @param stack Target stack image (a 3D volume)
    /// @param sliceIndex Index where the warped image will be added along the z-axis of the the 3D volume
    void CopyWarpedImageToStackImage(mitk::Image *warped, mitk::Image *stack, unsigned sliceIndex) const;


    unsigned int m_StackSize;
    double m_SpacingZ;

  public:
    // 
    void Insert(unsigned int sliceId ,std::shared_ptr<m2::ElxRegistrationHelper> transformer);

    const SpectrumImageStackVectorType & GetSliceTransformers() const {
      return m_SliceTransformers;
    }

    void InitializeGeometry() override;
    void InitializeProcessor() override;
    void InitializeImageAccess() override;
    void InitializeNormalizationImage(m2::NormalizationStrategyType /*type*/) override{}

    virtual void GetImage(double mz, double tol, const mitk::Image *mask, mitk::Image *img) const override;
    virtual void GetSpectrumFloat(unsigned int, std::vector<float> &, std::vector<float> &) const override{}
    virtual void GetIntensitiesFloat(unsigned int, std::vector<float> &) const override{}

    virtual void GetSpectrum(unsigned int, std::vector<double> &, std::vector<double> &) const override{}
    virtual void GetIntensities(unsigned int, std::vector<double> &) const override{}
  };

} // namespace m2
