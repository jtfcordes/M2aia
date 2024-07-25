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
#include <m2ISpectrumImageSource.h>
#include <m2SpectrumImage.h>
#include <signal/m2Baseline.h>
#include <signal/m2Smoothing.h>
#include <signal/m2Transformer.h>

namespace m2
{

  template <class MassAxisType, class IntensityType>
  class ImzMLSpectrumImageSource;

  /**
   * @class ImzMLSpectrumImage
   * @brief Represents an imzML spectrum image for interactive analysis in MSI applications.
   *
   * This class inherits from SpectrumImage that is a mitk::Image and provides functionality to work with imzML images.
   *
   *
   */
  class M2AIACORE_EXPORT ImzMLSpectrumImage final : public SpectrumImage
  {
  private:
    friend class Processor;

  public:
    typedef ImzMLSpectrumImage Self;
    typedef SpectrumImage Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    virtual std::vector<std::string> GetClassHierarchy() const override { return mitk::GetClassHierarchy<Self>(); }

    // Overwrite those methods to make MITK recognize this as default mitk::Image
    static const char *GetStaticNameOfClass() { return "Image"; }
    const char *GetNameOfClass() const override { return "Image"; }

    itkNewMacro(Self);

    itkGetConstReferenceMacro(ImzMLDataPath, std::string);
    itkSetMacro(ImzMLDataPath, std::string);

    itkGetConstReferenceMacro(BinaryDataPath, std::string);
    itkSetMacro(BinaryDataPath, std::string);

    
    /// @brief Type representing the offset from file start in bytes 
    using BinaryDataOffsetType = unsigned long long;
    
    /// @brief Type representing the data entry length in bytes
    using BinaryDataLengthType = unsigned long;
    
    /// @brief The BinarySpectrumMetaData structure holds meta data for a single spectrum. 
    struct BinarySpectrumMetaData
    {
      BinaryDataOffsetType mzOffset;
      BinaryDataOffsetType intOffset;
      BinaryDataLengthType mzLength;
      BinaryDataLengthType intLength;
      // size_t id;
      itk::Index<3> index;
      struct
      {
        float x, y, z;
      } world;
      // In imzML-file defined normalization constant
      m2::NormImagePixelType inFileNormalizationFactor = 1.0;
      // std::vector<m2::Peak> peaks;
    };

    /// @brief This type provides meta data objects to access spectra in the binary data file
    using SpectrumVectorType = std::vector<BinarySpectrumMetaData>;

    itkGetMacro(Spectra, SpectrumVectorType &);
    itkGetConstReferenceMacro(Spectra, SpectrumVectorType);

    void GetImage(double mz, double tol, const mitk::Image *mask, mitk::Image *img) const override;

    double GetXMin() const;
    double GetXMax() const;

    /**
     * @brief Initialize image access for processing.
     */
    void InitializeImageAccess() override;

    /**
     * @brief Initialize the image geometry.
     */
    void InitializeGeometry() override;

    /**
     * @brief Initialize the image processor.
     */
    void InitializeProcessor() override;

    /**
     *  @brief Initialize a normalization image for the given type
     */
    void InitializeNormalizationImage(m2::NormalizationStrategyType type) override;

    /**
     * @brief Get the spectrum as float values.
     * @param id The spectrum ID.
     * @param xs The vector to store x values.
     * @param ys The vector to store y values.
     * @param source The source index.
     */
    void GetSpectrumFloat(unsigned int id, std::vector<float> &xs, std::vector<float> &ys) const override;

    /**
     * @brief Get the intensities as float values.
     * @param id The spectrum ID.
     * @param xs The vector to store x values.
     * @param source The source index.
     */
    void GetIntensitiesFloat(unsigned int id, std::vector<float> &xs) const override;

    /**
     * @brief Get the spectrum as double values.
     * @param id The spectrum ID.
     * @param xs The vector to store x values.
     * @param ys The vector to store y values.
     * @param source The source index.
     */
    void GetSpectrum(unsigned int id, std::vector<double> &xs, std::vector<double> &ys) const override;

    /**
     * @brief Get the intensities as double values.
     * @param id The spectrum ID.
     * @param xs The vector to store x values.
     * @param source The source index.
     */
    void GetIntensities(unsigned int id, std::vector<double> &xs) const override;

    std::string GetMzGroupID() const {return m_MzGroupID;}
    std::string GetIntensityGroupID() const {return m_IntensityGroupID;}

  private:
    /// @brief path to the imzML file
    std::string m_ImzMLDataPath;

    /// @brief path to the ibd file
    std::string m_BinaryDataPath;

    /// @brief referenceableParameterGroupName
    std::string m_MzGroupID;
    
    /// @brief referenceableParameterGroupName
    std::string m_IntensityGroupID;

    /// @brief the source object is used to read/process data from the disk
    std::unique_ptr<m2::ISpectrumImageSource> m_SpectrumImageSource;

    /// @brief For each spectrum in the image exists a meta data object
    SpectrumVectorType m_Spectra;

    /// @brief Transformations are applied if available using elastix transformix
    std::vector<std::string> m_Transformations;

    mutable bool m_IsIonImageGenerated = true;

    ImzMLSpectrumImage();
    ~ImzMLSpectrumImage() override;
    using m2::SpectrumImage::InternalClone;
  };
} // namespace m2
