/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <itkImageRegionIterator.h>
#include <itkVectorImage.h>
#include <m2TSNEImageFilter.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>

#include <m2SpectrumImage.h> // should be removed, only used for image types
#include <tsne/tsne.h>

m2::TSNEImageFilter::TSNEImageFilter() {}

m2::TSNEImageFilter::~TSNEImageFilter() {}

itk::DataObject::Pointer m2::TSNEImageFilter::MakeOutput(const DataObjectIdentifierType &name)
{
  itkDebugMacro("MakeOutput(" << name << ")");
  if (this->IsIndexedOutputName(name))
  {
    return this->MakeOutput(this->MakeIndexFromOutputName(name));
  }
  return mitk::Image::New().GetPointer();
}

itk::DataObject::Pointer m2::TSNEImageFilter::MakeOutput(DataObjectPointerArraySizeType /*idx*/)
{
  // auto input = this->GetInput(0);
  auto ret = mitk::Image::New();

  return ret.GetPointer();
}

void m2::TSNEImageFilter::SetNumberOfOutputDimensions(unsigned int v)
{
  this->SetNumberOfIndexedOutputs(v);
  this->m_NumberOfOutputDimensions = v;
}

void m2::TSNEImageFilter::SetProgressFunc(std::function<void()> f)
{
  this->m_ProgressFunc = f;
}

void m2::TSNEImageFilter::GenerateData()
{
  MITK_INFO << "-----------------------";
  MITK_INFO << "ImageTotSNEImage Filter";

  auto data = this->GetIndexedInputs();

  this->GetValidIndices();
  m_NumberOfValidPixels = m_ValidIndices.size();
  m_NumberOfOutputDimensions = 3;
  m_NumberOfInputDimensions = data.size();

  MITK_INFO << "Input features " << m_NumberOfInputDimensions << " ~ Numbers of indexed input images";
  MITK_INFO << "Output features " << m_NumberOfOutputDimensions;

  auto IndexedLoopOverMaskImage = [&](auto lambda) {
    for (auto const &index : m_ValidIndices)
    {
      lambda(index);
    }
  };

  auto LoopOverMaskImage = [&](auto lambda) {
    for (unsigned int i = 0; i < m_ValidIndices.size(); ++i)
    {
      lambda();
    }
  };

  unsigned rows = m_ValidIndices.size();
  unsigned cols = this->GetNumberOfInputs();

  MITK_INFO << "Number of datapoints " << m_NumberOfValidPixels;
  MITK_INFO << "Matrix size " << rows << "x" << cols;

  std::vector<double> d(m_NumberOfValidPixels * m_NumberOfInputDimensions);

  for (unsigned int c = 0; c < m_NumberOfInputDimensions; ++c)
  {
    // get image values
    unsigned j = 0;
    mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(data[c].GetPointer());
    mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> iia(image);
    //		double max = 0;
    double sum = 0;
    double sq_sum = 0;
    IndexedLoopOverMaskImage([&](const itk::Index<3> &index) { sum += iia.GetPixelByIndex(index); });

    double mean = sum / double(m_NumberOfValidPixels);
    IndexedLoopOverMaskImage([&](const itk::Index<3> &index) {
      auto v = iia.GetPixelByIndex(index);
      sq_sum += (v - mean) * (v - mean);
    });
    double stdev = std::sqrt((1.0 / double((m_NumberOfValidPixels - 1))) * sq_sum);

    IndexedLoopOverMaskImage([&](const itk::Index<3> &index) {
      d[j * m_NumberOfInputDimensions + c] = (iia.GetPixelByIndex(index) - mean) / stdev;
      ++j;
    });
  }

  MITK_INFO << "After the loop " << m_Perplexity;
  std::shared_ptr<TSNE::TSNE> tsne(new TSNE::TSNE);
  int max_iter = m_Iterations;
  double perplexity = m_Perplexity, theta = m_Theta;;
  int rand_seed = -1;
  // Now fire up the SNE implementation

  std::vector<double> Y(m_NumberOfValidPixels * m_NumberOfOutputDimensions, 0);
  MITK_INFO << "Perplexity " << m_Perplexity << " Iterations " << max_iter;
  tsne->run(d.data(),
            m_NumberOfValidPixels,
            m_NumberOfInputDimensions,
            Y.data(),
            m_NumberOfOutputDimensions,
            perplexity,
            theta,
            rand_seed,
            false,
            max_iter,
            250,
            250);
  MITK_INFO << "Finished";
  mitk::Image::Pointer input = this->GetInput(0);
  std::vector<double> firstDimension;
  firstDimension.reserve(m_NumberOfValidPixels);
  std::vector<double> scndDimension;
  scndDimension.reserve(m_NumberOfValidPixels);
  std::vector<double> thrdDimension;
  thrdDimension.reserve(m_NumberOfValidPixels);

  unsigned j = 0;
  auto output = mitk::Image::New();
  output->Initialize(input->GetPixelType(), *(input->GetGeometry()));
  mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> outAc(output);

  std::vector<std::vector<double>> values;

  LoopOverMaskImage([&]() {
    firstDimension.push_back(Y[j * m_NumberOfOutputDimensions + 0]);
    scndDimension.push_back(Y[j * m_NumberOfOutputDimensions + 1]);
    thrdDimension.push_back(Y[j * m_NumberOfOutputDimensions + 2]);
    ++j;
  });

  values.push_back(firstDimension);
  values.push_back(scndDimension);
  values.push_back(thrdDimension);

  std::vector<double> minValues;
  std::vector<double> maxValues;
  for (const std::vector<double> &vector : values)
  {
    double min = *std::min_element(vector.begin(), vector.end());
    double max = *std::max_element(vector.begin(), vector.end());

    minValues.push_back(min);
    maxValues.push_back(max);
  }

  // auto itkImage = itk::Image<RGBPixel, 3>::New();
  // initializeItkImage(itkImage);
  // itkImage->Allocate();

  auto vectorImage = initializeItkVectorImage(m_NumberOfComponents);
    unsigned int k = 0;
  for (auto const &indeex : m_ValidIndices)
  {
    itk::VariableLengthVector<m2::DisplayImagePixelType> variableVector;
    variableVector.SetSize(3);
    for (unsigned int i = 0; i < minValues.size(); ++i)
    {
      variableVector[i] = ((values[i][k] - minValues[i]) / (maxValues[i] - minValues[i])) * 255;
    }
    vectorImage->SetPixel(indeex, variableVector);
    ++k;
  }

  mitk::Image::Pointer outputImage = this->GetOutput();

  MITK_INFO << "Cast";
  mitk::CastToMitkImage(vectorImage, outputImage);
  MITK_INFO << "Cast to MITK Done";
  outputImage->SetSpacing(this->GetInput()->GetGeometry()->GetSpacing());
  outputImage->SetOrigin(this->GetInput()->GetGeometry()->GetOrigin());
}
