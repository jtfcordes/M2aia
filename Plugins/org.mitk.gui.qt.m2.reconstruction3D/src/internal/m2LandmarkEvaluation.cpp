/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

// Blueberry
#include "m2LandmarkEvaluation.h"

#include <itkImageIteratorWithIndex.h>
#include <mitkImageAccessByItk.h>
#include <mitkImagePixelReadAccessor.h>
#include <numeric>

void m2::LandMarkEvaluation::PointSetToImage(const mitk::PointSet *input, mitk::Image *output)
{
  if (!input)
  {
    MITK_ERROR << "Input PointSet is null!";
    return;
  }

  if (!output)
  {
    MITK_ERROR << "Output Image is null!";
    return;
  }

  unsigned id = 1;

  AccessByItk(const_cast<mitk::Image *>(output), ([&](auto I) {
                using ImageType = itk::Image<typename std::remove_pointer<decltype(I)>::type::PixelType,
                                             std::remove_pointer<decltype(I)>::type::ImageDimension>;
                I->FillBuffer(0);

                for (auto it = input->Begin(); it != input->End(); ++it)
                {
                  typename ImageType::IndexType idx;
                  output->GetGeometry()->WorldToIndex(it->Value(), idx);
                  I->SetPixel(idx - typename ImageType::OffsetType{1, 0}, id);
                  I->SetPixel(idx + typename ImageType::OffsetType{1, 0}, id);
                  I->SetPixel(idx - typename ImageType::OffsetType{0, 1}, id);
                  I->SetPixel(idx + typename ImageType::OffsetType{0, 1}, id);
                  I->SetPixel(idx, id);
                  ++id;
                }
              }));
}

void m2::LandMarkEvaluation::ImageToPointSet(const mitk::Image *input, mitk::PointSet *output)
{
  if (!input)
  {
    MITK_ERROR << "Input Image is null!";
    return;
  }

  if (!output)
  {
    MITK_ERROR << "Output PointSet is null!";
    return;
  }

  output->SetGeometry(input->GetGeometry());

  if (input != nullptr && output != nullptr)
  {
    AccessByItk(const_cast<mitk::Image *>(input), ([&](auto I) {
                  using ImageType = itk::Image<typename std::remove_pointer<decltype(I)>::type::PixelType,
                                               std::remove_pointer<decltype(I)>::type::ImageDimension>;
                  ImageType *itkImage = I;
                  itk::ImageRegionConstIteratorWithIndex<ImageType> it(itkImage, itkImage->GetLargestPossibleRegion());
                  mitk::Point3D p;
                  std::map<unsigned, std::map<unsigned, std::vector<vnl_vector<double>>>> layers;
                  while (!it.IsAtEnd())
                  {
                    if (it.Value() != 0)
                    {
                      input->GetGeometry()->IndexToWorld(it.GetIndex(), p);
                      layers[it.GetIndex()[2]][it.Value() - 1].push_back(p.GetVnlVector());
                    }
                    ++it;
                  }

                  unsigned maxPointsInLayer = 0;
                  for (auto &layer : layers)
                  {
                    if (layer.second.size() > maxPointsInLayer)
                      maxPointsInLayer = layer.second.size();
                  }
                  for (auto &layer : layers)
                  {
                    const auto &layerId = layer.first;
                    const auto &layerPoints = layer.second;
                    for (auto &kv : layerPoints)
                    {
                      const auto &id = kv.first;
                      const auto &indices = kv.second;
                      vnl_vector<double> s(3, 0);
                      for (auto p : indices)
                        s += p;
                      s /= indices.size();
                      output->InsertPoint(id + layerId * maxPointsInLayer, mitk::Point3D(s.data_block()));
                    }
                  }
                }));
  }
}

m2::LandMarkEvaluation::Result m2::LandMarkEvaluation::EvaluatePointSet(const mitk::PointSet *input,
                                                                        const int &pointsPerSlice)
{
  Result r;
  r.errorsPoints.resize(pointsPerSlice);
  MITK_INFO << "Evaluate " << input->GetSize() << " Points";
  for (int i = 0; i < input->GetSize() - pointsPerSlice; ++i)
  {
    auto c0 = input->GetPoint((i + pointsPerSlice)); // point next layer
    auto c1 = input->GetPoint(i);                    // point current layer
    MITK_INFO << c0 << " " << c1 << "\n";
    c0[2] = c1[2] = 0; // look only in x,y space
    r.errorsPoints[i % pointsPerSlice].push_back((c0.GetVnlVector() - c1.GetVnlVector()).magnitude());
  }

  std::vector<double> errors;
  for (auto &m : r.errorsPoints)
  {
    auto val = std::accumulate(std::begin(m), std::end(m), 0.0) / float(m.size());
    errors.emplace_back(val);
  }

  const auto meanError = std::accumulate(std::begin(errors), std::end(errors), 0.0) / float(errors.size());
  auto stdDevError =
    std::sqrt(std::accumulate(std::begin(errors),
                              std::end(errors),
                              0.0,
                              [&meanError](auto s, auto x) { return s + (x - meanError) * (x - meanError); }) /
              float(errors.size()));
  MITK_INFO << meanError << " +- " << stdDevError << "\n";
  r.meanError = meanError;
  r.stdDevError = stdDevError;

  return r;
}
