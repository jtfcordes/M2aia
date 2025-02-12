/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#ifndef m2LandmarkEvaluation_h
#define m2LandmarkEvaluation_h

#include <mitkPointSet.h>
#include <mitkImage.h>


namespace m2
{
  class LandMarkEvaluation
  {
  public:
	  struct Result {
		  double meanError, stdDevError;
		  std::vector<std::vector<double>> errorsPoints;
	  };

    static void PointSetToImage(const mitk::PointSet *input, mitk::Image *output);
    static void ImageToPointSet(const mitk::Image *input, mitk::PointSet *output);
	static Result EvaluatePointSet(const mitk::PointSet *input, const int & pointsPerSlice);

  protected:
    LandMarkEvaluation() = delete;
  };
} // namespace m2

#endif // m2LandmarkEvaluation_h
