/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/m2aia/m2aia for details.

===================================================================*/

#include "m2SeriesDataProvider.h"

#include <mitkCoreServices.h>
#include <mitkIPreferencesService.h>
#include <mitkIPreferences.h>

#include <QColor>
#include <QHash>
#include <QLineSeries>
#include <QRandomGenerator>
#include <QScatterSeries>
#include <QVariant>
#include <QXYSeries>
#include <m2IntervalVector.h>
#include <mitkExceptionMacro.h>
#include <signal/m2Binning.h>

m2::SeriesDataProvider::SeriesDataProvider()
  : m_DefaultColor(QColor::fromHsv(QRandomGenerator::global()->generateDouble() * 255, 255, 255))
{
}

void m2::SeriesDataProvider::Initialize(const m2::IntervalVector *data, m2::SpectrumFormat format)
{
  if(format == m2::SpectrumFormat::None)
    m_Format = data->GetType();
  else
    m_Format = format;
  m_IntervalVector = data;

  InitializeSeries();
  Update();
  UpdateBoundaries();
}

m2::SeriesDataProvider::PointsVector m2::SeriesDataProvider::ConvertToQVector(const std::vector<double> &xs,
                                                                              const std::vector<double> &ys)
{
  PointsVector data;
  data.reserve(xs.size());
  auto inserter = std::back_inserter(data);
  std::transform(
    xs.begin(), xs.end(), ys.begin(), inserter, [](const auto &x, const auto &y) { return QPointF(x, y); });
  return data;
}

void m2::SeriesDataProvider::InitializeSeries()
{
  if (m_Format == m2::SpectrumFormat::Centroid)
  {
    m_Series = new QtCharts::QLineSeries();
    SetProfileSpectrumDefaultStyle(m_Series);
    m_Series->setPointsVisible(true);
  }
  else if (m_Format == m2::SpectrumFormat::Profile)
  {
    m_Series = new QtCharts::QLineSeries();
    SetProfileSpectrumDefaultStyle(m_Series);
    m_Series->setPointsVisible(false);
  }
}

void m2::SeriesDataProvider::Update()
{
  if (m_IntervalVector)
  {
    m_xs = m_IntervalVector->GetXMean();
    m_ys = m_IntervalVector->GetYMean(); // Critical?

    m_DataLoD.clear();
    if (m_Format == m2::SpectrumFormat::Centroid)
    {
      PointsVector target;
      using namespace std;
      transform(
        begin(m_xs), end(m_xs), begin(m_ys), back_inserter(target), [](auto a, auto b) { return QPointF(a, b); });
      m_DataLoD.push_back(target);
      return;
    }

    m_DataLoD.resize(m_Levels.size());
    unsigned int dataLodVectorIndex = 0;

    for (unsigned int level : m_Levels)
    {
      m_DataLoD[dataLodVectorIndex] = GenerateLoDData(m_xs, m_ys, level);
      ++dataLodVectorIndex;
    }
  }else{
    mitkThrow() << "Interval Vector not set correctly";
  }
}

void m2::SeriesDataProvider::UpdateBoundaries(double x1, double x2)
{
  using namespace std;

  if (m_xs.empty())
    return;

  if (x1 < 0 && x2 < 0)
  {
    x1 = m_xs.front();
    x2 = m_xs.back();
  }

  //   if (m_Format == Format::centorid)
  //   {
  //     m_Series->setVisible(true);

  //     auto lowerX = lower_bound(begin(m_xs), end(m_xs), x1);
  //     auto upperX = upper_bound(lowerX, end(m_xs), x2);

  //     auto lowerY = next(begin(m_ys), distance(begin(m_xs), lowerX));
  //     // auto upperY = next(lowerY, distance(begin(m_xs), upperX));

  //     PointsVector seriesData;
  //     seriesData.reserve(distance(lowerX, upperX));
  //     auto insert = back_inserter(seriesData);

  //     for (; lowerX != upperX; ++lowerX, ++lowerY)
  //       insert = QPointF{*lowerX, *lowerY};

  //     m_Series->replace(seriesData);
  //   }
  //   else
  {
    auto level = FindLoD(x1, x2);
    const auto &currentData = m_DataLoD[level];
    auto lower = std::lower_bound(std::begin(currentData),
                                  std::end(currentData),
                                  QPointF{x1, 0},
                                  [](const auto &a, const auto &b) { return a.x() < b.x(); });
    const auto upper = std::upper_bound(
      lower, std::end(currentData), QPointF{x2, 0}, [](const auto &a, const auto &b) { return a.x() < b.x(); });

    QVector<QPointF> seriesData;
    seriesData.reserve(std::distance(lower, upper));
    auto insert = std::back_inserter(seriesData);

    // marker highlighting
    unsigned int i = std::distance(std::begin(currentData), lower);

    for (; lower != upper; ++lower, ++i)
    {
      if (m_Format == m2::SpectrumFormat::Centroid)
      {
        insert = QPointF{lower->x(), -0.001};
        insert = QPointF{lower->x(), lower->y()};
        insert = QPointF{lower->x(), -0.001};
      }
      else
      {
        insert = QPointF{lower->x(), lower->y()};
      }
    }
    m_Series->replace(seriesData);

    using namespace QtCharts;
    // process series after points were added
  }
}

m2::SeriesDataProvider::PointsVector m2::SeriesDataProvider::GenerateLoDData(std::vector<double> &xs,
                                                                             std::vector<double> &ys,
                                                                             unsigned int level)
{
  if (xs.size() != ys.size())
    mitkThrow() << "Generate LODData faild: no equal xs size and ys size";

  const auto numBins = xs.size() / level;

  std::vector<m2::Interval> binnedData;
  m2::Signal::equidistantBinning(std::begin(xs), std::end(xs), std::begin(ys), std::back_inserter(binnedData), numBins);

  PointsVector target;
  target.reserve(binnedData.size());

  std::transform(std::begin(binnedData),
                 std::end(binnedData),
                 std::back_inserter(target),
                 [](const m2::Interval &p) {
                   return QPointF{p.x.mean(), p.y.max()};
                 });
  return target;
}

void m2::SeriesDataProvider::SetProfileSpectrumDefaultStyle(QtCharts::QXYSeries *series)
{
  auto p = series->pen();
  p.setWidthF(.75);
  p.setCapStyle(Qt::FlatCap);
  series->setPen(p);
}

void m2::SeriesDataProvider::SetCentroidSpectrumDefaultMarkerStyle(QtCharts::QXYSeries *series)
{
  if (auto scatterSeries = dynamic_cast<QtCharts::QScatterSeries *>(series))
  {
    scatterSeries->setMarkerShape(QtCharts::QScatterSeries::MarkerShape::MarkerShapeRectangle);
    scatterSeries->setMarkerSize(3);
  }
}

void m2::SeriesDataProvider::SetMarkerSpectrumDefaultMarkerStyle(QtCharts::QXYSeries *series)
{
  if (auto scatterSeries = dynamic_cast<QtCharts::QScatterSeries *>(series))
  {
    scatterSeries->setMarkerShape(QtCharts::QScatterSeries::MarkerShape::MarkerShapeCircle);
    scatterSeries->setMarkerSize(3);
    scatterSeries->setColor(Qt::GlobalColor::red);
  }
}

int m2::SeriesDataProvider::FindLoD(double xMin, double xMax) const
{
  if (m_Format == m2::SpectrumFormat::Centroid)
    return 0;

  auto* preferencesService = mitk::CoreServices::GetPreferencesService();
  auto* preferences = preferencesService->GetSystemPreferences();
  auto pointsWanted = preferences->GetInt("m2aia.view.spectrum.bins", 15000);

  int level, levelIndex = 0;
  const auto wantedDensity = pointsWanted / double(xMax - xMin);
  double delta = std::numeric_limits<double>::max();
  for (const auto &dataLodVector : m_DataLoD)
  {
    const auto lodDensity = dataLodVector.size() / double(dataLodVector.back().x() - dataLodVector.front().x());

    if (std::abs(wantedDensity - lodDensity) < delta)
    {
      level = levelIndex;
      delta = std::abs(wantedDensity - lodDensity);
    }
    ++levelIndex;
  }
  return level;
}

void m2::SeriesDataProvider::SetColor(QColor c)
{
  m_DefaultColor = c;
  if (m_Series)
  {
    m_Series->setColor(c);
  }
}

void m2::SeriesDataProvider::SetColor(qreal r, qreal g, qreal b, qreal a = 1.0)
{
  QColor c;
  c.setRgbF(r, g, b, a);
  SetColor(c);
}