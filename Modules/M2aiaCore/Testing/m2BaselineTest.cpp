/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <cppunit/TestAssert.h>
#include <signal/m2Baseline.h>
#include <m2TestingConfig.h>
#include <m2TestFixture.h>
#include <mitkTestingMacros.h>

class m2BaselineTestSuite : public m2::TestFixture
{
  CPPUNIT_TEST_SUITE(m2BaselineTestSuite);
  MITK_TEST(TestTopHatBaseline);
  MITK_TEST(TestMedianBaseline);
  CPPUNIT_TEST_SUITE_END();

public:
  void TestTopHatBaseline()
  {
    m2::Signal::BaselineFunctor<double> bl;
    bl.Initialize(m2::BaselineCorrectionType::TopHat, 1);

    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> after = {0.0, 0.0, 0.0, 0.0, 0.0};  
    bl(data.begin(),data.end(), after.begin());

    for( auto v : after){
      std::cout << v << std::endl;
    }
    std::cout << std::endl;

    // CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, after, 0.0001);
  }

  void TestMedianBaseline()
  {
    m2::Signal::BaselineFunctor<double> bl;
    bl.Initialize(m2::BaselineCorrectionType::Median, 1);

    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> after = {0.0, 0.0, 0.0, 0.0, 0.0};
    bl(data.begin(),data.end(), after.begin());

    for( auto v : after){
      std::cout << v << std::endl;
    }
    std::cout << std::endl;

    // CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, after, 0.0001);
  }
};

MITK_TEST_SUITE_REGISTRATION(m2Baseline)