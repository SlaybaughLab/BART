#include <gtest/gtest.h>

#include <cstdlib>
#include <exception>
#include <unordered_map>
#include <utility>
#include <sstream>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "../gmock_wrapper.h"

// Forward declarations

namespace btest {
double RandomDouble(double, double);
std::vector<double> RandomVector(size_t, double, double);
std::unordered_map<int, std::vector<double>> RandomIntVectorMap(size_t, size_t,
                                                                double,double);
dealii::FullMatrix<double> RandomMatrix(size_t, size_t, double, double);
}

class TestHelperFunctionTest : public ::testing::Test {
 protected:
};

TEST_F(TestHelperFunctionTest, RandomDoubleTest) {

  const std::vector<std::pair<double, double>> test_pairs = {
    {0, 100}, {-50, 50}, {-50, 0}, {0, 20.4}, {-102.4, 293}};

  for (auto p : test_pairs) {
    EXPECT_THAT(btest::RandomDouble(p.first,p.second),
                ::testing::AllOf(::testing::Ge(p.first),
                                 ::testing::Le(p.second)));
  }

  const std::vector<std::pair<double, double>> bad_pairs = {
    {100, 0}, {50, -50}, {0, -50}, {100, 100}, {0,0}, {-100, -100}};

  for (auto p : bad_pairs) {
    
    EXPECT_THROW(btest::RandomDouble(p.first,p.second),
                 std::runtime_error) << "Values: min = " << p.first << ", max = "
                                     << p.second;
  }
}

TEST_F(TestHelperFunctionTest, RandomVectorTest) {
  const std::vector<std::pair<double, double>> test_pairs = {
    {0, 100}, {-50, 50}, {-50, 0}, {0, 20.4}, {-102.4, 293}};
  const std::vector<int> test_vector_length = {1, 4, 5, 100};

  for (const auto vector_size : test_vector_length) {
    for (const auto p : test_pairs) {
      std::vector<double>
          test_vector = btest::RandomVector(vector_size, p.first, p.second);
      EXPECT_EQ(vector_size, test_vector.size());
      for (auto val : test_vector)
        EXPECT_THAT(val, ::testing::AllOf(::testing::Ge(p.first),
                                          ::testing::Le(p.second)));
      EXPECT_ANY_THROW(btest::RandomVector(0, p.first, p.second));
    }
  }
}

TEST_F(TestHelperFunctionTest, RandomIntVectorMapTest) {
  const std::vector<std::pair<double, double>> test_pairs = {
    {0, 100}, {-50, 50}, {-50, 0}, {0, 20.4}, {-102.4, 293}};
  const std::vector<int> test_vector_length = {1, 4, 5, 100};
  const std::vector<int> test_map_size = {1, 4, 5, 8};

  for (const auto map_size : test_map_size) {
    for (const auto vector_size : test_vector_length) {
      for (const auto p : test_pairs) {
        std::unordered_map<int, std::vector<double>> test_map =
            btest::RandomIntVectorMap(map_size, vector_size, p.first, p.second);
        
        std::ostringstream test_msg;
        test_msg << "Map Size: " << map_size << "; vector size: "
                  << vector_size << "; min/max: "
                  << p.first << "/" << p.second;                                
        
        EXPECT_EQ(test_map.size(), map_size) << test_msg.str(); // Correct number of materials
        
        for (const auto mat : test_map) {
          // Material ID is between min and max
          EXPECT_THAT(mat.first, ::testing::AllOf(::testing::Ge(0),
                                                  ::testing::Le(map_size*10)))
              << test_msg.str();
          // Vector has correct values
          for (auto val : mat.second)
            EXPECT_THAT(val, ::testing::AllOf(::testing::Ge(p.first),
                                              ::testing::Le(p.second)))
                << test_msg.str();
        }

        EXPECT_ANY_THROW(btest::RandomIntVectorMap(0, vector_size, p.first,
                                                   p.second));
      }
    }
  }
}

TEST_F(TestHelperFunctionTest, RandomMatrixTest) {
  // Generate 10 random matrices and verify properties
  for (int i = 0; i < 10; ++i) {
    size_t n = rand()%10 + 1;
    size_t m = rand()%10 + 1;
    
    double val1 = -200 + static_cast<double>(rand()) / RAND_MAX*(400);
    double val2 = -200 + static_cast<double>(rand()) / RAND_MAX*(400);
    double min = (val1 > val2) ? val2 : val1;
    double max = (val1 > val2) ? val1 : val2;

    std::ostringstream test_msg;
    test_msg << "Matrix Size: " << m << "x" << n << "; min/max: "
             << min << "/" << max;                                
    
    dealii::FullMatrix<double> matrix = btest::RandomMatrix(m, n, min, max);

    for (auto cit = matrix.begin(); cit < matrix.end(); ++cit)
      EXPECT_THAT((*cit).value(), ::testing::AllOf(::testing::Ge(min),
                                                   ::testing::Le(max)))
          << test_msg.str();

    EXPECT_EQ(matrix.m(), m) << test_msg.str();
    EXPECT_EQ(matrix.n(), n) << test_msg.str();
  }
}
