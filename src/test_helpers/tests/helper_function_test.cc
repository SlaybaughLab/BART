#include "test_helpers/test_helper_functions.h"

#include <cstdlib>
#include <exception>
#include <unordered_map>
#include <utility>
#include <sstream>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::AllOf, ::testing::Ge, ::testing::Le;

class TestHelperFunctionTest : public ::testing::Test {
 protected:
};

TEST_F(TestHelperFunctionTest, RandomIntTest) {
  const std::vector<std::pair<int, int>> test_pairs = {
      {0, 100}, {-50, 50}, {-50, 0}, {-50, 100}, {0, 1},
      {1, 2}
  };

  for (auto [min, max] : test_pairs) {
    EXPECT_THAT(test_helpers::RandomInt(min, max), AllOf(Ge(min), Le(max)));
  }
  const std::vector<std::pair<double, double>> bad_pairs = {
      {100, 0}, {50, -50}, {0, -50}, {100, 100}, {0,0}, {-100, -100}};

  for (auto [min, max] : bad_pairs) {
    EXPECT_THROW(test_helpers::RandomInt(min, max), std::runtime_error)
              << "Values: min = " << min << ", max = " << max;
  }
}

TEST_F(TestHelperFunctionTest, RandomDoubleTest) {
  const std::vector<std::pair<double, double>> test_pairs = {
      {0, 100}, {-50, 50}, {-50, 0}, {0, 20.4}, {-102.4, 293}};

  for (auto p : test_pairs) {
    EXPECT_THAT(test_helpers::RandomDouble(p.first,p.second),
                ::testing::AllOf(::testing::Ge(p.first),
                                 ::testing::Le(p.second)));
  }

  const std::vector<std::pair<double, double>> bad_pairs = {
      {100, 0}, {50, -50}, {0, -50}, {100, 100}, {0,0}, {-100, -100}};

  for (auto p : bad_pairs) {

    EXPECT_THROW(test_helpers::RandomDouble(p.first,p.second),
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
          test_vector = test_helpers::RandomVector(vector_size, p.first, p.second);
      EXPECT_EQ(vector_size, test_vector.size());
      for (auto val : test_vector)
        EXPECT_THAT(val, ::testing::AllOf(::testing::Ge(p.first),
                                          ::testing::Le(p.second)));
      EXPECT_ANY_THROW(test_helpers::RandomVector(0, p.first, p.second));
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
            test_helpers::RandomIntVectorMap(map_size, vector_size, p.first, p.second);

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

        EXPECT_ANY_THROW(test_helpers::RandomIntVectorMap(0, vector_size, p.first,
                                                   p.second));
      }
    }
  }
}

TEST_F(TestHelperFunctionTest, RandomMatrixTest) {
  // Generate 10 random matrices and verify properties
  for (int i = 0; i < 10; ++i) {
    std::size_t n = std::rand()%10 + 1;
    std::size_t m = std::rand()%10 + 1;

    double val1 = -200 + static_cast<double>(std::rand()) / RAND_MAX*(400);
    double val2 = -200 + static_cast<double>(std::rand()) / RAND_MAX*(400);
    double min_value = (val1 > val2) ? val2 : val1;
    double max_value = (val1 > val2) ? val1 : val2;

    std::ostringstream test_msg;
    test_msg << "Matrix Size: " << m << "x" << n << "; min/max: "
             << min_value << "/" << max_value;

    dealii::FullMatrix<double> matrix = test_helpers::RandomMatrix(m, n, min_value,
                                                            max_value);

    for (auto cit = matrix.begin(); cit < matrix.end(); ++cit)
      EXPECT_THAT((*cit).value(), ::testing::AllOf(::testing::Ge(min_value),
                                                   ::testing::Le(max_value)))
                << test_msg.str();

    EXPECT_EQ(matrix.m(), m) << test_msg.str();
    EXPECT_EQ(matrix.n(), n) << test_msg.str();
  }
}

TEST_F(TestHelperFunctionTest, RandomIntMatrixMapTest) {
  // Generate 20 random int to matrix maps and verify properties
  for (int i = 0; i < 20; ++i) {
    std::size_t n = std::rand()%10 + 1;
    std::size_t m = std::rand()%10 + 1;
    std::size_t map_size = std::rand()%10 + 1;

    double val1 = -200 + static_cast<double>(std::rand()) / RAND_MAX*(400);
    double val2 = -200 + static_cast<double>(std::rand()) / RAND_MAX*(400);
    double min_value = (val1 > val2) ? val2 : val1;
    double max_value = (val1 > val2) ? val1 : val2;

    std::ostringstream test_msg;
    test_msg << "Map size: " << map_size << "; Matrix Size: " << m << "x"
             << n << "; min/max: " << min_value << "/" << max_value;

    auto matrix_map = test_helpers::RandomIntMatrixMap(map_size, m, n, min_value,
                                                max_value);

    EXPECT_EQ(matrix_map.size(), map_size) << test_msg.str();

    for (const auto mat : matrix_map) {
      EXPECT_THAT(mat.first, ::testing::AllOf(::testing::Ge(0),
                                              ::testing::Le(map_size*10)))
                << test_msg.str();

      for (auto cit = mat.second.begin(); cit < mat.second.end(); ++cit)
        EXPECT_THAT((*cit).value(), ::testing::AllOf(::testing::Ge(min_value),
                                                     ::testing::Le(max_value)))
                  << test_msg.str();

      EXPECT_EQ(mat.second.m(), m) << test_msg.str();
      EXPECT_EQ(mat.second.n(), n) << test_msg.str();
    }
  }
}

} // namespace


