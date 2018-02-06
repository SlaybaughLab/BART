#ifndef BART_TEST_HELPER_H_
#define BART_TEST_HELPER_H_

#include <string>

namespace btest {

class BartTestHelper {
 public:
  BartTestHelper() : gold_files_directory_("/test_data"),
                     report_(false) {};
  BartTestHelper(std::string gold_files_directory, bool report)
      : gold_files_directory_(gold_files_directory),
        report_(report) {};
 private:
  const std::string gold_files_directory_;
  const bool report_;
};

}

#endif
