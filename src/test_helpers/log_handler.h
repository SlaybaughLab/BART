#ifndef LOG_HANDLER_H_
#define LOG_HANDLER_H_

#include <memory>
#include <iostream>
#include <fstream>
#include <string>

#include "log_handler_i.h"

namespace btest {

class LogHandler : LogHandlerI{
 public:
  LogHandler() {};
  
  void OpenLogFile(std::unique_ptr<std::ofstream> log_stream);
  void CloseLogFile() { log_stream_.reset(); };
  bool IsLogging() const { return log_stream_ != nullptr; };
  
 private:
  std::unique_ptr<std::ofstream> log_stream_;
};

}

#endif
