#ifndef LOG_HANDLER_H_I_
#define LOG_HANDLER_H_I_

#include <memory>
#include <iostream>
#include <fstream>
#include <string>

namespace btest {

class LogHandlerI {
 public:
  ~LogHandlerI() {}
  
  virtual void OpenLogFile(std::unique_ptr<std::ofstream> log_stream) = 0;
  virtual void CloseLogFile() = 0;
  virtual bool IsLogging() const =0;
  
};

}

#endif
