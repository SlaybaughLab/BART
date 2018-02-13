#include "log_handler.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

namespace btest {

void LogHandler::OpenLogFile(std::unique_ptr<std::ofstream> log_stream) {
  log_stream_ = std::move(log_stream);
  dealii::deallog.attach(*log_stream_, false);
}

}
