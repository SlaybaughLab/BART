lcov --directory . --capture --output-file coverage.info # capture coverage info
lcov --remove coverage.info '/usr/*' --output-file coverage.info # filter out system
lcov --remove coverage.info '*_test.*' --output-file coverage.info
lcov --list coverage.info #debug info
bash <(curl -s https://codecov.io/bash)
