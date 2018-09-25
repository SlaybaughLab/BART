find . -type f \( -iname \*.gcno -or -iname \*.gcda \) -exec cp {} . \;
lcov --directory . --capture --output-file coverage.info # capture coverage info
lcov --remove coverage.info '/usr/*' \
     '*_test.*' \
     '/home/dealii/BART/inc/*' \
     '/home/dealii/dealii.9/*' \
     '/home/dealii/libs/*' \
     --output-file coverage.info # filter out system and includes
lcov --list coverage.info #debug info
if [ $1 = "-r" ]
then
    genhtml coverage.info --branch-coverage --output-directory ./code_coverage_report/
fi
rm *.gcda *.gcno
