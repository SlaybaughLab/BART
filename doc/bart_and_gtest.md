# Integration of GTest with BART and use of xtrans_test

We are using GTest to develop unit testing for the BART code. Test
files should be in a `/tests/` directory in the folder of the code
that they are testing. In general, they should be named
`pronter_test.cc` if they are testing a class called `Pronter`. Any
`.cpp` or `.cc` file in a directory called `/tests` in the `/src`
directory will be identified as a file with tests, and will only be
compiled and linked with the test executable. Integration tests that
test functionality of multiple classes should be called
`directory_integration_test` where `directory` is the name of the
directory containing the classes that it is testing. Any mock classes
used by tests should be included in an `.h` file in the testing
directory where they are used. Finally, any gold files for the tests
(discussed below) should be in a `tests/data` directory.

For example, the file structure may look like this:

```
src
└── pronters
    ├── pronter.cc
    ├── pronter.h
    ├── scroncher.cc
    ├── scroncher.h
    └── tests
        ├── data
        │   └── pronter_test.gold
        ├── mock_scroncher.h
        ├── pronters_integration_test.cc
        └── pronter_test.cc
```

`pronter_test.cc` and `scroncher_test` provides _unit_ tests for the classes
`Pronter` and `Scroncher`, and `pronters_integration_test.cc` provides
_integration_ tests using multiple classes from the `pronters`
directory (both `Pronter`s and `Scroncher`s).


## Writing Tests

For information about writing unit tests, refer to the documentation
for
[google test](https://github.com/google/googletest/blob/master/googletest/docs/)
and
[google mock](https://github.com/google/googletest/tree/master/googlemock/docs).

## Writing Gold File Comparison Tests

Sometimes, unit testing is best accomplished by outputting the results
of a test routine to a text file and comparing to a known standard
text file, a *gold* file. The test will then only pass if both files
are identical. The process of developing a gold file comparison test
is outlined below.

The `BartTestHelper` class provides functions to use the
`dealii::deallog` functionality to accomplish this.

### Writing the test

A BART gold file test has three parts:

1. Initialize `dealii` logging.
2. Output test data to log.
3. Close log and compare to gold file.

A filename unique to the test (unique among all gold file tests) must
be chosen, and then the logging is initialized using:

```
std::string filename = "my_test";
btest::GoldTestInit(filename);
```

This will open the log, or throw a `runtime_error` if a log is already
open. Next, the test must output some data to the `dealii` log. For
example, we will output a test string using:

```
dealii::deallog << std::endl;
```

Finally, we can close the log and run the gold file test, using:

```
btest::GoldTestRun(filename);
```

All together, our test is the following:

```
std::string filename = "my_test";
btest::GoldTestInit(filename);

dealii::deallog << "test_me" << std::endl;

btest::GoldTestRun(filename);
```

Running `xtrans_test` will now generate the `my_test` file with our
logged string, and look for a file `test_data/my_test.gold` to compare
to. If the files are different or if there is no gold file to compare,
the test will fail.

### Generating a Report

By default if the test fails, the generated file is deleted, and no
information outside of the google test output is given. To see this
generated file, the `-r` or `--report` flag can be used:

```
./xtrans_test --report
```

This command will create a folder within `/test_data/` named with the
date and time the tests were run. If a test fails, the file generated
by the test will be moved here. Then if there was a valid gold file as
well, a `.diff` file will be generated showing the difference between
the gold file and the generated file (in unified format).

### Getting a Gold File

To get a gold file for new tests, follow this procedure, which follows
from the explained functionality above:

1. Run the new test with the `--report` flag, it will fail as there is
   no gold file.
2. Check the generated file in the report directory, this will be your
   new gold file.
2. Copy the generated file to the test data directory and add `.gold`,
   for example, `BART/test_data/my_test.gold`.
2. Run the test again and verify it passes.
2. Copy your new gold file to the `/tests/data` directory where your test is
   located.
