#!/bin/bash
## Bash script to call runUnitTests for GitLab CI.
cd "$(dirname "$0")/unit"
matlab -nojvm -r "test_results_json = runUnitTests();" -logfile matlab.log
exit $?