============================
Adding to the Test Suite
============================

Ideally, new capabilities added to the code should be accompanied by a new
test in the test-suite that checks the proper execution of the specific new
addition.

To add a test to the suite, do the following

* in the ``tests/`` directory copy the ``_Template/`` directory to a new directory with a name appropriate to the test.

* Put the appropriate **param.lua** file, model file and any comparison files needed for running your test in that directory.

* Edit the **run_test.py** script so that it will plot up any relevant comparison(s) of the |sedona| output to a reference solution.

* Determine some metric of quantitative success of the test comparison. If the test results fail, set the **failure** variable to a non-zero value.


To verify that your new test is working correctly, do the following 

* Copy the current **sedona6.ex** executable into your new test directory.

* Run the code with ``./sedona6 param.lua``.

* Run the comparison script with ``python run_test.py``.

The script should produce your plots and return a SUCCESS or FAILURE when done.

Finally, to make your test part of the standard test suite

* Edit the **test_description.rst** file to describe the nature of your test.

* Add the name of your directory to the **suite_test_lists** file in the ``tests/`` directory.



