# Extended Kalman Filter

My solution of CarND-Unscented-Kalman-Filter-Project assignment from Udacity Self Driving Car nanodegree course, Term 2. See project assignment starter code in https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project

---

## Dependencies

The project compilation and work have been verified under the following platforms: 
* Mac OX Sierra XCode 8.3.1 
* Windows 10 Visual Studio 2015 64-bit

I used cmake 3.7.2 to build project files.

## Code Style

To enforce [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html), I included Google's `cpplint.py` file available in `./src/Lint` folder. This tool expects installed python2.7. To check style of my code, run the following command line (from `./src/Lint`):

```
./cpplint.py ../*.h ../*.cpp
```

## Project Structure

The project consists of **UnscentedKF** and **UnscentedKFTests** applications, linking with **UnscentedKFLib** static library, which implements all the functionality.

**UnscentedKF** application has the same command line syntax as in the original [repo](https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project). Use the same build commands to generate and compile project files as described there.

For **UnscentedKFTests**, I used [GoogleTest](https://github.com/google/googletest) framework, cloned in `./lib` folder.

For more expressive and laconic code, I enforced using C++11 standard in CMake file by `-std=c++11` command.

The code consists of the following modules:
* `ukf.h`/`.cpp` implements Unscented Kalman Filter mathematics and CTRV car model with Radar and Lidar sensors
* `tools.h`/`.cpp` implements EvaluateRmse() function to calculate RMSE against data set
* `ground_truth_package.h`, `measurement_package.h` declares structures for ground truth and measurement data
* `main.cpp` is an entry point for **UnscentedKF** application
* `*_test.cpp` implement test cases for corresponding modules, run by **UnscentedKFTests** application