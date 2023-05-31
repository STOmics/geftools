Coding Style Guide {#styleguide}
==========================

<!-- https://google.github.io/styleguide/cppguide.html#Variable_Names
https://zh-google-styleguide.readthedocs.io/en/latest/google-cpp-styleguide/contents/
https://github.com/opencv/opencv/wiki/Coding_Style_Guide -->

## Files
- All the file names are written in lower case for better compatibility with both POSIX and Windows.
- C++ interface headers have .h extension.
- Implementation files have .cpp extension.
- Documentation is written in .h files and any additional files (bibliographic reference, images) can be put into docs/ directory.
- Accuracy tests are put to test/ directory, performance tests - to perf/.

--------------

## File Structure
- Code lines should not be very long. Normally, they should be limited to 120 characters.
- No tabulation should be used. Set your editor to use spaces instead.
- Indentation is 4 spaces.

--------------

## Naming conventions
- Uses mixed-case style identifiers for external functions, types and class methods.
- Class names start with a capital letter.
- Methods and functions names start with a small letter.
- Macros and enumeration constants are written with all capital letters. Words are separated by underscore.
- The names of variables (including function parameters) and data members are all lowercase, with underscores between words. Data members of classes (but not structs) additionally have trailing underscores. For instance: a_local_variable, a_struct_data_member, a_class_data_member_.

--------------

## Designing functions and class interfaces

It is important to design function interface in a way, consistent with the rest of the library. The elements of function interface include:

- Functionality
- Name
- Return value
- Type of arguments
- Order of arguments
- Default values for some arguments

### Functionality

The functionality must be well defined and non-redundant. The function should be easily embedded into different processing pipelines that use other functions.

### Name

The name should basically reflect the function purpose.
There are a few common naming patterns:

- Majority of function names have form:  &lt;actionName&gt;&lt;Object&gt;&lt;Modifiers&gt; e.g. calibrateCamera, calcOpticalFlowPyrLK.

### Return value

It should be chosen to simplify function usage. Generally, a function that creates/computes a value should return it. It is the good practice to do so for the functions returning scalar values. However, for vector values this would lead to frequent allocation/deallocation of large memory blocks. Instead of creating and returning vector values, pass them as parameters (by reference).

Functions should not use return value for signaling about critical errors, such as null pointers, division by zero, bad argument range etc. Instead, they should throw an exception. On the other hand, it is recommended to use a return value to report normal run-time situations that can happen in a correctly working system.

### Types of arguments

A consistent argument order is important because it becomes easier to remember the order and it helps programmer to avoid errors, connecting with wrong argument order. The usual order is: input parameters, output parameters, flags and optional parameters.

Input parameters usually have const qualifiers. Large objects are normally passed by a constant reference; primitive types and small structures are passed by value.

Optional arguments often simplify function usage. Because C++ allows optional arguments in the end of parameters list only, it also may affect decisions on argument order—the most important flags go first and less important—after.

--------------

## Portability, External Dependencies
Code written for geftools (master branch) must comply with the C++11 standard.

One should get rid of compiler-dependent or platform-dependent constructions and system calls, such as:

- Compiler pragma's
- Specific keywords, e.g. __stdcall, __inline, __int64.
- Compiler extensions, e.g. special macros for min and max, overloaded macros etc.
- Inline assembly
- Unix or Win32-specific calls, e.g. bcopy, readdir, CreateFile, WaitForSingleObject etc.
- Concrete data sizes instead of sizeof's (sizeof(int) rather than 4), byte order ( *(int*)"\x1\x2\x3\x4" is 0x01020304 or 0x04030201 or what?), simple char instead of signed char or unsigned char anywhere except for text strings. Use short forms uchar for unsigned char and schar for signed char. Use preprocessor directives for surrounding non-portable pieces of code.

--------------
## Writing documentation on functions

The documentation for contributed functions is written using inline [Doxygen](https://www.doxygen.nl/manual) comments.

Use the existing documentation as an example. You are also welcome to provide tutorials for large descriptive chunks of text with pictures, code samples etc.

--------------
## Implementing tests

For tests we use GTest framework. Please, check the documentation at the prect site.
