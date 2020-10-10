# Units
C++ Library for managing values with units

This standalone file (only Units.h is mandatory) allows you to use variables with units. The library is very useful for physical and chemical simulations. 

Very simple and intuitive to use, this library can save you a lot of trouble finding your errors in formulas. The error message are directly on the formulas in error. No need to check complicated template errors.

Units allow you to seamlessly convert unit from a valid type to another. Compile-time errors are thrown when units error occur in formulas.

Your comments are welcome to improve the library!

# Files
Units.h : Add units to value

LinearAlgebra.h : Define Points, Vectors and Matrices for use with units

Maths.h : Helper math functions that work with units

# Short example

```C++
#include <Units.h>
#include <iostream>

int main()
{
    Metre myMeter(5);
    Centimetre myCM(4.5);
    Decimetre dm = myMeter + myCM;

    Speed mySpeed = dm / Millisecond(4);
    std::cout << dm;
    return 0;
}
```

# To be added :
- Complete tutorials and documentation (but usage is pretty straightforward)

