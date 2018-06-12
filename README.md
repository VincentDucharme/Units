# Units
C++ Library for managing values with units


This standalone file allows you to use variables with units. The library is very useful for physical and chemical simulations. 

Very simple and intuitive to use, this library can save you a lot of trouble finding your errors in formulas.

Units allow you to seamlessly convert unit from a valid type to another. Compile-time errors are thrown when units error occur in formulas.

Your comments are welcome to improve the library!

# Short example

```C++
#include <Units.h>

int main()
{
    Metre myMeter(5);
    Centimetre myCM(4.5);
    Decimetre dm = myMeter + myCM;

    Speed mySpeed = dm / Millisecond(4);
    cout << dm;
    return 0;
}
```

# To be added :
- Complete tutorials and documentation

