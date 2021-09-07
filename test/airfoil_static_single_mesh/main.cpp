#include "tailor.h" 

void dummy(Tailor::Tailor& tailor)
{
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(dummy);

    return 0;
}
