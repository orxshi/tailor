#include "tailor.h" 

void move(Tailor::Tailor& tailor)
{
    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    Tailor::Vector3 vel(compo.u, compo.v, compo.w);
    Tailor::Vector3 distance;
    distance = tailor.solver()->dt() * vel;

    //tailor.move(Tailor::Tag(1), distance);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(&move);

    return 0;
}
