#include "tailor.h" 

void move(Tailor::Tailor& tailor)
{
    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    Vector3 v(compo.u, compo.v, compo.w);

    tailor.move(Tailor::Tag(1), v);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(move);

    return 0;
}
