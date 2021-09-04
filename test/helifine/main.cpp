#include "tailor.h" 

void rotate(Tailor::Tailor& tailor)
{
    Tailor::Component compo;
    compo.read(Tailor::Tag(1));

    double rpm = compo.rpm_;
    double om = rpm * 2. * Tailor::PI / 60.; // rad/s
    double azimuth = om * tailor.solver()->dt(); // rad/s * time step
    int axis = compo.rotaxis_;
    Tailor::Vector3 pivot(compo.pivot_(0), compo.pivot_(1), compo.pivot_(2));

    tailor.rotate(Tailor::Tag(1), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(2), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(3), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(4), azimuth, axis, pivot);
    tailor.rotate(Tailor::Tag(5), azimuth, axis, pivot);
}

int main()
{
    Tailor::Tailor tailor;
    tailor.make(rotate);

    return 0;
}
