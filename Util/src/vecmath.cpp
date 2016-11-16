
#include "Util/vecmath.h"

int vec3f::random_number = 3457734;

void vec3f::random_init()
{
    random_number = 3457734;
}

double vec3f::random_double(){
    static int t = random_number;
    t = (t * 2345633 + t*1245 + t*356 + t*35 + t/34 + t/325 - 8647445);
    random_number = t;

    return ( double (abs(t)%10000)/10000);
}




vec3f vec3f::random(){
    static vec3f rnd;
    rnd.x=random_double()*2-1;
    rnd.y=random_double()*2-1;
    rnd.z=random_double()*2-1;
    rnd.normalize();
    return rnd;
}

vec3f vec3f::normalize( vec3f a )
{
    double square = a.x*a.x + a.y*a.y + a.z*a.z;
    if (square <= 0.00001f )
    {
        a.x=1;a.y=0;a.z=0;
        return a;
    }
    double len = 1.0f / (double)sqrt(square);
    a.x*=len;a.y*=len;a.z*=len;
    return a;
}
