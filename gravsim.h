/*
    gravsim.h  -  by Don Cross

    Solar System gravity simulator.
    https://github.com/cosinekitty/gravsim

    MIT License

    Copyright (c) 2020 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
#ifndef __DDC_GRAVSIM_H
#define __DDC_GRAVSIM_H

#define CHECK(x)    do{if(0 != (error = (x))) goto fail;}while(0)
#define FAIL(...)   do{fprintf(stderr, __VA_ARGS__); error = 1; goto fail;}while(0)

#define SECONDS_PER_DAY           (24.0 * 3600.0)
#define AU_KM                     1.49597870691e+08
#define AU_M                      (AU_KM * 1000.0)
#define LIGHT_METERS_PER_SECOND   299792458.0
#define LIGHT_AU_PER_DAY          (LIGHT_METERS_PER_SECOND * (SECONDS_PER_DAY / AU_M))

#define MAX_BODIES   10

/*
    To use C in a functional programming style,
    I define the vector type 'vector_t' as an array inside a struct.
    Structs can be returned by value and passed by value,
    eliminating mutation and side-effects.
*/
typedef struct
{
    double c[3];        /* vector components: c[0]=x, c[1]=y, c[2]=z */
}
vector_t;


typedef struct
{
    vector_t    pos;    /* position vector [au] */
    vector_t    vel;    /* velocity vector [au/day] */
}
state_t;


typedef struct
{
    const char *name;       /* The name of the body, e.g. "Sun" or "Mars". */
    double      gm;         /* The product G*M (gravity constant times mass) for this body */
}
body_t;


typedef struct
{
    double   tt;                /* Terrestrial Time, relative to 1 January 2000 noon [days] */
    int      nbodies;
    body_t   body[MAX_BODIES];
    state_t  state[MAX_BODIES];
}
sim_t;

vector_t Vector(double x, double y, double z);
vector_t Sub(vector_t a, vector_t b);
vector_t Add(vector_t a, vector_t b);
vector_t Mul(double k, vector_t v);
double Dot(vector_t a, vector_t b);

double RelativeDiscrepancy(vector_t a, vector_t b);

void SimUpdate1(sim_t *sim, double dt);
void SimUpdate2(sim_t *sim, double dt);
void SimUpdate3(sim_t *sim, double dt);
void SimUpdate4(sim_t *sim, double dt);

#endif /* __DDC_GRAVSIM_H */
