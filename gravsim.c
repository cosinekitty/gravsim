/*
    gravsim.c  -  by Don Cross

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

#include <math.h>
#include <stdio.h>
#include "gravsim.h"


const vector_t ZeroVector = { {0.0, 0.0, 0.0} };


void PrintVector(vector_t v)
{
    printf("(%0.8lg, %0.8lg, %0.8lg)", v.c[0], v.c[1], v.c[2]);
}


vector_t Vector(double x, double y, double z)
{
    vector_t vector;
    vector.c[0] = x;
    vector.c[1] = y;
    vector.c[2] = z;
    return vector;
}


vector_t Sub(vector_t a, vector_t b)
{
    vector_t c;
    c.c[0] = a.c[0] - b.c[0];
    c.c[1] = a.c[1] - b.c[1];
    c.c[2] = a.c[2] - b.c[2];
    return c;
}


vector_t Add(vector_t a, vector_t b)
{
    vector_t c;
    c.c[0] = a.c[0] + b.c[0];
    c.c[1] = a.c[1] + b.c[1];
    c.c[2] = a.c[2] + b.c[2];
    return c;
}


vector_t Average(vector_t a, vector_t b)
{
    vector_t c;
    c.c[0] = (a.c[0] + b.c[0]) / 2.0;
    c.c[1] = (a.c[1] + b.c[1]) / 2.0;
    c.c[2] = (a.c[2] + b.c[2]) / 2.0;
    return c;
}


vector_t Mul(double k, vector_t v)
{
    vector_t p;
    p.c[0] = k * v.c[0];
    p.c[1] = k * v.c[1];
    p.c[2] = k * v.c[2];
    return p;
}


double Dot(vector_t a, vector_t b)
{
    return a.c[0]*b.c[0] + a.c[1]*b.c[1] + a.c[2]*b.c[2];
}


void Accelerations(
    int nbodies,
    const body_t body[],
    const state_t state[],
    vector_t acc[])
{
    int i, j;
    vector_t rv;
    double r2, r3;

    for (i=0; i < nbodies; ++i)
        acc[i] = ZeroVector;

    /* Explore every pair of distinct bodies. */
    for (i=0; i+1 < nbodies; ++i)
    {
        for (j=i+1; j < nbodies; ++j)
        {
            rv = Sub(state[i].pos, state[j].pos);
            r2 = Dot(rv, rv);
            r3 = r2 * sqrt(r2);

            /* acceleration = GM / r^2 */
            /* divide by r^3 to also convert rv into a unit vector. */
            acc[i] = Sub(acc[i], Mul(body[j].gm/r3, rv));
            acc[j] = Add(acc[j], Mul(body[i].gm/r3, rv));
        }
    }
}


void MoveBody(const state_t *instate, state_t *outstate, vector_t acc, double dt)
{
    /*
        Adjust body state's position and velocity,
        using its current position, velocity, and acceleration.

        pos' = pos + vel*dt + (1/2)acc*dt^2
        vel' = vel + acc*dt
    */
    vector_t dv = Mul(dt, acc);
    vector_t dr = Add(Mul(dt, instate->vel), Mul(dt/2.0, dv));
    outstate->vel = Add(instate->vel, dv);
    outstate->pos = Add(instate->pos, dr);
}


void MoveAllBodies(int nbodies, const state_t instates[], state_t outstates[], vector_t acc[], double dt)
{
    int b;

    for (b=0; b < nbodies; ++b)
        MoveBody(&instates[b], &outstates[b], acc[b], dt);
}


void CopyStates(int nbodies, const state_t instates[], state_t outstates[])
{
    int b;

    for (b=0; b < nbodies; ++b)
        outstates[b] = instates[b];
}


double RelativeDiscrepancy(vector_t a, vector_t b)
{
    vector_t diff = Sub(a, b);
    double numer2 = Dot(diff, diff);
    double denom2 = Dot(a, a);
    return sqrt(numer2 / denom2);
}


void SimUpdate1(sim_t *sim, double dt)
{
    vector_t acc[MAX_BODIES];

    /*
        This is a naive algorithm, just to get started.
        Apply the current state[0] accelerations as if they are
        constant over the interval dt.
    */
    Accelerations(sim->nbodies, sim->body, sim->state, acc);
    MoveAllBodies(sim->nbodies, sim->state, sim->state, acc, dt);
    sim->tt += dt;
}


void SimUpdate2(sim_t *sim, double dt)
{
    int b, i;
    vector_t curr_acc[MAX_BODIES];
    vector_t next_acc[MAX_BODIES];
    vector_t mean_acc[MAX_BODIES];
    state_t  next_state[MAX_BODIES];

    Accelerations(sim->nbodies, sim->body, sim->state, curr_acc);
    MoveAllBodies(sim->nbodies, sim->state, next_state, curr_acc, dt);

    for (i=0; i < 3; ++i)
    {
        Accelerations(sim->nbodies, sim->body, next_state, next_acc);

        for (b=0; b < sim->nbodies; ++b)
            mean_acc[b] = Average(curr_acc[b], next_acc[b]);

        MoveAllBodies(sim->nbodies, sim->state, next_state, mean_acc, dt);
    }

    CopyStates(sim->nbodies, next_state, sim->state);
    sim->tt += dt;
}
