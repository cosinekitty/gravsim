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
#include "gravsim.h"


const vector_t ZeroVector = { {0.0, 0.0, 0.0} };


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


void Accelerations(sim_t *sim, int s)
{
    int i, j;
    vector_t rv;
    double r2, r3;

    /* Calculate the accelerations for the given state index. */
    /* Start out by initialzing all accelerations to 0. */
    for (i=0; i < sim->nbodies; ++i)
        sim->body[i].state[s].acc = ZeroVector;

    /* Explore every pair of distinct bodies. */
    for (i=0; i+1 < sim->nbodies; ++i)
    {
        double igm = sim->body[i].gm;
        state_t *si = &sim->body[i].state[s];
        for (j=i+1; j < sim->nbodies; ++j)
        {
            double jgm = sim->body[j].gm;
            state_t *sj = &sim->body[j].state[s];

            rv = Sub(si->pos, sj->pos);        /* vector from body[j] to body[i] */
            r2 = Dot(rv, rv);
            r3 = r2 * sqrt(r2);

            /* acceleration = GM / r^2 */
            /* divide by r^3 to also convert rv into a unit vector. */
            si->acc = Sub(si->acc, Mul(jgm/r3, rv));
            sj->acc = Add(sj->acc, Mul(igm/r3, rv));
        }
    }
}


void MoveBody(state_t *state, double dt)
{
    /*
        Adjust body state's position and velocity,
        using its current position, velocity, and acceleration.

        pos' = pos + vel*dt + (1/2)acc*dt^2
        vel' = vel + acc*dt
    */
    vector_t dv = Mul(dt, state->acc);
    vector_t dr = Add(Mul(dt, state->vel), Mul(dt/2.0, dv));
    state->vel = Add(state->vel, dv);
    state->pos = Add(state->pos, dr);
}


void MoveBodyColumn(sim_t *sim, int s, double dt)
{
    int b;

    for (b=0; b < sim->nbodies; ++b)
        MoveBody(&sim->body[b].state[s], dt);

    /* Update acceleration vectors to match new positions. */
    Accelerations(sim, s);
}



void SimUpdate1(sim_t *sim, double dt)
{
    /*
        This is a naive algorithm, just to get started.
        Apply the current state[0] accelerations as if they are
        constant over the interval dt.
    */
    MoveBodyColumn(sim,  0, dt);
    sim->tt += dt;
}
