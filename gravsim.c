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


double PosError(int nbodies, state_t state1[], state_t state2[])
{
    int b;
    double sum;

    sum = 0.0;
    for (b = 0; b < nbodies; ++b)
        sum += RelativeDiscrepancy(state1[b].pos, state2[b].pos);

    return sum;
}


void SimUpdate1(sim_t *sim, double dt)
{
    vector_t acc[MAX_BODIES];

    /* Calculate the accelerations acting on the bodies at the current time. */
    Accelerations(sim->nbodies, sim->body, sim->state, acc);

    /* Naively assume that accerlation applies over the entire time increment. */
    MoveAllBodies(sim->nbodies, sim->state, sim->state, acc, dt);
    sim->tt += dt;
}


void ApproximateMovement(
    int nbodies,
    double dt,
    const body_t body[],
    const state_t curr_state[],
    state_t next_state[],
    vector_t curr_acc[],
    vector_t mean_acc[],
    vector_t next_acc[])
{
    int i, b;

    /* Calculate accelerations of each body at the current time. */
    Accelerations(nbodies, body, curr_state, curr_acc);

    /* Move the bodies as if current accerlation applies over the whole interval [0, dt]. */
    MoveAllBodies(nbodies, curr_state, next_state, curr_acc, dt);

    for (i = 0; i < 2; ++i)
    {
        /* Calculate accelerations of the estimated next location of the bodies. */
        Accelerations(nbodies, body, next_state, next_acc);

        /* Take the average of the beginning and ending accelerations */
        /* as estimates for mean acceleration. */
        for (b = 0; b < nbodies; ++b)
            mean_acc[b] = Average(curr_acc[b], next_acc[b]);

        /* Refine the estimate of where the bodies will be after dt. */
        MoveAllBodies(nbodies, curr_state, next_state, mean_acc, dt);
    }
}


void SimUpdate2(sim_t *sim, double dt)
{
    state_t next_state[MAX_BODIES];
    vector_t curr_acc[MAX_BODIES];
    vector_t mean_acc[MAX_BODIES];
    vector_t next_acc[MAX_BODIES];

    /* Find a time-reversible mean acceleration over the interval dt. */
    ApproximateMovement(sim->nbodies, dt, sim->body, sim->state, next_state, curr_acc, mean_acc, next_acc);

    /* Update the current state of each body to be the final refined estimate. */
    CopyStates(sim->nbodies, next_state, sim->state);
    sim->tt += dt;
}


void SimUpdate3(sim_t *sim, double dt)
{
    int b, k;
    double J, K, L, A, B, E, F, G, p;
    double dt2, dt3, dt4;
    state_t next_state[MAX_BODIES];
    state_t middle_state[MAX_BODIES];
    vector_t curr_acc[MAX_BODIES];
    vector_t middle_acc[MAX_BODIES];
    vector_t next_acc[MAX_BODIES];

    /* Find a time-reversible mean acceleration over the interval dt. */
    ApproximateMovement(sim->nbodies, dt, sim->body, sim->state, next_state, curr_acc, middle_acc, next_acc);

    /* Apply the mean acceleration for half the time (dt/2) to find middle state (position and velocity). */
    MoveAllBodies(sim->nbodies, sim->state, middle_state, middle_acc, dt / 2.0);

    p = 2.0 / dt;
    dt2 = dt * dt;
    dt3 = dt * dt2;
    dt4 = dt2 * dt2;

    /* Find the unique best-fit parabolas for the 3 acceleration components (x, y, z). */
    for (b = 0; b < sim->nbodies; ++b)
    {
        for (k = 0; k < 3; ++k)     /* iterate through the components of the vectors: 0=x, 1=y, 2=z */
        {
            J = curr_acc[b].c[k];
            K = middle_acc[b].c[k];
            L = next_acc[b].c[k];

            /* Find coefficients of the best fit parabola */
            A = (L + J)/2 - K;
            B = (L - J)/2;

            E = A*p*p;
            F = (B - 2*A)*p;
            G = J;

            /* acceleration = Et^2 + Ft + G */
            /* Integrating, we get: velocity = (1/3)Et^3 + (1/2)Ft^2 + Gt + V0 */
            next_state[b].vel.c[k] = (E/3)*dt3 + (F/2)*dt2 + G*dt + sim->state[b].vel.c[k];

            /* Integrating again, we get: position = (1/12)Et^4 + (1/6)Ft^3 + (1/2)Gt^2 + V0*t + r0 */
            next_state[b].pos.c[k] = (E/12)*dt4 + (F/6)*dt3 + (G/2)*dt2 + sim->state[b].vel.c[k]*dt + sim->state[b].pos.c[k];
        }
    }

    CopyStates(sim->nbodies, next_state, sim->state);
    sim->tt += dt;
}
