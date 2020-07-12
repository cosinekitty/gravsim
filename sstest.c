/*
    sstest.c  -  by Don Cross

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

#include <stdio.h>
#include <stdlib.h>
#include "gravsim.h"

static const size_t nbodies = 10;


static void SetBody(
    body_t *body, const char *name, double gm,
    double rx, double ry, double rz,
    double vx, double vy, double vz)
{
    body->name = name;
    body->gm = gm;
    body->state[0].pos = Vector(rx, ry, rz);
    body->state[0].vel = Vector(vx, vy, vz);
}


static int InitSolarSystem(sim_t *sim)
{
    int error;

    sim->tt = 0.0;
    sim->nbodies = (int)nbodies;
    sim->body = calloc(nbodies, sizeof(body_t));
    if (sim->body == NULL)
        FAIL("Out of memory!\n");

    SetBody(
        &sim->body[0], "Sun", 0.2959122082855911e-03,
        -7.1364589399065259e-03, -2.6470228609322332e-03, -9.2294970156656141e-04,
        +5.3784602410181226e-06, -6.7581870218649809e-06, -3.0328502580586604e-06);

    SetBody(
        &sim->body[1], "Mercury", 0.4912547451450812e-10,
        -1.3723006195467538e-01, -4.0324074408148058e-01, -2.0141225506190355e-01,
        +2.1371774112416420e-02, -4.9330574149022378e-03, -4.8504663531593545e-03);

    SetBody(
        &sim->body[2], "Venus", 0.7243452486162703e-09,
        -7.2543875484147236e-01, -4.8921273467320933e-02, +2.3717693023504526e-02,
        +8.0349602705784566e-04, -1.8498595719303294e-02, -8.3727680737444125e-03);

    SetBody(
        &sim->body[3], "Earth", 0.8997011346712499e-09,
        -1.8429524682327703e-01, +8.8475983851898110e-01, +3.8381376140494267e-01,
        -1.7197730582930743e-02, -2.9096002963053319e-03, -1.2615424279804276e-03);

    SetBody(
        &sim->body[4], "Mars", 0.9549535105779258e-10,
        +1.3835794628924982e+00, -1.2458004988146892e-03, -3.7883117515271375e-02,
        +6.7687793460626899e-04, +1.3807279375402957e-02, +6.3148674835543615e-03);

    SetBody(
        &sim->body[5], "Jupiter", 0.2825345909524226e-06,
        +3.9940404222298844e+00, +2.7339319061545413e+00, +1.0745894287353270e+00,
        -4.5629355212736143e-03, +5.8747037012365335e-03, +2.6292702270069392e-03);

    SetBody(
        &sim->body[6], "Saturn", 0.8459715185680659e-07,
        +6.3992748800141177e+00, +6.1720103478444583e+00, +2.2738496033938227e+00,
        -4.2869717425808437e-03, +3.5215864712979240e-03, +1.6388988371031218e-03);

    SetBody(
        &sim->body[7], "Uranus", 0.1292024916781969e-07,
        +1.4424723139268364e+01, -1.2508906775795596e+01, -5.6826051942721962e+00,
        +2.6834832774578900e-03, +2.4552472167487850e-03, +1.0373771677589703e-03);

    SetBody(
        &sim->body[8], "Neptune", 0.1524358900784276e-07,
        +1.6804919524159171e+01, -2.2982756707473023e+01, -9.8253477507922486e+00,
        +2.5846540556240267e-03, +1.6616650376509003e-03, +6.1578224469068194e-04);

    SetBody(
        &sim->body[9], "Pluto", 0.2188699765425970e-11,
        -9.8824799249935378e+00, -2.7981499149074953e+01, -5.7546082780601502e+00,
        +3.0341297634731501e-03, -1.1343428301178919e-03, -1.2681607296589918e-03);

    Accelerations(sim, 0);

    error = 0;
fail:
    return error;
}


int InitFinalState(sim_t *sim)
{
    int error;

    sim->tt = 36000.0;
    sim->nbodies = (int)nbodies;
    sim->body = calloc(nbodies, sizeof(body_t));
    if (sim->body == NULL)
        FAIL("Out of memory!\n");

    SetBody(
        &sim->body[0], "Sun", 0.2959122082855911e-03,
        +7.7442330999319582e-03, -2.8958174622971387e-03, -1.4843523935615082e-03,
        +3.7976242804768201e-06, +6.8873739539434805e-06, +2.8328030391439036e-06);

    SetBody(
        &sim->body[1], "Mercury", 0.4912547451450812e-10,
        +2.9998909445899702e-01, -2.5167075958321738e-01, -1.6463825444706792e-01,
        +1.4347702925469906e-02, +1.9275892909860873e-02, +8.8151240781442156e-03);

    SetBody(
        &sim->body[2], "Venus", 0.7243452486162703e-09,
        -1.2730466485862729e-01, -6.5678416128711048e-01, -2.8733612731354014e-01,
        +1.9741393720748273e-02, -3.0433142723558051e-03, -2.6179095787213476e-03);

    SetBody(
        &sim->body[3], "Earth", 0.8997011346712499e-09,
        +5.4268057037700779e-01, -7.9519201392980288e-01, -3.4479026527466322e-01,
        +1.4350563192874279e-02, +8.2605607839611600e-03, +3.5787087909035209e-03);

    SetBody(
        &sim->body[4], "Mars", 0.9549535105779258e-10,
        -1.3233061280808283e+00, +8.8734071813401461e-01, +4.4254312899760900e-01,
        -7.8463585498583580e-03, -9.1754296426277467e-03, -3.9988125767134418e-03);

    SetBody(
        &sim->body[5], "Jupiter", 0.2825345909524226e-06,
        -4.6210326953510954e+00, +2.4621350057089506e+00, +1.1674708023170912e+00,
        -3.9185718437625807e-03, -5.6887319737186108e-03, -2.3428130677038798e-03);

    SetBody(
        &sim->body[6], "Saturn", 0.8459715185680659e-07,
        -9.4886338896573026e+00, -3.3627229859393043e-01, +2.7058431810050654e-01,
        -1.8651175280703436e-04, -5.1701674264839478e-03, -2.1281654055284450e-03);

    SetBody(
        &sim->body[7], "Uranus", 0.1292024916781969e-07,
        +1.9467801910417869e+01, +4.3750090028448021e+00, +1.6411273005424660e+00,
        -9.4602209093205781e-04, +3.3311613626557600e-03, +1.4722799034082145e-03);

    SetBody(
        &sim->body[8], "Neptune", 0.1524358900784276e-07,
        -2.8549810690325550e+01, +8.8069648185486447e+00, +4.3155638698954348e+00,
        -1.0362322640330788e-03, -2.7392211213405310e-03, -1.0953738661569298e-03);

    SetBody(
        &sim->body[9], "Pluto", 0.2188699765425970e-11,
        +4.0183705547014910e+01, +2.7566070190543023e+01, -3.5042261894867393e+00,
        -9.2786393703440592e-04, +1.7551921708389899e-03, +8.2733972886151190e-04);

    error = 0;
fail:
    return error;
}


void Compare(sim_t *sim, sim_t *goal)
{

}


int main()
{
    int error  = 1;
    sim_t sim  = {0.0, 0, NULL};
    sim_t goal = {0.0, 0, NULL};
    const int nsteps = 1000;
    double dt;
    int n;

    CHECK(InitSolarSystem(&sim));
    CHECK(InitFinalState(&goal));
    dt = (goal.tt - sim.tt) / nsteps;
    for (n=0; n < nsteps; ++n)
        SimUpdate1(&sim, dt);

    Compare(&sim, &goal);

fail:
    if (sim.body != NULL) free(sim.body);
    if (goal.body != NULL) free(goal.body);
    return error;
}
