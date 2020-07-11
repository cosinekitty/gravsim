/* sstest.js -- Test of gravsim using actual Solar System data. */

const fs = require('fs');
const gravsim = require('./gravsim.js');

/*
    The following table of masses was calculated from
    reciprocal masses of solar system bodies from DE-405,
    copied from NOVAS C 3.1 file novascon.c.
    These reciprocals were all of the form (mass of Sun) / (mass of body).
    Using Earth mass = 5.9722e+24 kg, Moon mass = 7.3420e+22 kg,
    we get Earth + Moon = 6.04562e+24 kg. All other bodies were scaled
    accordingly.
*/

const masses = {
    Sun:        1.988407812011068e+30,
    Mercury:    3.301028972725725e+23,
    Venus:      4.867300877129182e+24,
    Earth:      6.04562e+24,                // actually the Earth and Moon combined
    Mars:       6.41689314388793e+23,
    Jupiter:    1.8985157492081126e+27,
    Saturn:     5.684579173009241e+26,
    Uranus:     8.681873764947043e+25,
    Neptune:    1.0243062171140825e+26,
    Pluto:      1.470715837286293e+22
};


function Test(label, update) {
    const ss = JSON.parse(fs.readFileSync('ephemeris.json'));
    const sim = gravsim.MakeSimulator(masses, ss.data[0].body);
    const nsteps = 10000;
    const dt = ss.dt / nsteps;

    console.log(`${label}: Simulating ${nsteps} steps of ${dt} days per step.`);
    for (let i=0; i < nsteps; ++i) {
        update(sim, dt);
    }

    const se = sim.state.Earth;
    const ce = ss.data[1].body.Earth;
    console.log(`Simulated Earth pos = ${se.pos}`);
    console.log(`Correct   Earth pos = ${ce.pos}`);
    console.log(`Error = ${gravsim.VectorError(se.pos, ce.pos)} AU`);
    console.log();
}

Test('Naive',    (sim, dt) => sim.NaiveUpdate(dt))
Test('Improved', (sim, dt) => sim.Update(dt))
