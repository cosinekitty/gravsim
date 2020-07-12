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


function Test(ss, method) {
    const sim = gravsim.MakeSimulator(masses, ss.data[0].body);
    const nsteps = 1000;
    const dt = ss.dt / nsteps;
    const AU = 1.4959787069098932e+11;      // astronomical unit [m/au]

    console.log(`BEGINNING: ${method}`);

    const p1 = sim.Momentum();
    console.log(`Initial momentum: ${p1}`);

    for (let n=0; n+1 < ss.data.length; ++n) {
        for (let i=0; i < nsteps; ++i) {
            sim[method](dt);
        }
    }

    const se = sim.state.Earth;
    const ce = ss.data[ss.data.length-1].body.Earth;
    const err_au = gravsim.VectorError(se.pos, ce.pos);
    const err_km = err_au * (AU / 1000);
    const p2 = sim.Momentum();

    console.log(`Simulated Earth pos = ${se.pos}`);
    console.log(`Correct   Earth pos = ${ce.pos}`);
    console.log(`Error = ${err_au} AU (${err_km} km)`);
    console.log(`Final momentum: ${p2}`);
    console.log(`Momentum error: ${gravsim.VectorError(p2, p1)}`);
    console.log();
}


function main() {
    const ss = JSON.parse(fs.readFileSync('ephemeris.json'));
    Test(ss, 'Update1');
    Test(ss, 'Update2');
}


main();
