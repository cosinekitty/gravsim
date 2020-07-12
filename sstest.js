/* sstest.js -- Test of gravsim using actual Solar System data. */

const fs = require('fs');
const gravsim = require('./gravsim.js');

/*
    https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf

    Page 10 in the above document describes the constants used in the DE405 ephemeris.
    The following are G*M values (gravity constant * mass) in [au^3 / day^2].
    This side-steps issues of not knowing the exact values of G and masses M[i];
    the products GM[i] are known extremely accurately.
*/
const gm = {
    Sun:        0.2959122082855911e-03,
    Mercury:    0.4912547451450812e-10,
    Venus:      0.7243452486162703e-09,
    Earth:      0.8997011346712499e-09,     // Earth/Moon barycenter (EMB)
    Mars:       0.9549535105779258e-10,
    Jupiter:    0.2825345909524226e-06,
    Saturn:     0.8459715185680659e-07,
    Uranus:     0.1292024916781969e-07,
    Neptune:    0.1524358900784276e-07,
    Pluto:      0.2188699765425970e-11
};


function Test(ss, method) {
    const sim = gravsim.MakeSimulator(gm, ss.data[0].body);
    const nsteps = 1000;
    const dt = ss.dt / nsteps;
    const AU = 1.49597870691e+11;           // astronomical unit [m/au]

    console.log();
    console.log(method);

    for (let n=0; n+1 < ss.data.length; ++n) {
        for (let i=0; i < nsteps; ++i) {
            sim[method](dt);
        }
    }

    for (let name in sim.state) {
        const se = sim.state[name];
        const ce = ss.data[ss.data.length-1].body[name];
        const err_au = gravsim.VectorError(se.pos, ce.pos);
        const err_km = err_au * (AU / 1000);
        console.log(`${name.padEnd(10)} error = ${err_au.toFixed(6).padStart(10)} AU (${err_km.toFixed(1).padStart(12)} km)`);
    }
}


function main() {
    const ss = JSON.parse(fs.readFileSync('ephemeris.json'));
    Test(ss, 'Update1');
    Test(ss, 'Update2');
    Test(ss, 'Update3');
}


main();
