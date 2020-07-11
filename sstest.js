/* sstest.js -- Test of gravsim using actual Solar System data. */

const fs = require('fs');
const gravsim = require('./gravsim.js');
const ss = JSON.parse(fs.readFileSync('ephemeris.json'));

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

function main() {
    const sim = gravsim.MakeSimulator(masses, ss.data[0].body);
    const dt = 1.0;     // time increment in days
    sim.Update(1.0);
}

main();
