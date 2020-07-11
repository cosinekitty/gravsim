/* sstest.js -- Test of gravsim using actual Solar System data. */

const fs = require('fs');
const gravsim = require('./gravsim.js');
const ss = JSON.parse(fs.readFileSync('ephemeris.json'));

/*
    The following are the reciprocal masses of the solar system bodies,
    from the DE-405 model. These data were copied from NOVAS C 3.1 file novascon.c.
    Each number is the mass of the Sun divided by the mass of the given body.
*/
const recipMasses = {
    Sun:              1.0,
    Mercury:    6023600.0,
    Venus:       408523.71,
    Earth:       328900.5614,      // actually the Earth/Moon system as a whole
    Mars:       3098708.0,
    Jupiter:       1047.3486,
    Saturn:        3497.898,
    Uranus:       22902.98,
    Neptune:      19412.24,
    Pluto:    135200000.0,
};

function main() {
    // Convert (Sun/body) mass ratios to (body/Earth) mass ratios.
    for (let name in recipMasses) {
        console.log(`${name}: ${recipMasses.Earth / recipMasses[name]}`)
    }
}

main();
