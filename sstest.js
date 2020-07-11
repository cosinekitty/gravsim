/* sstest.js -- Test of gravsim using actual Solar System data. */

const fs = require('fs');
const gravsim = require('./gravsim.js');
const ss = JSON.parse(fs.readFileSync('ephemeris.json'));

// Masses of Solar System bodies relative to the Earth/Moon system
const masses = {
    Sun:        328900.5614,
    Mercury:         0.05460199239657348,
    Venus:           0.8050954041321127,
    Earth:           1.0,                   // actually the Earth/Moon system as a whole
    Mars:            0.10614119220010404,
    Jupiter:       314.03160456795376,
    Saturn:         94.02805953747078,
    Uranus:         14.360601170677354,
    Neptune:        16.94294740843921,
    Pluto:           0.0024326964600591716
};

const EarthMassKg = 5.9722e+24;
const MoonMassKg  = 7.3420e+22;


function main() {
    // Convert masses to kilograms
    for (let name in masses) {
        console.log(`${name}: ${masses[name] * (EarthMassKg + MoonMassKg)},`);
    }
}

main();
