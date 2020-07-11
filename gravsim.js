'use strict';
/**
    @preserve

    Solar System gravity simulator for JavaScript (browser and Node.js).
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
(function(gravsim){
    'use strict';

    const C  = 299792458.0;                 // speed of light [m/s]
    const AU = 1.4959787069098932e+11;      // astronomical unit [m/au]
    const SecondsPerDay = 24 * 3600;

    // gravitation constant [AU^3 * kg^(-1) * day^(-2)]
    const G  = 6.67430e-11 * (SecondsPerDay * SecondsPerDay) / (AU*AU*AU);

    function Dot(a, b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    function Subtract(a, b) {
        return [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
    }

    function Add(a, b) {
        return [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
    }

    function Multiply(k, v) {
        return [k*v[0], k*v[1], k*v[2]];
    }

    gravsim.VectorError = function(a, b) {
        const x = Subtract(a, b);
        return Math.sqrt(Dot(x, x));
    }

    class Simulator {
        constructor(masses, initialStates) {
            this.mass = masses;
            this.state = initialStates;
        }

        Accelerations(state) {
            // Calculate the net acceleration vectors experienced by
            // each body due to the gravitational pull of all other bodies.
            const acc = {};
            for (let [name, body] of Object.entries(state)) {
                acc[name] = [0.0, 0.0, 0.0];
                for (let [otherName, otherBody] of Object.entries(state)) {
                    if (body !== otherBody) {
                        let dr = Subtract(otherBody.pos, body.pos);
                        let r2 = Dot(dr, dr);
                        let r = Math.sqrt(r2);
                        let g = G * this.mass[otherName] / r2;      // acceleration in [AU/day^2]
                        let a = Multiply(g/r, dr);      // acceleration vector toward other body
                        acc[name] = Add(acc[name], a);
                    }
                }
            }
            return acc;
        }

        Movement(state, acc, dt) {
            // Apply movement to every body in 'state' using provided acceleration vectors.
            // Return a new state object containing bodies with updated
            // position and velocity vectors.
            let newState = {};
            for (let [name, body] of Object.entries(state)) {
                // pos' = pos + vel*dt + (1/2)acc*dt^2
                // vel' = vel + acc*dt
                let dv = Multiply(dt, acc[name]);
                let dr = Add(Multiply(dt/2, dv), Multiply(dt, body.vel));
                newState[name] = {
                    vel: Add(body.vel, dv),
                    pos: Add(body.pos, dr)
                };
            }
            return newState;
        }

        NaiveUpdate(dt) {
            let acc = this.Accelerations(this.state);
            this.state = this.Movement(this.state, acc, dt);
            return this.state;
        }
    }

    gravsim.MakeSimulator = function(masses, initialStates) {
        return new Simulator(masses, initialStates);
    }
})(typeof exports==='undefined' ? (this.gravsim={}) : exports);
