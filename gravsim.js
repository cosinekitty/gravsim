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
    const AU = 1.49597870691e+11;           // astronomical unit [m/au]
    const SecondsPerDay = 24 * 3600;

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

    function Average(a, b) {
        return [(a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2];
    }

    function EvalPoly(poly, t) {
        let tpower = 1;
        let sum = 0;
        for (let coeff of poly) {
            sum += coeff * tpower;
            tpower *= t;
        }
        return sum;
    }

    gravsim.VectorError = function(a, b) {
        const x = Subtract(a, b);
        return Math.sqrt(Dot(x, x));
    }

    class Simulator {
        constructor(gm, initialStates) {
            this.gm = gm;
            this.state = initialStates;
        }

        Accelerations(state) {
            // Calculate the net acceleration vectors experienced by
            // each body due to the gravitational pull of all other bodies.
            const acc = {};
            for (let [name, body] of Object.entries(state)) {
                acc[name] = [0, 0, 0];
                for (let [otherName, otherBody] of Object.entries(state)) {
                    if (body !== otherBody) {
                        let dr = Subtract(otherBody.pos, body.pos);
                        let r2 = Dot(dr, dr);
                        let r = Math.sqrt(r2);
                        let g = this.gm[otherName] / r2;    // acceleration in [AU/day^2]
                        let a = Multiply(g/r, dr);          // acceleration vector toward other body
                        acc[name] = Add(acc[name], a);
                    }
                }
            }
            return acc;
        }

        AverageAccelerations(a, b) {
            let c = {}
            for (let name in a) {
                c[name] = Average(a[name], b[name]);
            }
            return c;
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

        Update1(dt) {
            // A naive update uses the accelerations of each body
            // at the start of the increment as if they are constant
            // over the time interval in that increment.
            // In reality, the accelerations change slightly as
            // the bodies move along their respective orbits.
            let acc = this.Accelerations(this.state);
            this.state = this.Movement(this.state, acc, dt);
            return this.state;
        }

        Update2(dt) {
            // This is a much more accurate update algorithm than Update1.
            // It starts the same way as Update1, but refines the answer
            // by using the average of the acceleration vectors for the
            // beginning and ending of the increment.
            // This causes the estimated changes in velocity and position
            // to change, which in turn causes the acceleration estimates to change.
            // The process is iterated until the position vector converges.
            let acc1 = this.Accelerations(this.state);
            let guess = this.Movement(this.state, acc1, dt);
            let diff;
            do {
                // Calculate the accelerations at the guessed position.
                let acc2 = this.Accelerations(guess);

                // Assume the overall acceleration effect is approximated by the
                // average of the beginning acceleration and ending acceleration.
                let acc = this.AverageAccelerations(acc1, acc2);

                // Update the guess using the refined acceleration.
                let refined = this.Movement(this.state, acc, dt);

                // Mercury experiences the most acceleration of all the bodies.
                // Use it as a gauge of convergence.
                diff = gravsim.VectorError(guess.Mercury.pos, refined.Mercury.pos);

                guess = refined;
            } while (diff > 1.0e-15);
            this.state = guess;
            return this.state;
        }

        Update3(dt) {
            // A refinement of the idea in Update2.
            // Find a best-fit quadratic function for the accelerations
            // at times 0, dt/2, and dt.
            // Integrate that acceleration to find a cubic velocity function.
            // Integrate the velocity to find a quartic position function.
            // Calculate the resulting positions at times 0, dt/2, and dt;
            // this yields three refined estimates for accelerations.
            // Iterate until convergence!
            let acc1 = this.Accelerations(this.state);
            let state2 = this.Movement(this.state, acc1, dt/2);
            let acc2 = this.Accelerations(state2);
            let state3 = this.Movement(state2, acc2, dt/2);
            let acc3 = this.Accelerations(state3);

            const p = 2 / dt;
            for (let n=0; n<2; ++n) {        // assume enough loops for convergence
                for (let [name, body] of Object.entries(this.state)) {
                    let pos2 = [];
                    let pos3 = [];
                    let vel3 = [];
                    // Solve each component of acceleration, velocity, and position independently.
                    for (let i=0; i<3; ++i) {
                        // Find the best-fit parabola through the three consecutive acceleration component values.
                        let A = (acc3[name][i] + acc1[name][i])/2 - acc2[name][i];
                        let B = (acc3[name][i] - acc1[name][i])/2;
                        let E = A*p*p;
                        let F = (B - 2*A)*p;
                        let G = acc1[name][i];

                        // Now we have a parabola that approximates acceleration:
                        // a(t) = Et^2 + Ft + G
                        // Integrate to get velocity:
                        // v(t) = (1/3)Et^3 + (1/2)Ft^2 + Gt + v(0)
                        let vel_poly = [body.vel[i], G, F/2, E/3];
                        vel3[i] = EvalPoly(vel_poly, dt);

                        // Integrate again to get position:
                        // r(t) = (1/12)Et^4 + (1/6)Ft^3 + (1/2)Gt^2 + v(0)t + r(0)
                        let pos_poly = [body.pos[i], body.vel[i], G/2, F/6, E/12];
                        pos2[i] = EvalPoly(pos_poly, dt/2);
                        pos3[i] = EvalPoly(pos_poly, dt);
                    }
                    // Now we have updated position and velocity vectors for this body
                    // at times dt/2 and dt. Replace them inside the estimated state objects.
                    // Tricky: we don't bother to update state2 velocity, because it has no effect
                    // on the final results. We *do* need state3 velocity because it will end up
                    // being part of the complete system state when we are done converging.
                    state2[name].pos = pos2;
                    state3[name].pos = pos3;
                    state3[name].vel = vel3;
                }
                acc2 = this.Accelerations(state2);
                acc3 = this.Accelerations(state3);
            }
            this.state = state3;
            return this.state;
        }
    }

    gravsim.MakeSimulator = function(gm, initialStates) {
        return new Simulator(gm, initialStates);
    }
})(typeof exports==='undefined' ? (this.gravsim={}) : exports);
