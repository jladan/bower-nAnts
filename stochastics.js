/** nAnts stochastic simulation module
 */
var stochastics;
(function (stochastics) {
    (function (Method) {
        Method[Method["Euler"] = 0] = "Euler";
        Method[Method["Milstein"] = 1] = "Milstein";
        Method[Method["Colour"] = 2] = "Colour";
    })(stochastics.Method || (stochastics.Method = {}));
    var Method = stochastics.Method;
    ;
    var Solution = (function () {
        function Solution(t, result, n) {
            this.t = t;
            this.N = n;
            this.result = result;
        }
        // Helper functions
        Solution.prototype.getDimension = function (d) {
            return this.result.subarray(d * this.N, (d + 1) * this.N);
        };
        Solution.prototype.getTrail = function (d) {
            var x = this.getDimension(d);
            var i;
            var result = new Array();
            for (i = 0; i < this.N; i++)
                result.push([this.t[i], x[i]]);
            return result;
        };
        Solution.prototype.getPhase = function (d1, d2) {
            var x = this.getDimension(d1);
            var y = this.getDimension(d2);
            var i;
            var result = new Array();
            for (i = 0; i < this.N; i++)
                result.push([x[i], y[i]]);
            return result;
        };
        return Solution;
    })();
    stochastics.Solution = Solution;
    /** Euler-Maruyama method
     *
     * Numerically simulates the Langevin equation,
     *
     *      dy/dt = A(y,t) + D(y,t) * noise(t)
     *
     * For Gaussian white noise. The equation is actually interpreted in
     * the Ito sense.
     *
     * Arguments:
     * A - function of (y,t,[c]) provides the non-stochastic rates
     * D - function of (y,t,[c]) provides the stochastic rates
     * parameters - array of parameters to be sent to the functions A,D
     * the rest should be self-explanatory
     */
    function euler(A, D, initial, dt, t_final, pA, pD) {
        var d = initial.length;
        var N = Math.floor(t_final / dt);
        var result = new Float32Array(N * d);
        var t = new Float32Array(N);
        var n; // the most random number
        var sdt = Math.sqrt(dt);
        // initial conditions
        var i, j;
        var y_cur = initial.slice();
        for (j = 0; j < d; j++)
            result[j * N] = initial[j];
        for (i = 0; i < N - 1; i++) {
            var A_cur = A(y_cur, t[i], pA);
            var D_cur = D(y_cur, t[i], pD);
            t[i + 1] = (i + 1) * dt;
            for (j = 0; j < d; j++) {
                n = stochastics.gaussian();
                result[i + 1 + j * N] = y_cur[j] = y_cur[j] + dt * A_cur[j] + n * D_cur[j] * sdt;
            }
        }
        return new Solution(t, result, N);
    }
    stochastics.euler = euler;
    /** Milstein method
     *
     * Numerically simulates the Langevin equation,
     *
     *      dy/dt = A(y,t) + D(y,t) * noise(t)
     *
     * For Gaussian white noise. The equation is actually interpreted in
     * the Ito sense.
     *
     * Arguments:
     * A - function of (y,t,[c]) provides the non-stochastic rates
     * D - function of (y,t,[c]) provides the stochastic rates
     * Dy - derivative of D in y
     * parameters - array of parameters to be sent to the functions A,D
     * the rest should be self-explanatory
     */
    function milstein(A, D, Dy, initial, dt, t_final, pA, pD) {
        var d = initial.length;
        var N = Math.floor(t_final / dt);
        var result = new Float32Array(N * d);
        var t = new Float32Array(N);
        var sdt = Math.sqrt(dt);
        // initial conditions
        var i, j, n;
        var y_cur = initial.slice();
        for (j = 0; j < d; j++)
            result[j * N] = initial[j];
        for (i = 0; i < N - 1; i++) {
            var A_cur = A(y_cur, t[i], pA);
            var D_cur = D(y_cur, t[i], pD);
            var Dy_cur = Dy(y_cur, t[i], pD);
            t[i + 1] = (i + 1) * dt;
            for (j = 0; j < d; j++) {
                n = boxMuller2();
                result[i + 1 + j * N] = y_cur[j] = y_cur[j] + dt * A_cur[j] + n[0] * D_cur[j] * sdt - dt / 2 * D_cur[j] * Dy_cur[j] * (1 - n[1] * n[1]);
            }
        }
        return new Solution(t, result, N);
    }
    stochastics.milstein = milstein;
    /** Coloured Noise method
     *
     * Numerically simulates the Langevin equation,
     *
     *      dy/dt = A(y,t) + D(y,t) * noise(t)
     *
     * where the noise term is an Ornstein-Uhlenbeck process,
     * with autocorrelation:
     *
     *      <F(t)F(t-tau)> = sigma^2 exp(-|tau|/tau_c) .
     *
     * tau_c being the correlation time.
     *
     * Because the SDE is differentiable with coloured noise, the
     * Heun method is used to advance each time-step.
     *
     * Arguments:
     * A - function of (y,t,[c]) provides the non-stochastic rates
     * D - function of (y,t,[c]) provides the stochastic rates
     * sigma - standard deviation of gaussian noise
     * tau - correlation time of the noise
     * parameters - array of parameters to be sent to the functions A,D
     * the rest should be self-explanatory
     */
    function colour(A, D, sigma, tau, initial, dt, t_final, pA, pD) {
        var d = initial.length;
        var N = Math.floor(t_final / dt);
        var result = new Float32Array(N * d);
        var t = new Float32Array(N);
        var n = 4; // the most random number
        // coefficients to advance noise
        var rho = Math.exp(-dt / tau);
        var rhoc = sigma * Math.sqrt(1 - rho * rho); // XXX why sigma and not previous
        // Required variables
        var i, j;
        var k1 = new Array(d), k2 = new Array(d);
        var noiseCur = new Array(d), noiseNext = new Array(d);
        // Initial Conditions
        var y_cur = initial.slice();
        for (j = 0; j < d; j++) {
            result[j * N] = initial[j];
            noiseNext[j] = 0;
        }
        var y_next = y_cur.slice();
        for (i = 0; i < N - 1; i++) {
            // advance the noise and time
            t[i + 1] = (i + 1) * dt;
            for (j = 0; j < d; j++) {
                noiseCur[j] = noiseNext[j];
                n = boxMuller();
                noiseNext[j] = noiseCur[j] * rho + n * rhoc;
            }
            // use heun's method to advance the system
            var A_cur = A(y_cur, t[i], pA);
            var D_cur = D(y_cur, t[i], pD);
            for (j = 0; j < d; j++) {
                k1[j] = A_cur[j] + D_cur[j] * noiseCur[j];
                y_next[j] = y_cur[j] + dt * k1[j];
            }
            var A_next = A(y_next, t[i], pA);
            var D_next = D(y_next, t[i], pD);
            for (j = 0; j < d; j++) {
                k2[j] = A_next[j] + D_next[j] * noiseNext[j];
                result[i + 1 + j * N] = y_cur[j] = y_cur[j] + dt / 2 * (k1[j] + k2[j]);
            }
        }
        return new Solution(t, result, N);
    }
    stochastics.colour = colour;
    /** The Box Muller transform
     */
    function boxMuller() {
        var v1, v2, s, x;
        do {
            var u1 = Math.random();
            var u2 = Math.random();
            v1 = 2 * u1 - 1;
            v2 = 2 * u2 - 1;
            s = v1 * v1 + v2 * v2;
        } while (s > 1);
        x = Math.sqrt(-2 * Math.log(s) / s) * v1;
        return x;
    }
    stochastics.boxMuller = boxMuller;
    /** Box Muller function that returns both Gaussians
     */
    function boxMuller2() {
        var v1, v2, s, x, y;
        do {
            var u1 = Math.random();
            var u2 = Math.random();
            v1 = 2 * u1 - 1;
            v2 = 2 * u2 - 1;
            s = v1 * v1 + v2 * v2;
        } while (s > 1);
        x = Math.sqrt(-2 * Math.log(s) / s) * v1;
        y = Math.sqrt(-2 * Math.log(s) / s) * v2;
        return [x, y];
    }
    stochastics.boxMuller2 = boxMuller2;
    stochastics.gaussian = boxMuller;
})(stochastics || (stochastics = {}));
;
