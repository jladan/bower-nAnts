/** nAnts stochastic simulation module
 */
export declare enum Method {
    Euler = 0,
    Milstein = 1,
    Colour = 2,
}
export declare class Solution {
    t: Float32Array | Float64Array;
    result: Float32Array | Float64Array;
    N: number;
    constructor(t: any, result: any, n: any);
    getDimension(d: number): Float32Array | Float64Array;
    /** Extracts a 'trail' from the solution: pairs of (x,t)
     * d is the index of the dimension (starting from 0)
     * outStep is the number of time-steps between samples. Default to 1
     */
    getTrail(d: number, outStep?: number): [number, number][];
    /** Extracts a trail in phase-space, coordinate (x_1, x_2)
     * d1, and d2  are the indices of the coordinates (starting from 0)
     * outStep is the number of time-steps between samples. Default to 1
     */
    getPhase(d1: number, d2: number, outStep?: number): [number, number][];
}
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
export declare function euler(A: (x: number[], t: number, p: number[]) => number[], D: (x: number[], t: number, p: number[]) => number[], initial: number[], dt: number, t_final: number, pA: number[], pD: number[]): Solution;
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
export declare function milstein(A: (x: number[], t: number, p: number[]) => number[], D: (x: number[], t: number, p: number[]) => number[], Dy: (x: number[], t: number, p: number[]) => number[], initial: number[], dt: number, t_final: number, pA: number[], pD: number[]): Solution;
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
export declare function colour(A: (x: number[], t: number, p: number[]) => number[], D: (x: number[], t: number, p: number[]) => number[], sigma: number, tau: number, initial: number[], dt: number, t_final: number, pA: number[], pD: number[]): Solution;
/** The Box Muller transform
 */
export declare function boxMuller(): number;
/** Box Muller function that returns both Gaussians
 */
export declare function boxMuller2(): number[];
export declare var gaussian: typeof boxMuller;
