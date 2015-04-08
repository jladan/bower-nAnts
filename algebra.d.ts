/** nAnjs Matrix and vector algebra library
 *
 * This module contains all of the matrix functions and algorithms,
 * notably:
 *  - Addition, multiplication
 *  - Inversion
 *  - Decomposition (LU, QR, etc.).
 *  - Determinants
 *
 * Also included are matrix and vector creation algorithms
*/
declare module Algebra {
    class Size {
        n: number;
        m: number;
        constructor(n: number, m: number);
        copy(): Size;
    }
    /** Matrix object
     *
     * Creates a matrix of size `size`.
     * If `buffer` is an array, the matrix points to it, instead of creating a new array
     *
     * There are no guarantees of the values being initialized.
     * The matrix is in row-major format (because of row-pivoting)
     */
    class Matrix {
        array: number[];
        size: Size;
        constructor(size: [number, number] | Size, buffer?: number[]);
        /** return a copy of the matrix */
        copy(): Matrix;
        /**  A matrix of size `size` with ones on the diagonal */
        static identity(size: number | Size): Matrix;
        /** A matrix of size `size` filled with random entries */
        static random(size: Size): Matrix;
        /** Return the sum of two matrices
         */
        add(m: Matrix): Matrix;
        /** Return the difference of two matrices
         */
        subtract(m: any): Matrix;
        /** Multiply the matrix on the *right* by m */
        multiply(m: any): Matrix;
        /** transpose a matrix
         */
        transpose(m: any): Matrix;
        /** A solver for upper triangular matrices, also known as a back-solve.
         *
         * The Matrix is not checked to be upper triangular, so no error will be
         * thrown if a non-zero element is below the diagonal. This makes it
         * possible to use on an LU decomposition that's stored in a single matrix.
         */
        utSolve(b: any): Matrix;
        /** A solver for lower triangular matrices, also known as a forward-solve.
         *
         * The Matrix is not checked to be lower triangular, so no error will be
         * thrown if a non-zero element is below the diagonal. This makes it
         * possible to use on an LU decomposition that's stored in a single matrix.
         *
         * Because it's designed to be used for LU decomposion, it is assumed that
         * the diagonal elements are all one.
         */
        ltSolve(b: any): Matrix;
        /** LU decomposition of a matrix
         *
         * L is stored below the diagonal.
         */
        luDecompose(): Matrix;
        /** Extract L from the LU decomposition
         */
        grabL(): Matrix;
        /** Extract U from the LU decomposition
         */
        grabU(): Matrix;
        /** A solver for general matrices by LU decomposition
         *
         * When the forward and back solves are designed to accept matrices rather
         * than just vectors, this can be used to invert a matrix.
         */
        solve(b: any): Matrix;
    }
}
