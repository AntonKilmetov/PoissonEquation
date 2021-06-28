using System;
using System.IO;
using System.Text;

namespace PoissonEquation
{
    class Program
    {
        static void Main(string[] args)
        {
            var galerkinMethod = new GalerkinMethod(0.01);
            var solution = galerkinMethod.ComputePhi();
        }
    }

    class GalerkinMethod
    {
        readonly double[,] A; // матрица А
        readonly double[,] L; // матрица L
        readonly double[,] U; // матрица U
        readonly double[,] B; // матрица B
        readonly double[] X; // разбиение x
        readonly double Step; // параметр дискретизации h
        const double pi = Math.PI;
        public GalerkinMethod(double step)
        {
            Step = step;
            int n = (int)(1.0 / step) + 1;
            A = new double[n, n];
            L = new double[n, n];
            U = new double[n, n];
            B = new double[n, 1];
            X = new double[n];
            for (int i = 0; i < n; i++)
                X[i] = i * step;
            SetAMatrixElemets();
            SetBMatrixElemets();
            SetLUMatrixElemets();
        }
        /// <summary>
        /// Заполнение матрицы А
        /// </summary>
        private void SetAMatrixElemets()
        {
            for (int i = 0; i < A.GetLength(0); i++)
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    if (i == j)
                        A[i, j] = 2.0 / Step;
                    else if (i == j + 1 || i == j - 1)
                        A[i, j] = -1.0 / Step;
                    else
                        A[i, j] = 0;
                }
        }
        /// <summary>
        /// Заполнение матрицы B
        /// </summary>
        private void SetBMatrixElemets()
        {
            for (int i = 1; i < B.GetLength(0) - 1; i++)
            {
                if (i <= B.GetLength(0) / 2)
                    B[i, 0] = (-pi / (3 * Step)) * (Math.Pow(X[i - 1], 4) - 2 * Math.Pow(X[i], 4) + Math.Pow(X[i + 1], 4));
                //else
                //    B[i, 0] = ((-4 * pi) / Step) * ((Math.Pow(X[i - 1], 4) - 2 * Math.Pow(X[i], 4) + Math.Pow(X[i + 1], 4)) / 12
                //                               - (Math.Pow(X[i - 1], 3) - 2 * Math.Pow(X[i], 3) + Math.Pow(X[i + 1], 3)) / 3
                //                               + (Math.Pow(X[i - 1], 2) - 2 * Math.Pow(X[i], 2) + Math.Pow(X[i + 1], 2)) / 2);
                else 
                    B[i, 0] = B[B.GetLength(0) - 1 - i, 0]; //Можно заменить при симметриии относительно середины ро(х)
            }
        }
        /// <summary>
        /// Заполнение матриц L и U
        /// </summary>
        private void SetLUMatrixElemets()
        {
            var omega = new double[X.Length];
            var u = new double[X.Length];
            for (int i = 0; i < X.Length; i++)
                for (int j = 0; j < X.Length; j++)
                {
                    if (i == 0)
                        omega[i] = A[i, i];
                    else
                        omega[i] = A[i, i] - A[i, i - 1] * u[i - 1];
                    if (i != X.Length - 1)
                        u[i] = A[i, i + 1] / omega[i];
                }
            for (int i = 0; i < X.Length; i++)
                for (int j = 0; j < X.Length; j++)
                {
                    if (i == j)
                    {
                        U[i, j] = 1;
                        L[i, j] = omega[i];
                    }
                    else if (i == j - 1)
                    {
                        U[i, j] = u[i];
                        L[i, j] = 0;
                    }
                    else if (i == j + 1) 
                    {
                        U[i, j] = 0;
                        L[i, j] = A[i, j];
                    }
                    else
                    {
                        U[i, j] = 0;
                        L[i, j] = 0;
                    }
                }
        }
        /// <summary>
        /// Вычисление y; Ly=B
        /// </summary>
        /// <returns>Матрица y</returns>
        private double[] ComputeY()
        {
            var y = new double[X.Length];
            for (int i = 0; i < X.Length; i++)
            {
                if (i == 0)
                    y[i] = B[i, 0] / L[i, i];
                else
                    y[i] = (B[i, 0] - L[i, i - 1] * y[i - 1]) / L[i, i];
            }
            return y;
        }
        /// <summary>
        /// Вычисление a; U*a=y
        /// </summary>
        /// <returns>матрица а = phi</returns>
        public double[] ComputePhi()
        {
            var y = ComputeY();
            var phi = new double[X.Length];
            for (int i = X.Length - 1; i >= 0; i--)
            {
                if (i == X.Length - 1)
                    phi[i] = y[i];
                else
                    phi[i] = y[i] - U[i, i + 1] * phi[i + 1];
            }
            return phi;
        }
    }
}
