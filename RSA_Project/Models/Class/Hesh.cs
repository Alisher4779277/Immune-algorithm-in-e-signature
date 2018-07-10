using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace RSA_Project.Controllers
{
    public class Hesh
    {
        double a, abcd;
        int b = 0, g = 0;
        double[] abc;
        public double[,] U_T;
        public double[] Hesh_result;
        public double Hesh_result_min = 0;

        public void U(double[,] U)
        {
            U_T = new double[U.GetLength(1), U.GetLength(1)];
            double[,] U_ans = new double[U.GetLength(1), U.GetLength(1)];

            for (int i = 0; i < U.GetLength(1); i++)
            {
                for (int j = 0; j < U.GetLength(1); j++)
                {
                    U_T[i, j] = U[j, i];
                }
            }


            Console.WriteLine("U ans =");
            for (int i = 0; i < U.GetLength(1); i++)
            {
                for (int j = 0; j < U.GetLength(1); j++)
                {
                    int mm = 0;
                    while (mm < U.GetLength(1))
                    {
                        U_ans[i, j] += U[j, mm] * U_T[mm, i];
                        mm++;
                    }
                    Console.Write(Math.Round(U_ans[i, j], 1) + "\t");
                }
                Console.WriteLine();
            }
        }
        public void V(double[,] V)
        {
            Console.WriteLine("V ans =");
            double[,] V_T = new double[V.GetLength(1), V.GetLength(1)];


            for (int i = 0; i < V.GetLength(1); i++)
            {
                for (int j = 0; j < V.GetLength(1); j++)
                {
                    V_T[i, j] = V[j, i];
                }
            }

            double[,] V_ans = new double[V.GetLength(1), V.GetLength(1)];
            for (int i = 0; i < V.GetLength(1); i++)
            {
                for (int j = 0; j < V.GetLength(1); j++)
                {
                    int mm = 0;
                    while (mm < V.GetLength(1))
                    {
                        V_ans[i, j] += V[j, mm] * V_T[mm, i];
                        mm++;
                    }
                    Console.Write(Math.Round(V_ans[i, j], 1) + "\t");
                }
                Console.WriteLine();
            }
        }
        public double Hesh_Calc(double[,] U_T, double[,] V, double[,] B)
        {
            Hesh_result = new double[B.GetLength(1)];
            abc = new double[B.GetLength(1)];

            while (b < B.GetLength(1))
            {

                for (int j = 0; j < B.GetLength(1); j++)
                {
                    abc[j] = 0;
                    for (int i = 0; i < B.GetLength(1); i++)
                    {
                        abc[j] += U_T[i, b] * B[i, j]; //abc=U*B
                    }
                }

                abcd = 0;
                for (int i = 0; i < B.GetLength(1); i++)
                {
                    abcd += abc[i] * V[i, b];//abcd=abc*V               
                }
                Hesh_result[b] = abcd;
                b++;
            }

            Hesh_result_min = Hesh_result[0];
            for (int i = 1; i < B.GetLength(1); i++)
            {
                if (Hesh_result[i] > Hesh_result_min)
                {
                    Hesh_result_min = Hesh_result[i];               // hesh minimumni topamiz u hesh qiymat hisoblanadi
                }
            }
            return Hesh_result_min;
        }


    }
}