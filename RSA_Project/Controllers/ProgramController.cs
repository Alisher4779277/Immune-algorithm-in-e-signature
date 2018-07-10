using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Web.Mvc;

namespace RSA_Project.Controllers
{

    public class Programm
    {

        Matrix mat = new Matrix();

        double[,] A = new double[8, 8];
        double[,] B = new double[8, 8];

        //A = Matrix.matrix_A();
        //B = mat.matrix_B();


        //var matrix = DenseMatrix.OfArray(A);

        //var svd = matrix.Svd();

        //var reconstruct = svd.U * svd.W * svd.VT;
        //Console.WriteLine("A=");
        //Console.WriteLine(reconstruct.ToString());
        //Console.WriteLine();


        //Console.WriteLine("U=");
        double[,] U = new double[8, 8];
        //for (int i = 0; i < 8; i++)
        //{
        //    for (int j = 0; j < 8; j++)
        //    {
        //        U[i, j] = svd.U[i, j];
        //        Console.Write(Math.Round(U[i, j], 4) + "\t");
        //    }
        //    Console.WriteLine();
        //}

        //Console.WriteLine("S=");
        //Console.WriteLine(svd.W.ToString());
        //Console.WriteLine();

        //Console.WriteLine("V=");
        //double[,] V = new double[8, 8];
        //for (int i = 0; i < 8; i++)
        //{
        //    for (int j = 0; j < 8; j++)
        //    {
        //        V[i, j] = svd.VT[j, i];
        //        Console.Write(Math.Round(V[i, j], 4) + "\t");
        //    }
        //    Console.WriteLine();
        //}

        //Hesh h = new Hesh();
        //h.U(U);
        //h.V(V);
        //h.Hesh_Calc(U, V, B);

        //Console.ReadKey();


    }

}