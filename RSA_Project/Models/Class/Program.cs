using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Globalization;

namespace RSA_Project.Controllers
{
    public class Program
    {
            public int Start(double[,] A0, double[,] B)
        {
            List<double[,]> ArrayToList = new List<double[,]>();
            List<double[]> U = new List<double[]>();
            List<double[]> V = new List<double[]>();

            ArrayToList.Add(A0);


            int width = A0.GetLength(1);     // matritsani eniga o`lchami
            int height = A0.GetLength(0);    // matritsani bo`yiga o`lchami


            double[,] temp_matrix;               // vaqtinchalik matritsa qo`llab turamiz

            int r = 0;
            if (width > height) { r = height; }     // r minimumga teng
            else { r = width; }

            double[,] A = new double[height, width];                // A matritsa ni e'lon qilib olamiz 


            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    A[i, j] = A0[i, j];                     //A ni A0 teng deb olamiz
                }
            }

            double[] U0;
            double[] V0;



            double[] V1;
            V1 = new double[width];
            for (int i = 0; i < width; i++)             // V1 ni 1 ga tenglashtirib olamiz
            {
                V1[i] = 1;
            }

            double[] U1;
            U1 = new double[height];
            for (int i = 0; i < height; i++)           // U1 ni 1 ga tenglashtirib olamiz
            {
                U1[i] = 1;
            }

            int aylanish = 1;



            while (aylanish <= r)
            {              //r minimumga teng bo`guncha ishlaydi

                temp_matrix = new double[height, width];

                U0 = new double[height];
                V0 = new double[width];

                double V1_norm = 0;
                double U1_norm = 0;         // o`zgaruvchilarni e'lon qilish

                double[] temp;
                double[] temp2;         // o`zgaruvchilarni e'lon qilish

                double eps = 0.01;        //Epsilon 0,01 ga teng
                double ans = 1;           // ans 1 db olamiz

                double S0 = 0;
                double S1 = 0;

                int s = 1;          // sikl necha marta aylanganini tekshirib beradi



                while (ans >= eps)  // ans epsilondan kichkina yoki teng bo`lgunicha ishlaydi
                {
                    S0 = 0;
                    S1 = 0;



                    for (int i = 0; i < height; i++)
                    {
                        U0[i] = U1[i];
                    }



                    for (int i = 0; i < width; i++)
                    {
                        V0[i] = V1[i];
                    }



                    V1 = new double[width];
                    for (int i = 0; i < width; i++)
                    {
                        for (int j = 0; j < height; j++)
                        {
                            V1[i] += U0[j] * A[j, i];        //V1 ni topish  V1=U0*A
                        }
                    }



                    V1_norm = 0;
                    for (int i = 0; i < width; i++)
                    {
                        V1_norm += Math.Pow(V1[i], 2);       // v1 normni topish V1=sqrt(v1^2)
                    }
                    V1_norm = Math.Sqrt(V1_norm);




                    for (int i = 0; i < width; i++)
                    {
                        V1[i] = V1[i] / V1_norm;                 // v1 ni topish v1=v1/norm(v1)
                    }




                    U1 = new double[height];
                    for (int i = 0; i < height; i++)
                    {
                        for (int j = 0; j < width; j++)
                        {
                            U1[i] += V1[j] * A[i, j];        //U1=A*V1'   U1 ni topish
                        }
                    }



                    U1_norm = 0;
                    for (int i = 0; i < height; i++)
                    {
                        U1_norm += Math.Pow(U1[i], 2);        // norma(U1) ni topish 
                    }
                    U1_norm = Math.Sqrt(U1_norm);



                    for (int i = 0; i < height; i++)
                    {
                        U1[i] = U1[i] / U1_norm;                    //U1 ni topish U1=U1/norm(U1)
                    }




                    temp = new double[width];                        //S0 ni topish S0=U0*A*V0'
                    for (int i = 0; i < width; i++)
                    {
                        for (int j = 0; j < height; j++)
                        {
                            temp[i] += U0[j] * A[j, i];         // oldin U0' ni A ga ko`paytiramiz
                        }
                        S0 += temp[i] * V0[i];                // U0'*A dan hosil bo`lgan sonni v0' ga ko`paytiramiz
                    }





                    temp2 = new double[width];                            //S1 ni topish S1=U1'*A*V1'
                    for (int i = 0; i < width; i++)
                    {
                        for (int j = 0; j < height; j++)
                        {
                            temp2[i] += U1[j] * A[j, i];        // oldin U1' ni A ga ko`paytiramiz
                        }
                        S1 += temp2[i] * V1[i];             // U1'*A dan hosil bo`lgan sonni v1' ga ko`paytiramiz
                    }
                    ans = Math.Abs(S1 - S0);                     //ans ni topish
                    s++;


                }

                for (int i = 0; i < height; i++)
                {
                    for (int j = 0; j < width; j++)
                    {
                        temp_matrix[i, j] = U1[i] * V1[j] * S1;             // A1,2,3 ni topish = U1*V1'*S1
                    }
                }


                for (int i = 0; i < height; i++)
                {
                    for (int j = 0; j < width; j++)
                    {
                        A[i, j] = A[i, j] - temp_matrix[i, j];               // A-A[] ayriladi
                    }
                }


                U.Add(U0);
                V.Add(V0);
                ArrayToList.Add(temp_matrix);

                aylanish++;
            }



            double[,] sum = new double[height, width];
            int ppp = 0;

            foreach (var item in ArrayToList)
            {
                if (ppp != 0)
                {
                    for (int i = 0; i < item.GetLength(0); i++)
                    {
                        for (int j = 0; j < item.GetLength(1); j++)
                        {
                            sum[i, j] += item[i, j];
                        }
                    }
                    ppp++;
                }
                ppp++;
            }


            

            ppp = 0;
            double[,] U_matrix = new double[height, height];

            foreach (var item in U)
            {
                while (ppp != height)
                {
                    for (int j = 0; j < height; j++)                // U larni massivga yigib oldik
                    {
                        U_matrix[j, ppp] += item[j];
                    }
                    ppp++;
                    break;
                }

            }

            ppp = 0;
            double[,] V_matrix = new double[width, width];

            foreach (var item in V)
            {
                while (ppp != width)
                {
                    for (int j = 0; j < width; j++)
                    {
                        V_matrix[j, ppp] += item[j];            //V larni massivga yigib oldik
                    }
                    ppp++;
                    break;
                }
            }

            double[,] U_ans = new double[height, height];

            double[,] U_T = new double[height, height];
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    U_T[i, j] = U_matrix[j, i];
                }
            }


            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    int mm = 0;
                    while (mm < height)
                    {
                        U_ans[i, j] += U_matrix[j, mm] * U_T[mm, i];
                        mm++;
                    }
                }
            }



            double[,] V_T = new double[width, width];
            for (int i = 0; i < width; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    V_T[i, j] = V_matrix[j, i];
                }
            }


            double[,] V_ans = new double[r, r];

            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < r; j++)
                {
                    int mm = 0;
                    while (mm < r)
                    {
                        V_ans[i, j] += V_matrix[j, mm] * V_T[mm, i];
                        mm++;
                    }
                }
            }


            Hesh h = new Hesh();                // heshni obyektini yaratamiz
            double result = h.Hesh_Calc(U_T, V_matrix, B);           //resultga h qiymatni qabul qilamiz

            return Convert.ToInt32(result);
        }
    }
}


//var matrix = DenseMatrix.OfArray(A);

//var svd = matrix.Svd();

//var reconstruct = svd.U * svd.W * svd.VT;

//double[,] U = new double[A.GetLength(1), A.GetLength(1)];
//            for (int i = 0; i<A.GetLength(1); i++)
//            {
//                for (int j = 0; j<A.GetLength(1); j++)
//                {
//                    U[i, j] = svd.U[i, j];              // U matritsani topib olamiz
//                }
//            }

//            double[,] V = new double[A.GetLength(1), A.GetLength(1)];
//            for (int i = 0; i<A.GetLength(1); i++)
//            {
//                for (int j = 0; j<A.GetLength(1); j++)
//                {
//                    V[i, j] = svd.VT[j, i];         // Vmatritsani topib olamiz
//                }
//            }
