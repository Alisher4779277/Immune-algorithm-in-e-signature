using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace RSA_Project.Models.Class
{
    public class Create_matrix
    {
        public double[,] Create(string crypto_text,out double[,] resultB)
        {
            double[,] A = null;                                 // A matritsani e'lon qilish
            double[,] B = null;                                 // B matritsani e'lon qilish

            int size_text = crypto_text.Count();                // Xabarning o`lchamini olish
            bool check = true;

        //pp:

            int temp_number = size_text;                        // vaqtinchalik o`zgaruvchiga xabar o`lchamini olamiz
            int size_array = 0;
            double sqrt_number = 0;
            int EKUB = 1, temp = temp_number, text_num = 0;       // o`zgaruvchilarni e'lon qilib olamiza

            while (check == true) {

                sqrt_number = Math.Sqrt(size_text);
                double numberresult=(sqrt_number % 1);

                if (numberresult == 0)
            {
                    size_array = Convert.ToInt32(Math.Sqrt(size_text));
                    check = false;
            }

                else {
                    size_text++;
                    crypto_text += " "; }
            }
            //p:


            //if ((Math.Abs(EKUB - temp_number)) < 5)
            //{
                A = new double[size_array, size_array];                 // A matritsa yaratiladi

                B = new double[size_array, size_array];        // B matritsa yaratiladi

                for (int i = 0; i < size_array; i++)
                {
                    for (int j = 0; j < size_array; j++)
                    {
                        A[i, j] = crypto_text[text_num];            // maxfiy xabarni raqamlarga o`girib A matritsaga yoziladi
                        text_num++;                                 // maxfiy xabarni 1 taga oshirib boriladi
                        Console.Write(A[i, j] + "   ");
                    }
                    Console.WriteLine("\n");
                }

                //Random rd = new Random();
                //int p=rd.Next(1,10);
                int k = 0;
                int p = 500;

                for (int i = 0; i < size_array; i++)
                {
                    for (int j = 0; j < size_array; j++)
                    {
                        B[i, j] = Math.Pow((p + k), 2) / 100;           // formula yordamida B matritsa yaratiladi
                        k=k+5;                                                             // maxfiy xabar 1 taga oshirib boriladi
                        Console.Write(B[i, j] + "   ");
                    }
                    Console.WriteLine("\n");
                }

            //}
            //else { size_text++; crypto_text += " "; goto pp; }

            resultB = B;
            return A;

        }
    }
}


 //while (EKUB<temp)                                 // matritsanin o`lchami |i-j|<=5 bo`lgunicha ishlaydi
 //           {
 //               for (int i = 2; i <= temp_number;)
 //               {
 //                   if (i != temp_number)
 //                   {
 //                       if ((Math.Abs(EKUB - temp_number)) > 5)     //matritsanin o`lchami |i-j|<=5 bo`lgunicha ishlaydi
 //                       {
 //                           if (temp_number % i == 0)               // kelgan son butun songa bo`linadimi tekshiradi
 //                           {
 //                               temp_number = temp_number / i;
 //                               temp = temp_number;
 //                               EKUB *= i;
 //                           }

 //                           else { i++; }                           // aks holda bo`linuvchini 1 taga oshiradi   
 //                       }
 //                       else { goto p; }                           // aks holda xabarning oxiriga " " qo`shiladi va o`lchami 1 taga oshadi   

 //                   }
 //                   else
 //                   {
 //                       temp = temp_number / i;
 //                       i++;
 //                   }
 //               }
 //           }