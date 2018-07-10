﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;

namespace RSA_Project.Controllers
{
    public class Matrix
    {

        public double[,] matrix_A()
        {
            double[,] A = {
               { 149, 18 , 54 , 48 , 42 , 124, 53 , 158 },
                { 82 , 149, 82 , 53 , 124, 81 , 33 , 18  },     //A0 berilgan matritsa
                { 54 , 18 , 33 , 53 , 1  , 7  , 53 , 93  },
                { 39 , 124, 114, 48 , 93 , 53 , 124, 53  },
                { 155, 42 , 51 , 30 , 124, 57 , 82 , 53  },
                { 149, 48 , 42 , 72 , 57 , 53 , 124, 51  },
                { 54 , 48 , 18 , 54 , 27 , 54 , 84 , 82  },
                { 53 , 152, 124, 114, 39 , 42 , 124, 82  }

        };

            return A;
        }

        public double[,] matrix_B()
        {
            double[,] getMatr_B =
               {
                { 0.25 , 0.36 , 0.49 , 0.64 , 0.81 ,   1  , 1.21 , 1.44  },
                { 1.69 , 1.96 , 2.25 , 2.56 , 2.89 , 3.24 , 3.61 ,   4   },     //A0 berilgan matritsa
                { 4.41 , 4.84 , 5.29 , 5.76 , 6.25 , 6.76 , 7.29 , 7.84  },
                { 8.41 ,   9  , 9.61 , 10.24, 10.89, 11.56, 12.25, 12.96 },
                { 13.69, 14.44, 15.21,  16  , 16.81, 17.64, 18.49, 19.36 },
                { 20.25, 21.16, 22.09, 23.04, 24.01,  25  , 26.01, 27.04 },
                { 28.09, 29.16, 30.25, 31.36, 32.49, 33.64, 34.81,  36   },
                { 37.21, 38.44, 39.69, 40.96, 42.25, 43.56, 44.89, 46.24 }
            };

            return getMatr_B;
        }
        //    var matrix = DenseMatrix.OfArray(new double[,] {
        //        { 33.2545, 17.0094, 15.9764, 15.5646 },
        //        { 8.0898 , 31.1955, 0.6476 , 21.5401 },     //A0 berilgan matritsa
        //        { 21.2395, 26.6734, 28.7493, 27.7178 }

        //});
    }
}