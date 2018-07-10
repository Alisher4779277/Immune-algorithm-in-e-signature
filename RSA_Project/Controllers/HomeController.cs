using System;
using System.Collections.Generic;
using System.Linq;
using System.Web;
using System.Web.Mvc;
using RSA_Project.Models.Class;
using System.Numerics;

namespace RSA_Project.Controllers
{
    public class HomeController : Controller
    {

        //numbers num = new numbers();
        //string mess;
        public ActionResult Index()
        {
            return View();
        }

        [HttpPost]
        public ActionResult Index(string message)
        {
            if (message == null || message == "")
            {
                return RedirectToAction("index");
            }
            //if (mess == message)
            //{
            //    return View();
            //}

                Create_matrix cr = new Create_matrix();
                double[,] B;
                double[,] A = cr.Create(message, out B);

                Program prog = new Program();
                BigInteger h = prog.Start(A, B);              // hesh qiymatini topadi

                ViewBag.message = message;
                ViewBag.h = h;
                string hexValue_h = h.ToString("X");
                ViewBag.h_hex = hexValue_h;


            Helper help = new Helper();
                Random rd = new Random();
                List<int> shortprimes = Helper.Primes.GetRange(20, 70);

                int p = shortprimes[rd.Next(shortprimes.Count)];
                ViewBag.p = p;

                int q = shortprimes[rd.Next(shortprimes.Count)];
                ViewBag.q = q;

                BigInteger n = p * q;
                ViewBag.n = n;
                string hexValue_n = n.ToString("X");
                ViewBag.n_hex = hexValue_n;

            int fi = (p - 1) * (q - 1);
                ViewBag.fi = fi;

                int e;
                e = shortprimes.OrderBy(x => Guid.NewGuid()).Where(y => y < fi).FirstOrDefault();
                ViewBag.E = e;
                string hexValue_e = e.ToString("X");
                ViewBag.E_hex = hexValue_e;


            BigInteger d = Helper.MultiplicativeInverse(e, fi);
            ViewBag.D = d;
            string hexValue_d = d.ToString("X");
            ViewBag.d_hex = hexValue_d;

            BigInteger sign = Pow(h, d) % n;
            ViewBag.sign = sign;
            string hexValue_sign = sign.ToString("X");
            ViewBag.sign_hex = hexValue_sign;

            BigInteger h_m = Pow(sign, e) % n;
            ViewBag.h_m = h_m;
            string hexValue_h_m = h_m.ToString("X");
            ViewBag.h_m_hex = hexValue_h_m;

            if (h == h_m)
            {
                ViewBag.verna = "Sizning imzoingiz to`g`ri";
            }
            else
            {
                ViewBag.neverna = "Sizning imzoingiz noto`g`ri";
            }

            //mess = message;
            return View();
        }
        public BigInteger Pow(BigInteger value, BigInteger exponent)
        {
            BigInteger originalValue = value;
            while (exponent-- > 1)
                value = BigInteger.Multiply(value, originalValue);
            return value;
        }
        public ActionResult About()
        {
            ViewBag.Message = "Your application description page.";

            return View();
        }

        public ActionResult Contact()
        {
            ViewBag.Message = "Your contact page.";

            return View();
        }
    }
}











