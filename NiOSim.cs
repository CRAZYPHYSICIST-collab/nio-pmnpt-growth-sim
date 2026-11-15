// Surface growth simulation for NiO wedge films on PMN-PT substrates
// License: MIT (see LICENSE)
//   dh/dtau = a*Lap(h) - b*Lap^2(h) + (c/2)|grad h|^2 + xi ,  <xi xi> = OMEGA * delta(x) delta(tau)
// tau = mean film thickness (nm). a = nu/F [nm], b = kappa/F [nm^3], c = lambda/F [-]
// Noise: shot noise of randomly depositing NiO formula units (OMEGA = a0^3/4).
// Build: Windows -> run.ps1 ; Linux/macOS -> run.sh (Mono)
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Drawing.Text;
using System.Globalization;
using System.IO;
using System.Text;
using System.Threading.Tasks;

namespace NiOSim
{
    public class Fft2
    {
        public int N; int[] rev; double[] ct, st;
        public Fft2(int n)
        {
            N = n; rev = new int[n]; int bits = 0; while ((1 << bits) < n) bits++;
            for (int i = 0; i < n; i++) { int r = 0; for (int b = 0; b < bits; b++) if (((i >> b) & 1) == 1) r |= 1 << (bits - 1 - b); rev[i] = r; }
            ct = new double[n / 2]; st = new double[n / 2];
            for (int k = 0; k < n / 2; k++) { ct[k] = Math.Cos(2 * Math.PI * k / n); st[k] = Math.Sin(2 * Math.PI * k / n); }
        }
        void T1(double[] re, double[] im, bool inv)
        {
            int n = N;
            for (int i = 0; i < n; i++) { int j = rev[i]; if (j > i) { double t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
            for (int size = 2; size <= n; size <<= 1)
            {
                int half = size >> 1, step = n / size;
                for (int i = 0; i < n; i += size)
                    for (int k = 0; k < half; k++)
                    {
                        double wr = ct[k * step], wi = inv ? st[k * step] : -st[k * step];
                        int a = i + k, b = a + half;
                        double tr = re[b] * wr - im[b] * wi, ti = re[b] * wi + im[b] * wr;
                        re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti;
                    }
            }
        }
        public void Run(double[] re, double[] im, bool inv)
        {
            int n = N;
            Parallel.For(0, n, r =>
            {
                double[] br = new double[n], bi = new double[n]; int o = r * n;
                Array.Copy(re, o, br, 0, n); Array.Copy(im, o, bi, 0, n);
                T1(br, bi, inv);
                Array.Copy(br, 0, re, o, n); Array.Copy(bi, 0, im, o, n);
            });
            Parallel.For(0, n, c =>
            {
                double[] br = new double[n], bi = new double[n];
                for (int r = 0; r < n; r++) { br[r] = re[r * n + c]; bi[r] = im[r * n + c]; }
                T1(br, bi, inv);
                for (int r = 0; r < n; r++) { re[r * n + c] = br[r]; im[r * n + c] = bi[r]; }
            });
            if (inv) { double s = 1.0 / ((double)n * n); for (int i = 0; i < n * n; i++) { re[i] *= s; im[i] *= s; } }
        }
    }

    public class Rng
    {
        ulong s; bool has; double spare;
        public Rng(ulong seed) { s = seed * 2685821657736338717UL + 88172645463325252UL; if (s == 0) s = 1; }
        public double U() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return ((s >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
        public double G()
        {
            if (has) { has = false; return spare; }
            double u1 = U(), u2 = U(); double r = Math.Sqrt(-2 * Math.Log(u1));
            spare = r * Math.Sin(2 * Math.PI * u2); has = true; return r * Math.Cos(2 * Math.PI * u2);
        }
    }

    public class Grid
    {
        public int N; public double dx; public double[] L, kmag, kx, ky; public bool[] keep;
        public Grid(int n, double d)
        {
            N = n; dx = d; int m = n * n;
            L = new double[m]; kmag = new double[m]; kx = new double[m]; ky = new double[m]; keep = new bool[m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    int iy = i < n / 2 ? i : i - n, jx = j < n / 2 ? j : j - n;
                    double qy = 2 * Math.PI * iy / (n * d), qx = 2 * Math.PI * jx / (n * d);
                    int id = i * n + j;
                    double sx = Math.Sin(qx * d / 2), sy = Math.Sin(qy * d / 2);
                    L[id] = 4.0 / (d * d) * (sx * sx + sy * sy);
                    kmag[id] = Math.Sqrt(qx * qx + qy * qy);
                    kx[id] = (j == n / 2) ? 0 : qx; ky[id] = (i == n / 2) ? 0 : qy;
                    keep[id] = Math.Abs(iy) < n / 3 && Math.Abs(jx) < n / 3;
                }
        }
        public double Lsys { get { return N * dx; } }
    }

    public class Model
    {
        public string Name; public double a, b, c, W0, xi;
        public Model(string name, double a_, double b_, double c_, double W0_, double xi_) { Name = name; a = a_; b = b_; c = c_; W0 = W0_; xi = xi_; }
        public static double BFromLD(double lD) { return Math.Pow(lD, 4) / 100.0; }   // kappa/F; lD = (b*100nm)^(1/4)
        public static double AFromLnu(double ln) { return ln * ln / 100.0; }         // nu/F;   lnu = sqrt(a*100nm)
    }

    // ---------------------------------------------------------------- charts
    public class Ser
    {
        public string Name; public double[] X, Y, E; public Color Col; public bool Line = true; public int Marker = 0; public float Width = 3f; public DashStyle Dash = DashStyle.Solid;
        public Ser(string n, double[] x, double[] y, Color c) { Name = n; X = x; Y = y; Col = c; }
    }
    public class RefLine { public double V; public string Label; public Color Col; public RefLine(double v, string l, Color c) { V = v; Label = l; Col = c; } }

    public class Chart
    {
        public string Title = "", XLabel = "", YLabel = ""; public bool LogX, LogY; public double X0, X1, Y0, Y1;
        public List<Ser> S = new List<Ser>(); public List<RefLine> H = new List<RefLine>(); public List<RefLine> V = new List<RefLine>();
        public int LegendPos = 0; public string Note = null;
        float ml = 150, mr = 40, mt = 80, mb = 115, pw, ph;
        double TX(double x) { return LogX ? Math.Log10(x) : x; }
        double TY(double y) { return LogY ? Math.Log10(y) : y; }
        float PX(double x) { return ml + (float)((TX(x) - TX(X0)) / (TX(X1) - TX(X0)) * pw); }
        float PY(double y) { return mt + ph - (float)((TY(y) - TY(Y0)) / (TY(Y1) - TY(Y0)) * ph); }
        bool Ok(double x, double y) { return !double.IsNaN(x) && !double.IsNaN(y) && !double.IsInfinity(x) && !double.IsInfinity(y) && (!LogX || x > 0) && (!LogY || y > 0); }
        static string Fmt(double v) { return v.ToString("0.####", CultureInfo.InvariantCulture); }
        static List<double[]> Ticks(double a, double b, bool log)
        {
            List<double[]> t = new List<double[]>();
            if (log)
            {
                int e0 = (int)Math.Floor(Math.Log10(a)), e1 = (int)Math.Ceiling(Math.Log10(b));
                bool fewDec = Math.Log10(b / a) < 1.6;
                for (int e = e0; e <= e1; e++)
                    for (int m = 1; m <= 9; m++)
                    {
                        double v = m * Math.Pow(10, e);
                        if (v < a * 0.9999 || v > b * 1.0001) continue;
                        bool major = m == 1 || (fewDec && (m == 2 || m == 5));
                        t.Add(new double[] { v, major ? 1 : 0 });
                    }
            }
            else
            {
                double raw = (b - a) / 6, p = Math.Pow(10, Math.Floor(Math.Log10(raw))), st = p;
                if (raw / p > 5) st = 10 * p; else if (raw / p > 2) st = 5 * p; else if (raw / p > 1) st = 2 * p;
                for (double v = Math.Ceiling(a / st - 1e-9) * st; v <= b + 1e-9 * st; v += st) t.Add(new double[] { Math.Abs(v) < 1e-12 ? 0 : v, 1 });
            }
            return t;
        }
        public void Save(string path)
        {
            int W = 1500, Hh = 1000; pw = W - ml - mr; ph = Hh - mt - mb;
            using (Bitmap bmp = new Bitmap(W, Hh))
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias; g.TextRenderingHint = TextRenderingHint.AntiAlias; g.Clear(Color.White);
                Font ft = new Font("Segoe UI", 30f, FontStyle.Bold, GraphicsUnit.Pixel);
                Font fl = new Font("Segoe UI", 28f, FontStyle.Regular, GraphicsUnit.Pixel);
                Font fk = new Font("Segoe UI", 24f, FontStyle.Regular, GraphicsUnit.Pixel);
                Font fg = new Font("Segoe UI", 21f, FontStyle.Regular, GraphicsUnit.Pixel);
                Pen gMaj = new Pen(Color.FromArgb(210, 210, 210), 1.5f), gMin = new Pen(Color.FromArgb(236, 236, 236), 1f);
                foreach (double[] t in Ticks(X0, X1, LogX))
                {
                    float x = PX(t[0]); g.DrawLine(t[1] > 0 ? gMaj : gMin, x, mt, x, mt + ph);
                    g.DrawLine(Pens.Black, x, mt + ph, x, mt + ph + (t[1] > 0 ? 10 : 5));
                    if (t[1] > 0) { string s = Fmt(t[0]); SizeF z = g.MeasureString(s, fk); g.DrawString(s, fk, Brushes.Black, x - z.Width / 2, mt + ph + 12); }
                }
                foreach (double[] t in Ticks(Y0, Y1, LogY))
                {
                    float y = PY(t[0]); g.DrawLine(t[1] > 0 ? gMaj : gMin, ml, y, ml + pw, y);
                    g.DrawLine(Pens.Black, ml - (t[1] > 0 ? 10 : 5), y, ml, y);
                    if (t[1] > 0) { string s = Fmt(t[0]); SizeF z = g.MeasureString(s, fk); g.DrawString(s, fk, Brushes.Black, ml - 14 - z.Width, y - z.Height / 2); }
                }
                g.SetClip(new RectangleF(ml, mt, pw, ph));
                foreach (RefLine r in H)
                {
                    Pen p = new Pen(r.Col, 2.5f); p.DashStyle = DashStyle.Dash; float y = PY(r.V); g.DrawLine(p, ml, y, ml + pw, y);
                    if (!string.IsNullOrEmpty(r.Label)) { SizeF z = g.MeasureString(r.Label, fg); g.DrawString(r.Label, fg, new SolidBrush(r.Col), ml + pw - z.Width - 8, y - z.Height - 2); }
                }
                foreach (RefLine r in V)
                {
                    Pen p = new Pen(r.Col, 2.5f); p.DashStyle = DashStyle.Dot; float x = PX(r.V); g.DrawLine(p, x, mt, x, mt + ph);
                    if (!string.IsNullOrEmpty(r.Label)) g.DrawString(r.Label, fg, new SolidBrush(r.Col), x + 4, mt + 6);
                }
                foreach (Ser s in S)
                {
                    Pen p = new Pen(s.Col, s.Width); p.DashStyle = s.Dash; p.LineJoin = LineJoin.Round;
                    if (s.Line)
                    {
                        List<PointF> pts = new List<PointF>();
                        for (int i = 0; i < s.X.Length; i++)
                        {
                            if (Ok(s.X[i], s.Y[i])) pts.Add(new PointF(PX(s.X[i]), PY(s.Y[i])));
                            else { if (pts.Count > 1) g.DrawLines(p, pts.ToArray()); pts.Clear(); }
                        }
                        if (pts.Count > 1) g.DrawLines(p, pts.ToArray());
                    }
                    if (s.E != null)
                    {
                        Pen pe = new Pen(s.Col, 2f);
                        for (int i = 0; i < s.X.Length; i++)
                        {
                            double lo = s.Y[i] - s.E[i], hi = s.Y[i] + s.E[i]; if (LogY && lo <= 0) lo = s.Y[i] * 0.5;
                            if (!Ok(s.X[i], s.Y[i])) continue; float x = PX(s.X[i]);
                            g.DrawLine(pe, x, PY(lo), x, PY(hi)); g.DrawLine(pe, x - 7, PY(lo), x + 7, PY(lo)); g.DrawLine(pe, x - 7, PY(hi), x + 7, PY(hi));
                        }
                    }
                    if (s.Marker > 0)
                        for (int i = 0; i < s.X.Length; i++)
                        {
                            if (!Ok(s.X[i], s.Y[i])) continue; DrawMarker(g, s.Marker, PX(s.X[i]), PY(s.Y[i]), s.Col);
                        }
                }
                g.ResetClip();
                g.DrawRectangle(new Pen(Color.Black, 2f), ml, mt, pw, ph);
                SizeF tz = g.MeasureString(Title, ft); g.DrawString(Title, ft, Brushes.Black, ml + pw / 2 - tz.Width / 2, 20);
                SizeF xz = g.MeasureString(XLabel, fl); g.DrawString(XLabel, fl, Brushes.Black, ml + pw / 2 - xz.Width / 2, mt + ph + 55);
                GraphicsState gs = g.Save(); g.TranslateTransform(22, mt + ph / 2); g.RotateTransform(-90);
                SizeF yz = g.MeasureString(YLabel, fl); g.DrawString(YLabel, fl, Brushes.Black, -yz.Width / 2, 0); g.Restore(gs);
                // legend
                List<Ser> ls = new List<Ser>(); foreach (Ser s in S) if (!string.IsNullOrEmpty(s.Name)) ls.Add(s);
                if (ls.Count > 0)
                {
                    float lw = 0, lh = 34; foreach (Ser s in ls) lw = Math.Max(lw, g.MeasureString(s.Name, fg).Width);
                    float bw = lw + 90, bh = lh * ls.Count + 14;
                    float bx = (LegendPos == 1 || LegendPos == 2) ? ml + pw - bw - 14 : ml + 14;
                    float by = (LegendPos == 2 || LegendPos == 3) ? mt + ph - bh - 14 : mt + 14;
                    g.FillRectangle(new SolidBrush(Color.FromArgb(235, 255, 255, 255)), bx, by, bw, bh); g.DrawRectangle(Pens.Gray, bx, by, bw, bh);
                    for (int i = 0; i < ls.Count; i++)
                    {
                        Ser s = ls[i]; float yy = by + 7 + lh * i + lh / 2;
                        if (s.Line) { Pen p = new Pen(s.Col, s.Width); p.DashStyle = s.Dash; g.DrawLine(p, bx + 10, yy, bx + 62, yy); }
                        if (s.Marker > 0) DrawMarker(g, s.Marker, bx + 36, yy, s.Col);
                        g.DrawString(s.Name, fg, Brushes.Black, bx + 74, yy - 14);
                    }
                }
                if (Note != null) { SizeF nz = g.MeasureString(Note, fg); g.DrawString(Note, fg, Brushes.DimGray, ml + pw - nz.Width, Hh - 34); }
                bmp.Save(path, ImageFormat.Png);
            }
        }
        static void DrawMarker(Graphics g, int m, float x, float y, Color c)
        {
            float r = 9; Brush b = new SolidBrush(c); Pen w = new Pen(Color.White, 1.5f);
            if (m == 1) { g.FillEllipse(b, x - r, y - r, 2 * r, 2 * r); g.DrawEllipse(w, x - r, y - r, 2 * r, 2 * r); }
            else if (m == 2) { g.FillRectangle(b, x - r, y - r, 2 * r, 2 * r); g.DrawRectangle(w, x - r, y - r, 2 * r, 2 * r); }
            else { PointF[] p = { new PointF(x, y - r - 2), new PointF(x + r + 1, y + r), new PointF(x - r - 1, y + r) }; g.FillPolygon(b, p); g.DrawPolygon(w, p); }
        }
    }

    // ---------------------------------------------------------------- main program
    public static class Program
    {
        public const double OMEGA = 0.01822;      // nm^3, volume of one NiO formula unit, a0^3/4 (a0 = 0.4177 nm)
        public static readonly CultureInfo IC = CultureInfo.InvariantCulture;
        public static readonly double[] T6 = { 10, 16, 27, 44, 73, 120 };   // six AFM sampling points along the wedge (nm)
        // values collected for the summary table
        static double[] Tctrl = new double[3], Tpix = new double[3];
        static double TbW, TbWsd, TbS, TbSsd, TbW2, TbW2sd, TbS2, TbS2sd, Rid, Rsu, Kstar, NLmin = 9, NLmax = -9;
        static bool haveExact, haveVirtual, haveNL;
        public const double TIP_R = 10.0;          // AFM tip radius (nm)
        static string Out; static StringBuilder Sum = new StringBuilder();
        static readonly Color CW = Color.FromArgb(214, 39, 40), CS = Color.FromArgb(31, 119, 180), CG = Color.FromArgb(44, 160, 44), CO = Color.FromArgb(255, 127, 14), CP = Color.FromArgb(148, 103, 189);

        static void Log(string s) { Console.WriteLine(s); Sum.AppendLine(s); }
        static string F(double v) { return v.ToString("0.####", IC); }
        static string F2(double v) { return v.ToString("0.00", IC); }
        static string F3(double v) { return v.ToString("0.000", IC); }

        public static void Main(string[] args)
        {
            string mode = args[0]; Out = args[1]; Directory.CreateDirectory(Out);
            System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew();
            if (mode == "test") Test();
            if (mode == "exact" || mode == "all") Exact();
            if (mode == "virtual" || mode == "all") Virtual();
            if (mode == "nonlinear" || mode == "all") Nonlinear();
            if (mode == "all") SummaryTable();
            Log("Total time: " + F(sw.Elapsed.TotalSeconds) + " s");
            File.WriteAllText(Path.Combine(Out, "summary_" + mode + ".txt"), Sum.ToString(), new UTF8Encoding(true));
        }

        // ---------- basic functions
        public static double PsdExact(Model m, double L, double k, double tau)
        {
            double w = m.a * L + m.b * L * L, x = 2 * w * tau, gg;
            if (x < 1e-6) gg = tau * (1 - x / 2 + x * x / 6); else gg = (1 - Math.Exp(-x)) / (2 * w);
            double p = OMEGA * gg;
            if (m.W0 > 0) p += Math.PI * m.xi * m.xi * m.W0 * m.W0 * Math.Exp(-k * k * m.xi * m.xi / 4 - x);
            return p;   // nm^4
        }
        public static double WExact(Grid g, Model m, double tau)
        {
            int n = g.N * g.N; double s = 0;
            for (int i = 1; i < n; i++) s += PsdExact(m, g.L[i], g.kmag[i], tau);
            return Math.Sqrt(s / ((double)n * g.dx * g.dx));
        }
        static double[] Sample(Grid g, Fft2 f, Rng r, Model m, double tau)
        {
            int n = g.N * g.N; double[] re = new double[n], im = new double[n];
            for (int i = 0; i < n; i++) re[i] = r.G();
            f.Run(re, im, false); re[0] = 0; im[0] = 0;
            for (int i = 1; i < n; i++) { double s = Math.Sqrt(PsdExact(m, g.L[i], g.kmag[i], tau) / (g.dx * g.dx)); re[i] *= s; im[i] *= s; }
            f.Run(re, im, true);
            return re;
        }
        static double Rms(double[] u)
        {
            double m = 0; foreach (double v in u) m += v; m /= u.Length;
            double s = 0; foreach (double v in u) s += (v - m) * (v - m); return Math.Sqrt(s / u.Length);
        }
        static double[] Tip(double[] u, int N, double dx, double R)
        {
            int rp = (int)Math.Ceiling(12.0 / dx); double[] o = new double[N * N];
            List<int[]> off = new List<int[]>(); List<double> pen = new List<double>();
            for (int di = -rp; di <= rp; di++) for (int dj = -rp; dj <= rp; dj++) if (di * di + dj * dj <= rp * rp) { off.Add(new int[] { di, dj }); pen.Add((di * di + dj * dj) * dx * dx / (2 * R)); }
            int[][] offs = off.ToArray(); double[] pens = pen.ToArray();
            Parallel.For(0, N, i =>
            {
                for (int j = 0; j < N; j++)
                {
                    double best = double.NegativeInfinity;
                    for (int q = 0; q < offs.Length; q++)
                    {
                        int ii = (i + offs[q][0] + N) % N, jj = (j + offs[q][1] + N) % N;
                        double v = u[ii * N + jj] - pens[q]; if (v > best) best = v;
                    }
                    o[i * N + j] = best;
                }
            });
            return o;
        }
        // radially averaged PSD (nm^4): bin width 2pi/L
        static double[][] RadialPsdFromField(Grid g, Fft2 f, double[] u)
        {
            int n = g.N * g.N; double[] re = (double[])u.Clone(), im = new double[n];
            double m = 0; foreach (double v in re) m += v; m /= n; for (int i = 0; i < n; i++) re[i] -= m;
            f.Run(re, im, false);
            double[] p = new double[n]; for (int i = 0; i < n; i++) p[i] = g.dx * g.dx * (re[i] * re[i] + im[i] * im[i]) / n;
            return Bin(g, p);
        }
        static double[][] RadialPsdExact(Grid g, Model m, double tau)
        {
            int n = g.N * g.N; double[] p = new double[n]; for (int i = 0; i < n; i++) p[i] = PsdExact(m, g.L[i], g.kmag[i], tau);
            return Bin(g, p);
        }
        static double[][] Bin(Grid g, double[] p)
        {
            int nb = g.N / 2; double dk = 2 * Math.PI / g.Lsys; double[] s = new double[nb + 1]; int[] c = new int[nb + 1];
            for (int i = 1; i < p.Length; i++) { int b = (int)Math.Round(g.kmag[i] / dk); if (b >= 1 && b <= nb) { s[b] += p[i]; c[b]++; } }
            List<double> kk = new List<double>(), pp = new List<double>();
            for (int b = 1; b <= nb; b++) if (c[b] > 0) { kk.Add(b * dk); pp.Add(s[b] / c[b]); }
            return new double[][] { kk.ToArray(), pp.ToArray() };
        }
        // weighted linear fit: [slope, standard error, intercept]
        static double[] Fit(double[] x, double[] y, double[] sy)
        {
            int n = x.Length; double Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
            for (int i = 0; i < n; i++) { double w = sy == null ? 1 : 1 / (sy[i] * sy[i]); Sw += w; Sx += w * x[i]; Sy += w * y[i]; Sxx += w * x[i] * x[i]; Sxy += w * x[i] * y[i]; }
            double D = Sw * Sxx - Sx * Sx, sl = (Sw * Sxy - Sx * Sy) / D, ic = (Sxx * Sy - Sx * Sxy) / D, chi = 0;
            for (int i = 0; i < n; i++) { double w = sy == null ? 1 : 1 / (sy[i] * sy[i]); double r = y[i] - ic - sl * x[i]; chi += w * r * r; }
            double red = chi / (n - 2), se = Math.Sqrt(Sw / D) * (sy == null ? Math.Sqrt(red) : Math.Sqrt(Math.Max(1, red)));
            return new double[] { sl, se, ic };
        }
        static double BetaExact(Grid g, Model m)
        {
            double[] x = new double[T6.Length], y = new double[T6.Length];
            for (int i = 0; i < T6.Length; i++) { x[i] = Math.Log(T6[i]); y[i] = Math.Log(WExact(g, m, T6[i])); }
            return Fit(x, y, null)[0];
        }
        static double[] LogSpace(double a, double b, int n) { double[] r = new double[n]; for (int i = 0; i < n; i++) r[i] = a * Math.Pow(b / a, i / (double)(n - 1)); return r; }
        static void Csv(string name, string header, List<string> rows) { StringBuilder sb = new StringBuilder(); sb.AppendLine(header); foreach (string r in rows) sb.AppendLine(r); File.WriteAllText(Path.Combine(Out, name), sb.ToString(), new UTF8Encoding(false)); }

        // ---------- test: FFT and sampling verification
        // ---------- summary table: key results of this run
        static void SummaryTable()
        {
            if (!(haveExact && haveVirtual && haveNL)) return;
            Log("");
            Log("=== SUMMARY TABLE ===");
            Log("Boundary cases, beta (RD / nu only / kappa only): " + F2(Tctrl[0]) + " / " + F2(Tctrl[1]) + " / " + F2(Tctrl[2]) + "  (theory: 1/2, 0, 1/4)");
            Log("Effective beta, weak / strong diffusion (10-120 nm, 8 repeats): " + F2(TbW) + " +- " + F2(TbWsd) + " / " + F3(TbS) + " +- " + F3(TbSsd));
            Log("Same, with 0.2 nm substrate roughness + AFM tip: " + F2(TbW2) + " +- " + F2(TbW2sd) + " / " + F3(TbS2) + " +- " + F3(TbS2sd));
            Log("Roughness ratio at 100 nm, weak / strong: " + F2(Rid) + " (ideal) - " + F2(Rsu) + " (substrate + tip)");
            Log("Spatial frequency where the PSDs separate: " + F2(Kstar) + " nm^-1 (tip limit ~ " + F2(2 * Math.PI / TIP_R) + " nm^-1)");
            Log("Effect of the nonlinear term, strong diffusion (lambda/F = -1...3): beta = " + F2(NLmin) + " - " + F2(NLmax));
            Log("Pixel-size sensitivity of beta, weak diffusion (1 / 2 / 4 nm): " + F2(Tpix[0]) + " / " + F2(Tpix[1]) + " / " + F2(Tpix[2]));
        }

        static void Test()
        {
            Grid g = new Grid(256, 2.0); Fft2 f = new Fft2(256); Rng r = new Rng(1);
            int n = 256 * 256; double[] a = new double[n], b = new double[n], a0 = new double[n];
            for (int i = 0; i < n; i++) { a[i] = r.G(); a0[i] = a[i]; }
            System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew();
            f.Run(a, b, false); f.Run(a, b, true);
            double err = 0; for (int i = 0; i < n; i++) err = Math.Max(err, Math.Abs(a[i] - a0[i]));
            Log("FFT round-trip max error: " + err.ToString("E2", IC) + "  (2 FFT 256^2: " + sw.ElapsedMilliseconds + " ms)");
            Model[] ms = { new Model("RD", 0, 0, 0, 0, 0), new Model("MH", 0, Model.BFromLD(10), 0, 0, 0), new Model("EW", Model.AFromLnu(10), 0, 0, 0, 0) };
            foreach (Model m in ms)
            {
                double we = WExact(g, m, 100), ws = 0; int K = 20; for (int k = 0; k < K; k++) ws += Rms(Sample(g, f, r, m, 100)) / K;
                Log(m.Name + ": W_exact(100nm)=" + F3(we) + "  W_sampled=" + F3(ws) + "   RD analytic=" + F3(Math.Sqrt(100 * OMEGA / 4.0)));
            }
            Model mh = new Model("MH20", 0, Model.BFromLD(20), 0, 0, 0);
            double cont = Math.Sqrt(OMEGA * Math.Sqrt(Math.PI) / (4 * Math.PI) * Math.Sqrt(100 / (2 * mh.b)));
            Log("MH lD=20: W_exact=" + F3(WExact(g, mh, 100)) + "  continuum analytic=" + F3(cont));
        }

        // ---------- 1) exact (linear) solution and parameter scan
        static void Exact()
        {
            Grid g = new Grid(512, 2.0);
            Log("=== EXACT LINEAR SOLUTION (N=512, dx=2 nm, L=1024 nm) ===");
            Log("Noise: OMEGA=" + F(OMEGA) + " nm^3; RD limit W=sqrt(tau*OMEGA/dx^2) -> at 100 nm: " + F3(Math.Sqrt(100 * OMEGA / 4)) + " nm");
            double[] lDs = LogSpace(0.5, 64, 29); double[] lnus = { 0, 2, 5, 10, 20 };
            int NL = lDs.Length, NN = lnus.Length;
            double[,] beta = new double[NN + 1, NL], w60 = new double[NN + 1, NL], w100 = new double[NN + 1, NL];
            Parallel.For(0, (NN + 1) * NL, id =>
            {
                int iv = id / NL, il = id % NL;
                Model m = iv < NN ? new Model("", Model.AFromLnu(lnus[iv]), Model.BFromLD(lDs[il]), 0, 0, 0)
                                  : new Model("", 0, Model.BFromLD(lDs[il]), 0, 0.2, 20);
                beta[iv, il] = BetaExact(g, m); w60[iv, il] = WExact(g, m, 60); w100[iv, il] = WExact(g, m, 100);
            });
            List<string> rows = new List<string>();
            for (int iv = 0; iv <= NN; iv++) for (int il = 0; il < NL; il++)
                    rows.Add((iv < NN ? F(lnus[iv]) : "0+substrate") + "," + F(lDs[il]) + "," + F(Model.BFromLD(lDs[il])) + "," + F3(beta[iv, il]) + "," + F(w60[iv, il]) + "," + F(w100[iv, il]));
            Csv("scan.csv", "l_nu_nm,l_D_nm,b_nm3,beta_10_120,W60_nm,W100_nm", rows);

            // boundary cases: show that the code reproduces the theoretical growth exponents (Barabasi & Stanley 1995)
            Model[] ctrl = { new Model("RD", 0, 0, 0, 0, 0), new Model("EW", Model.AFromLnu(20), 0, 0, 0, 0), new Model("MH", 0, Model.BFromLD(20), 0, 0, 0) };
            for (int q = 0; q < 3; q++) Tctrl[q] = BetaExact(g, ctrl[q]);
            Log("--- Boundary cases (beta, 10-120 nm fit) ---");
            Log("random deposition only: " + F3(Tctrl[0]) + " (theory 1/2) | surface tension only, l_nu=20 nm: " + F3(Tctrl[1]) + " (theory 0, logarithmic) | surface diffusion only, l_D=20 nm: " + F3(Tctrl[2]) + " (theory 1/4)");
            haveExact = true;

            Log("--- beta (10-120 nm fit) and W100 ---");
            for (int iv = 0; iv <= NN; iv++)
            {
                string nm = iv < NN ? "l_nu=" + F(lnus[iv]) + " nm" : "l_nu=0 + substrate W0=0.2nm";
                StringBuilder sb = new StringBuilder(nm + ": ");
                for (int il = 0; il < NL; il += 4) sb.Append("[lD=" + F2(lDs[il]) + " b=" + F3(beta[iv, il]) + " W=" + F3(w100[iv, il]) + "] ");
                Log(sb.ToString());
            }

            // resolution (pixel size) sensitivity: same physics, different dx (L fixed = 1024 nm)
            Log("--- Pixel-size sensitivity (nu=0, L=1024 nm) ---");
            foreach (double ld in new double[] { 1, 10 })
            {
                StringBuilder sb = new StringBuilder("lD=" + F(ld) + " nm: ");
                foreach (int nn in new int[] { 1024, 512, 256 })
                {
                    Grid gg = new Grid(nn, 1024.0 / nn); Model m = new Model("", 0, Model.BFromLD(ld), 0, 0, 0);
                    double bb = BetaExact(gg, m); if (ld == 1) Tpix[nn == 1024 ? 0 : nn == 512 ? 1 : 2] = bb;
                    sb.Append("[dx=" + F(gg.dx) + " nm: beta=" + F3(bb) + ", W100=" + F3(WExact(gg, m, 100)) + "] ");
                }
                Log(sb.ToString());
            }

            // Figures 3 and 4
            Color[] cols = { Color.Black, CO, CG, CS, CP, CW };
            Chart c3 = new Chart(); c3.Title = "Effective growth exponent β (10–120 nm, 6-point fit)"; c3.XLabel = "Surface diffusion length ℓ_D = (κ/F · 100 nm)^1/4  [nm]"; c3.YLabel = "β_eff";
            c3.LogX = true; c3.X0 = 0.5; c3.X1 = 64; c3.Y0 = 0; c3.Y1 = 0.55; c3.LegendPos = 3;
            Chart c4 = new Chart(); c4.Title = "RMS roughness at 100 nm thickness (shot noise, 2 nm pixel)"; c4.XLabel = c3.XLabel; c4.YLabel = "W(100 nm)  [nm]";
            c4.LogX = true; c4.LogY = true; c4.X0 = 0.5; c4.X1 = 64; c4.Y0 = 0.005; c4.Y1 = 5; c4.LegendPos = 3;
            for (int iv = 0; iv <= NN; iv++)
            {
                double[] yb = new double[NL], yw = new double[NL]; for (int il = 0; il < NL; il++) { yb[il] = beta[iv, il]; yw[il] = w100[iv, il]; }
                string nm = iv < NN ? (lnus[iv] == 0 ? "ν = 0 (pure Mullins–Herring)" : "ℓ_ν = " + F(lnus[iv]) + " nm") : "ν = 0 + substrate roughness 0.2 nm";
                Ser s = new Ser(nm, lDs, yb, cols[iv]); if (iv == NN) s.Dash = DashStyle.Dash; c3.S.Add(s);
                Ser s2 = new Ser(nm, lDs, yw, cols[iv]); if (iv == NN) s2.Dash = DashStyle.Dash; c4.S.Add(s2);
            }
            c3.H.Add(new RefLine(0.5, "random deposition 0.5", Color.Silver)); c3.H.Add(new RefLine(0.25, "MH asymptote 0.25", Color.Silver));
            c3.V.Add(new RefLine(2.0, "AFM pixel (2 nm)", Color.DimGray)); c4.V.Add(new RefLine(2.0, "AFM pixel (2 nm)", Color.DimGray));
            c4.Note = "AFM noise floor ≈ 0.03–0.05 nm";
            c3.Save(Path.Combine(Out, "fig3_beta_scan.png")); c4.Save(Path.Combine(Out, "fig4_W100_scan.png"));
        }

        // ---------- 2) virtual experiment: 6 thicknesses x 3 regions, tip and substrate effects
        static void Virtual()
        {
            Grid g = new Grid(512, 2.0); Fft2 f = new Fft2(512); int N = 512;
            double lw = 1, ls = 10; int M = 8;
            Log("=== VIRTUAL EXPERIMENT (N=512, dx=2 nm, equivalent to a 1x1 um AFM image; weak lD=" + F(lw) + " nm, strong lD=" + F(ls) + " nm; nu=0) ===");
            Log("Arrhenius equivalent: kappa ratio = (lD ratio)^4 = " + F(Math.Pow(ls / lw, 4)) + " ; between 25 and 300 C, Ed ~ " + F2(Math.Log(Math.Pow(ls / lw, 4)) / ((1 / 298.15 - 1 / 573.15) / 8.617e-5)) + " eV (1/T prefactor and F neglected)");
            Model[] base_ = { new Model("weak", 0, Model.BFromLD(lw), 0, 0, 0), new Model("strong", 0, Model.BFromLD(ls), 0, 0, 0) };
            Model[] subst = { new Model("weak+substrate", 0, Model.BFromLD(lw), 0, 0.2, 20), new Model("strong+substrate", 0, Model.BFromLD(ls), 0, 0.2, 20) };
            string[] vn = { "ideal", "tip (R=10nm)", "substrate (0.2nm)", "substrate+tip" };
            double[,] bMean = new double[2, 4], bSd = new double[2, 4], bSe = new double[2, 4];
            double[,,] Wm = new double[2, 4, 6], We = new double[2, 4, 6];
            double[][][] psdMeas = new double[2][][], psdTip = new double[2][][];
            double[][] maps = new double[2][];
            double[,] w100v = new double[2, 4], w60v = new double[2, 4];
            List<string> rows = new List<string>();
            for (int ci = 0; ci < 2; ci++)
                for (int sv = 0; sv < 2; sv++)
                {
                    Model m = sv == 0 ? base_[ci] : subst[ci];
                    Rng r = new Rng((ulong)(1000 + 10 * ci + sv));
                    double[][] betas = { new double[M], new double[M] }, ses = { new double[M], new double[M] };
                    double[,,] Wacc = new double[2, 6, 2];
                    for (int rep = 0; rep < M; rep++)
                    {
                        double[,] mean = new double[2, 6], sd = new double[2, 6];
                        for (int ti = 0; ti < 6; ti++)
                        {
                            double[,] w = new double[2, 3];
                            for (int reg = 0; reg < 3; reg++)
                            {
                                double[] u = Sample(g, f, r, m, T6[ti]); double[] ut = Tip(u, N, g.dx, TIP_R);
                                w[0, reg] = Rms(u); w[1, reg] = Rms(ut);
                                if (rep == 0 && reg == 0 && ti == 5 && sv == 0) { psdMeas[ci] = RadialPsdFromField(g, f, u); psdTip[ci] = RadialPsdFromField(g, f, ut); }
                            }
                            for (int t = 0; t < 2; t++)
                            {
                                double mu = (w[t, 0] + w[t, 1] + w[t, 2]) / 3, v = 0; for (int q = 0; q < 3; q++) v += (w[t, q] - mu) * (w[t, q] - mu);
                                mean[t, ti] = mu; sd[t, ti] = Math.Sqrt(v / 2); Wacc[t, ti, 0] += mu / M; Wacc[t, ti, 1] += mu * mu / M;
                            }
                        }
                        for (int t = 0; t < 2; t++)
                        {
                            double[] x = new double[6], y = new double[6], sy = new double[6];
                            for (int ti = 0; ti < 6; ti++) { x[ti] = Math.Log(T6[ti]); y[ti] = Math.Log(mean[t, ti]); sy[ti] = Math.Max(sd[t, ti] / Math.Sqrt(3), 1e-4 * mean[t, ti]) / mean[t, ti]; }
                            double[] ft = Fit(x, y, sy); betas[t][rep] = ft[0]; ses[t][rep] = ft[1];
                        }
                    }
                    for (int t = 0; t < 2; t++)
                    {
                        int vi = sv * 2 + t; double mu = 0, s2 = 0, se = 0;
                        for (int k = 0; k < M; k++) { mu += betas[t][k] / M; se += ses[t][k] / M; }
                        for (int k = 0; k < M; k++) s2 += (betas[t][k] - mu) * (betas[t][k] - mu);
                        bMean[ci, vi] = mu; bSd[ci, vi] = Math.Sqrt(s2 / (M - 1)); bSe[ci, vi] = se;
                        for (int ti = 0; ti < 6; ti++) { Wm[ci, vi, ti] = Wacc[t, ti, 0]; We[ci, vi, ti] = Math.Sqrt(Math.Max(0, Wacc[t, ti, 1] - Wacc[t, ti, 0] * Wacc[t, ti, 0])); }
                    }
                    // 60 and 100 nm
                    double[] u60 = Sample(g, f, r, m, 60), u100 = Sample(g, f, r, m, 100);
                    w60v[ci, sv * 2] = Rms(u60); w60v[ci, sv * 2 + 1] = Rms(Tip(u60, N, g.dx, TIP_R));
                    double[] u100t = Tip(u100, N, g.dx, TIP_R);
                    w100v[ci, sv * 2] = Rms(u100); w100v[ci, sv * 2 + 1] = Rms(u100t);
                    if (sv == 0) maps[ci] = u100;
                }
            string[] cn = { "weak", "strong" };
            Log("--- beta (mean +- SD across repeats; SE = mean fit standard error) ---");
            for (int vi = 0; vi < 4; vi++)
                Log(vn[vi].PadRight(18) + " weak: " + F3(bMean[0, vi]) + " +- " + F3(bSd[0, vi]) + " (SE " + F3(bSe[0, vi]) + ")   strong: " + F3(bMean[1, vi]) + " +- " + F3(bSd[1, vi]) + " (SE " + F3(bSe[1, vi]) + ")");
            Log("--- W(60 nm) / W(100 nm) [nm] ---");
            for (int vi = 0; vi < 4; vi++)
                Log(vn[vi].PadRight(18) + " weak: " + F3(w60v[0, vi]) + " / " + F3(w100v[0, vi]) + "   strong: " + F3(w60v[1, vi]) + " / " + F3(w100v[1, vi]) + "   ratio(100nm): " + F2(w100v[0, vi] / w100v[1, vi]));
            TbW = bMean[0, 0]; TbWsd = bSd[0, 0]; TbS = bMean[1, 0]; TbSsd = bSd[1, 0];
            TbW2 = bMean[0, 3]; TbW2sd = bSd[0, 3]; TbS2 = bMean[1, 3]; TbS2sd = bSd[1, 3];
            Rid = w100v[0, 0] / w100v[1, 0]; Rsu = w100v[0, 3] / w100v[1, 3];
            {   // PSD separation frequency: first k where the ratio of the noise-free exact PSDs at 120 nm exceeds 3
                double[][] pz = RadialPsdExact(g, base_[0], 120), pg = RadialPsdExact(g, base_[1], 120); Kstar = double.NaN;
                for (int i = 0; i < pz[0].Length; i++) if (pz[1][i] / pg[1][i] > 3) { Kstar = pz[0][i]; break; }
                Log("PSD separation frequency (120 nm, weak/strong > 3): k* = " + F3(Kstar) + " nm^-1 ; tip limit ~ 2pi/R = " + F2(2 * Math.PI / TIP_R) + " nm^-1");
            }
            haveVirtual = true;
            Log("Exact solution W100: weak=" + F3(WExact(g, base_[0], 100)) + " strong=" + F3(WExact(g, base_[1], 100)) + "; W60: weak=" + F3(WExact(g, base_[0], 60)) + " strong=" + F3(WExact(g, base_[1], 60)));
            for (int ci = 0; ci < 2; ci++) for (int vi = 0; vi < 4; vi++) for (int ti = 0; ti < 6; ti++)
                        rows.Add(cn[ci] + "," + vn[vi] + "," + F(T6[ti]) + "," + F(Wm[ci, vi, ti]) + "," + F(We[ci, vi, ti]));
            Csv("virtual_experiment_W.csv", "state,variant,thickness_nm,W_mean_nm,W_sd_nm", rows);

            // Figure 1: W(tau)
            double[] td = LogSpace(5, 150, 50);
            Chart c1 = new Chart(); c1.Title = "Roughness vs. thickness (linear model, shot noise)"; c1.XLabel = "Film thickness h  [nm]"; c1.YLabel = "RMS roughness W  [nm]";
            c1.LogX = true; c1.LogY = true; c1.X0 = 5; c1.X1 = 150; c1.Y0 = 0.01; c1.Y1 = 2; c1.LegendPos = 3;
            Model[] all = { base_[0], base_[1], subst[0], subst[1] }; Color[] cc = { CW, CS, CW, CS };
            string[] nmL = { "weak (ℓ_D=1 nm), exact solution", "strong (ℓ_D=10 nm), exact solution", "weak + substrate 0.2 nm", "strong + substrate 0.2 nm" };
            for (int k = 0; k < 4; k++)
            {
                double[] y = new double[td.Length]; for (int i = 0; i < td.Length; i++) y[i] = WExact(g, all[k], td[i]);
                Ser s = new Ser(nmL[k], td, y, cc[k]); if (k >= 2) { s.Dash = DashStyle.Dash; s.Width = 2.5f; } c1.S.Add(s);
            }
            for (int ci = 0; ci < 2; ci++)
            {
                double[] y = new double[6], e = new double[6], yt = new double[6];
                for (int ti = 0; ti < 6; ti++) { y[ti] = Wm[ci, 0, ti]; e[ti] = We[ci, 0, ti]; yt[ti] = Wm[ci, 1, ti]; }
                Ser s = new Ser((ci == 0 ? "weak" : "strong") + " virtual AFM, β=" + F2(bMean[ci, 0]), T6, y, cc[ci]); s.Line = false; s.Marker = 1; s.E = e; c1.S.Add(s);
                Ser s2 = new Ser((ci == 0 ? "weak" : "strong") + " virtual AFM + tip (R=10 nm), β=" + F2(bMean[ci, 1]), T6, yt, cc[ci]); s2.Line = false; s2.Marker = 3; c1.S.Add(s2);
            }
            c1.Save(Path.Combine(Out, "fig1_W_thickness.png"));

            // Figure 2: PSD
            Chart c2 = new Chart(); c2.Title = "Radial PSD: weak and strong surface diffusion"; c2.XLabel = "Spatial frequency k  [nm⁻¹]"; c2.YLabel = "PSD  [nm⁴]";
            c2.LogX = true; c2.LogY = true; c2.X0 = 0.006; c2.X1 = 1.6; c2.Y0 = 1e-4; c2.Y1 = 10; c2.LegendPos = 2;
            double[] tp = { 10, 44, 120 }; DashStyle[] ds = { DashStyle.Dot, DashStyle.Dash, DashStyle.Solid };
            for (int ci = 0; ci < 2; ci++) for (int k = 0; k < 3; k++)
                {
                    double[][] p = RadialPsdExact(g, base_[ci], tp[k]);
                    Ser s = new Ser((ci == 0 ? "weak" : "strong") + ", h = " + F(tp[k]) + " nm", p[0], p[1], cc[ci]); s.Dash = ds[k]; s.Width = 2.5f; c2.S.Add(s);
                }
            for (int ci = 0; ci < 2; ci++)
            {
                Ser s = new Ser((ci == 0 ? "weak" : "strong") + ", 120 nm, measured with tip", psdTip[ci][0], psdTip[ci][1], ci == 0 ? CO : CP); s.Width = 2f; c2.S.Add(s);
            }
            c2.V.Add(new RefLine(1.0 / TIP_R, "k ≈ 1/R_tip", Color.DimGray));
            c2.Save(Path.Combine(Out, "fig2_PSD.png"));
            List<string> pr = new List<string>();
            for (int ci = 0; ci < 2; ci++) for (int i = 0; i < psdMeas[ci][0].Length; i++) pr.Add(cn[ci] + "," + F(psdMeas[ci][0][i]) + "," + psdMeas[ci][1][i].ToString("E4", IC) + "," + psdTip[ci][1][i].ToString("E4", IC));
            Csv("psd_120nm.csv", "state,k_nm-1,PSD_raw_nm4,PSD_tip_nm4", pr);

            // Figure 5: surface maps (256x256 crop = 512x512 nm)
            SaveMaps(Path.Combine(Out, "fig5_surfaces.png"), maps, new string[] { "Weak diffusion (ℓ_D=1 nm), h=100 nm, W=" + F2(w100v[0, 0]) + " nm", "Strong diffusion (ℓ_D=10 nm), h=100 nm, W=" + F3(w100v[1, 0]) + " nm" }, N, g.dx, 256);
        }

        static void SaveMaps(string path, double[][] us, string[] titles, int N, double dx, int crop)
        {
            int S = 620, gap = 70, top = 90, W = 2 * S + 3 * gap + 150, Hh = S + top + 80;
            using (Bitmap bmp = new Bitmap(W, Hh))
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias; g.TextRenderingHint = TextRenderingHint.AntiAlias; g.Clear(Color.White); g.InterpolationMode = InterpolationMode.NearestNeighbor; g.PixelOffsetMode = PixelOffsetMode.Half;
                Font ft = new Font("Segoe UI", 24f, FontStyle.Bold, GraphicsUnit.Pixel), fk = new Font("Segoe UI", 20f, FontStyle.Regular, GraphicsUnit.Pixel);
                for (int q = 0; q < 2; q++)
                {
                    double[] u = us[q]; double mu = 0; foreach (double v in u) mu += v; mu /= u.Length;
                    double sd = Rms(u), lo = -3 * sd, hi = 3 * sd;
                    using (Bitmap im = new Bitmap(crop, crop))
                    {
                        for (int i = 0; i < crop; i++) for (int j = 0; j < crop; j++) im.SetPixel(j, i, Hot((u[i * N + j] - mu - lo) / (hi - lo)));
                        float x0 = gap + q * (S + gap);
                        g.DrawImage(im, new RectangleF(x0, top, S, S));
                        g.DrawString(titles[q], ft, Brushes.Black, x0, 30);
                        float bar = (float)(100.0 / (crop * dx) * S);
                        g.FillRectangle(Brushes.White, x0 + 20, top + S - 40, bar, 10); g.DrawString("100 nm", fk, Brushes.White, x0 + 20, top + S - 70);
                        g.DrawString("color scale: ±3W = ±" + F2(3 * sd) + " nm", fk, Brushes.Black, x0, top + S + 12);
                    }
                }
                float cx = 3 * gap + 2 * S - 20;
                for (int y = 0; y < S; y++) { Color c = Hot(1 - y / (double)S); g.DrawLine(new Pen(c), cx, top + y, cx + 30, top + y); }
                g.DrawRectangle(Pens.Black, cx, top, 30, S);
                g.DrawString("+3W", fk, Brushes.Black, cx + 36, top - 8); g.DrawString("0", fk, Brushes.Black, cx + 36, top + S / 2 - 12); g.DrawString("−3W", fk, Brushes.Black, cx + 36, top + S - 16);
                bmp.Save(path, ImageFormat.Png);
            }
        }
        static Color Hot(double t)
        {
            t = Math.Max(0, Math.Min(1, t));
            int r = (int)(255 * Math.Min(1, t / 0.4)), gg = (int)(255 * Math.Max(0, Math.Min(1, (t - 0.3) / 0.45))), b = (int)(255 * Math.Max(0, Math.Min(1, (t - 0.7) / 0.3)));
            return Color.FromArgb(r, gg, b);
        }

        // ---------- 3) nonlinear (KPZ term) 2D simulation -- ETD2RK (Cox-Matthews) + exact OU noise
        static void NLTerm(Grid g, Fft2 f, double c, double[] Ur, double[] Ui, double[] outR, double[] outI, double[] br, double[] bi)
        {
            int n = g.N * g.N;
            for (int i = 0; i < n; i++)
            {
                if (!g.keep[i]) { outR[i] = outI[i] = br[i] = bi[i] = 0; continue; }
                outR[i] = -g.kx[i] * Ui[i]; outI[i] = g.kx[i] * Ur[i]; br[i] = -g.ky[i] * Ui[i]; bi[i] = g.ky[i] * Ur[i];
            }
            f.Run(outR, outI, true); f.Run(br, bi, true);
            for (int i = 0; i < n; i++) { double s2 = outR[i] * outR[i] + br[i] * br[i]; outR[i] = 0.5 * c * (s2 < 1e-6 ? s2 : 1 - Math.Exp(-s2)); outI[i] = 0; }  // Dasgupta-Das Sarma-Kim (1996) regularization
            f.Run(outR, outI, false);
            for (int i = 0; i < n; i++) if (!g.keep[i]) { outR[i] = 0; outI[i] = 0; }
            outR[0] = 0; outI[0] = 0;
        }
        static double[] RunNL(Grid g, Fft2 f, Model m, double dt, ulong seed, double[] rec, out bool blew)
        {
            int N = g.N, n = N * N; Rng r = new Rng(seed); blew = false;
            double[] Ur = new double[n], Ui = new double[n], E = new double[n], P = new double[n], P2 = new double[n], Q = new double[n];
            for (int i = 0; i < n; i++)
            {
                double w = m.a * g.L[i] + m.b * g.L[i] * g.L[i], z = w * dt;
                E[i] = Math.Exp(-z);
                if (z < 1e-5) { P[i] = dt * (1 - z / 2); P2[i] = dt * (0.5 - z / 6); }
                else { P[i] = (1 - E[i]) / w; P2[i] = (E[i] - 1 + z) / (w * w * dt); }
                double gg = 2 * z < 1e-8 ? dt : (1 - Math.Exp(-2 * z)) / (2 * w); Q[i] = Math.Sqrt(OMEGA / (g.dx * g.dx) * gg);
            }
            double[] n1r = new double[n], n1i = new double[n], n2r = new double[n], n2i = new double[n], br = new double[n], bi = new double[n];
            double[] Ar = new double[n], Ai = new double[n], nr = new double[n], ni = new double[n];
            double[] W = new double[rec.Length]; int ri = 0; double tau = 0; int steps = (int)Math.Round(rec[rec.Length - 1] / dt);
            bool nl = m.c != 0;
            for (int st = 1; st <= steps; st++)
            {
                if (nl)
                {
                    NLTerm(g, f, m.c, Ur, Ui, n1r, n1i, br, bi);
                    for (int i = 0; i < n; i++) { Ar[i] = E[i] * Ur[i] + P[i] * n1r[i]; Ai[i] = E[i] * Ui[i] + P[i] * n1i[i]; }
                    NLTerm(g, f, m.c, Ar, Ai, n2r, n2i, br, bi);
                }
                for (int i = 0; i < n; i++) { nr[i] = r.G(); ni[i] = 0; }
                f.Run(nr, ni, false);
                for (int i = 1; i < n; i++)
                {
                    if (nl) { Ur[i] = Ar[i] + P2[i] * (n2r[i] - n1r[i]) + Q[i] * nr[i]; Ui[i] = Ai[i] + P2[i] * (n2i[i] - n1i[i]) + Q[i] * ni[i]; }
                    else { Ur[i] = E[i] * Ur[i] + Q[i] * nr[i]; Ui[i] = E[i] * Ui[i] + Q[i] * ni[i]; }
                }
                Ur[0] = 0; Ui[0] = 0; tau = st * dt;
                if (st % 50 == 0 || (ri < rec.Length && tau >= rec[ri] - 1e-9))
                {
                    double s = 0; for (int i = 1; i < n; i++) s += Ur[i] * Ur[i] + Ui[i] * Ui[i];
                    double w = Math.Sqrt(s) / n;
                    if (double.IsNaN(w) || w > 50) { blew = true; Log("   blow-up: tau=" + F(tau) + " nm"); return W; }
                    while (ri < rec.Length && tau >= rec[ri] - 1e-9) { W[ri] = w; ri++; }
                }
            }
            return W;
        }
        static void Nonlinear()
        {
            Grid g = new Grid(256, 2.0); Fft2 f = new Fft2(256);
            double[] rec = { 5, 7, 10, 13, 16, 20, 27, 35, 44, 55, 60, 73, 90, 100, 120 };
            int[] i6 = { 2, 4, 6, 8, 11, 14 };
            Log("=== NONLINEAR 2D SIMULATION (N=256, dx=2 nm, L=512 nm, ETD2RK + 2/3 de-aliasing + Dasgupta regularization (1-exp(-|grad h|^2))) ===");
            double[] lD = { 1, 10 }; double[] cs = { 0, 1, -1, 3 }; int R = 3;
            Chart c6 = new Chart(); c6.Title = "Effect of the nonlinear (λ) term"; c6.XLabel = "Film thickness h  [nm]"; c6.YLabel = "RMS roughness W  [nm]";
            c6.LogX = true; c6.LogY = true; c6.X0 = 5; c6.X1 = 150; c6.Y0 = 0.01; c6.Y1 = 2; c6.LegendPos = 3;
            Color[] cc = { CW, CS }; DashStyle[] dsc = { DashStyle.Solid, DashStyle.Dash, DashStyle.DashDot, DashStyle.Dot };
            List<string> rows = new List<string>();
            for (int ci = 0; ci < 2; ci++)
            {
                Model lin = new Model("", 0, Model.BFromLD(lD[ci]), 0, 0, 0);
                double[] wex = new double[rec.Length]; for (int i = 0; i < rec.Length; i++) wex[i] = WExact(g, lin, rec[i]);
                for (int k = 0; k < cs.Length; k++)
                {
                    // with weak diffusion (l_D < pixel) the lambda != 0 case is numerically unstable; only c = 0 is solved
                    if (ci == 0 && cs[k] != 0) continue;
                    Model m = new Model("", 0, Model.BFromLD(lD[ci]), cs[k], 0, 0);
                    double[] wm = new double[rec.Length]; double[] bs = new double[R]; bool anyBlow = false;
                    for (int rr = 0; rr < R; rr++)
                    {
                        bool blew; double[] w = RunNL(g, f, m, 0.1, (ulong)(77 + 100 * ci + 10 * k + rr), rec, out blew);
                        if (blew) { anyBlow = true; break; }
                        for (int i = 0; i < rec.Length; i++) wm[i] += w[i] / R;
                        double[] x = new double[6], y = new double[6]; for (int q = 0; q < 6; q++) { x[q] = Math.Log(rec[i6[q]]); y[q] = Math.Log(w[i6[q]]); }
                        bs[rr] = Fit(x, y, null)[0];
                    }
                    if (anyBlow) { Log((ci == 0 ? "weak" : "strong") + " c=" + F(cs[k]) + ": NUMERICAL BLOW-UP"); continue; }
                    double bm = 0; foreach (double b in bs) bm += b / R; double bsd = 0; foreach (double b in bs) bsd += (b - bm) * (b - bm); bsd = Math.Sqrt(bsd / (R - 1));
                    if (ci == 1) { NLmin = Math.Min(NLmin, bm); NLmax = Math.Max(NLmax, bm); haveNL = true; }
                    Log((ci == 0 ? "weak" : "strong") + " c=lambda/F=" + F(cs[k]) + ": beta=" + F3(bm) + " +- " + F3(bsd) + "  W60=" + F3(wm[10]) + "  W100=" + F3(wm[13]) + "  (exact linear W100=" + F3(wex[13]) + ")");
                    for (int i = 0; i < rec.Length; i++) rows.Add((ci == 0 ? "weak" : "strong") + "," + F(cs[k]) + "," + F(rec[i]) + "," + F(wm[i]));
                    Ser s = new Ser((ci == 0 ? "weak" : "strong") + ", λ/F = " + F(cs[k]) + ", β=" + F2(bm), rec, wm, cc[ci]); s.Dash = dsc[k]; s.Width = 2.5f; if (k == 0) { s.Marker = 1; s.Line = false; } c6.S.Add(s);
                }
                Ser se = new Ser((ci == 0 ? "weak" : "strong") + ", exact linear solution", rec, wex, Color.Gray); se.Width = 1.5f; c6.S.Add(se);
            }
            // time-step check
            bool b1, b2; Model mt = new Model("", 0, Model.BFromLD(10), 3, 0, 0);
            double[] wa = RunNL(g, f, mt, 0.1, 5, rec, out b1), wb = RunNL(g, f, mt, 0.05, 5, rec, out b2);
            Log("Time-step check (strong, c=3, same seed): W100 dt=0.1 -> " + F3(wa[13]) + ", dt=0.05 -> " + F3(wb[13]));
            Csv("nonlinear_W.csv", "state,c,thickness_nm,W_nm", rows);
            c6.Save(Path.Combine(Out, "fig6_nonlinear.png"));
        }
    }
}
