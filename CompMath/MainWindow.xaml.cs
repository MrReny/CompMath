// Программа по выч мату.

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.Remoting.Messaging;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using LiveCharts;
using LiveCharts.Configurations;
using LiveCharts.Defaults;
using LiveCharts.Wpf;

namespace CompMath
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow
    {
        private double _lX = 0; // lower edge
        private double _uX = 2; // upper edge

        private double _y0 = 15; // y(0) 
        private int _n = 41; // number of steps
        
        private double _h = 0.05; // value of step

        private double[] _xMass { get; set; }
    
        private double[] _yMass { get; set; }

        private double[] _bvxMass;
        private double[] _bvyMass;
        
        private double[] _r;
        
        public Func<double, string> XFormatter { get; set; }
        
        public SeriesCollection SeriesCollection { get; set; }

        public double H // property for _h
        {
            set
            {
                if (value > 0.1)
                {
                    H = _h;
                }
                else
                {
                    _h = value;
                }

                _n = (int)((_uX - _lX) / _h)+1;
            }
        }

        public List<XYtable> XYlist;
        public List<XYdYtable> XYdYlist;
        
        public MainWindow()
        {
            InitializeComponent();
        }

        private double _Function(double x)
        {
            double fx = 25 * Math.Cos(7.5 * x) / (x + 2);
            return fx;
        }

        // Euler
        private void EulerMethod()
        {
            _xMass = new double[_n];
            _yMass = new double[_n];

            _xMass[0] = 0;
            _yMass[0] = _y0;
            
            XYlist = new List<XYtable>(_n);
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));
            
            for (int i = 1; i < _n; i++)
            {
                _yMass[i] = _yMass[i - 1] + _h * _Function(_xMass[i - 1]);
                _xMass[i] = _xMass[i - 1] + _h;
                
                XYlist.Add(new XYtable(_xMass[i],_yMass[i],i));
            }

            Table.ItemsSource = XYlist;
            
            DrawGraph("Euler");
        }
        
        // dont even work, because of differentiating on one argument
        private void EulerMethod_mod2()
        {
            _xMass = new double[_n+1];
            _yMass = new double[_n+1];

            _xMass[0] = 0;
            _yMass[0] = _y0;
            
            XYlist = new List<XYtable>(_n);
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));

            double dY = 0;
            
            for (int i = 1; i < _n+1; i++)
            {
                dY = _h * _Function(_xMass[i - 1] * _h / 2);
                _yMass[i] = _yMass[i - 1] + dY;
                _xMass[i] = _xMass[i - 1] + _h;
                
                XYlist.Add(new XYtable(_xMass[i],_yMass[i],i));
            }

            Table.ItemsSource = XYlist;
            
            DrawGraph("Euler 2 - dont work");
        }
        
        private void RungeKuttMethod()
        {
            _xMass = new double[_n];
            _yMass = new double[_n];

            _xMass[0] = 0;
            _yMass[0] = _y0;
            
            XYlist = new List<XYtable>(_n);
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));

            double dY = 0;

            double K1 = 0;
            double K2 = 0;
            double K3 = 0;
            double K4 = 0;
            
            for (int i = 1; i < _n; i++)
            {
                K1 = _Function(_xMass[i - 1]);
                K2 = _Function(_xMass[i - 1] + _h / 2);
                K3 = _Function(_xMass[i - 1] + _h / 2);
                K4 = _Function(_xMass[i - 1] + _h);
                
                dY = _h / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
                _yMass[i] = _yMass[i - 1] + dY;
                _xMass[i] = _xMass[i - 1] + _h;
                
                XYlist.Add(new XYtable(_xMass[i],_yMass[i],i));
            }

            Table.ItemsSource = XYlist;
            
            DrawGraph("Runge\n-Kutt");
        }

        private void DrawGraph(string s)
        {
            var obsPlist = new List<ObservablePoint>();
            for(int i =0; i<_n; i++)
            {
                obsPlist.Add(new ObservablePoint(_xMass[i],_yMass[i]));
            }
            SeriesCollection = new SeriesCollection()
            {
                new LineSeries()
                {
                    Title = s,
                    Values = new ChartValues<ObservablePoint>(obsPlist),
                    Fill=Brushes.Transparent
                }
            };

            CartesianChart1.Series = SeriesCollection;
            
            XFormatter = value => value.ToString("C");
            
            DataContext = this; // very important
        }
        
        private double[][] NewtonPolynomial(double[][] dy=null, int n=0)
        {
            double[][] dY = new double[_n][];
            int i = 0;
            
            if (dy == null)
            {
                dy = new double[_n][];
                i = 0;
                foreach (var y in _yMass)
                {
                    dy[i] = new[] {y};
                    i++;
                }
                for(i = 0;i<_n;i++)
                {
                    dY[i] = new double[n+1];
                }
                for (i = 0; i < _n-1 ; i++)
                {
                    dY[i][n] = Math.Abs(dy[i + 1][n] - dy[i][n]);
                }
            }
            else
            {
                for(i = 0;i<_n;i++)
                {
                    dY[i] = new double[n+1];
                    for (var j = 0; j < n; j++)
                    {
                        dY[i][j] = dy[i][j];
                    }
                }
                for (i = 0; i < _n-1 ; i++)
                {
                    dY[i][n] = Math.Abs(dy[i+1][n-1] - dy[i][n-1]);
                    if (0.000000001 >= dY[i+1][n - 1] && dY[i+1][n - 1] >= -0.000000001) dY[i][n] = 0;
                }
            }

            if (n == 7) return dY;
            if (Math.Abs(dY[8][n] - dY[9][n]) > 0.001) dY = NewtonPolynomial(dY, n + 1);

            return dY;
        }

        private double[] GaussMethod(double[][] a, double[] b)
        {
            double dAd;
            double[] r = new double[b.Length];
            for (int k = 0; k < a.Length; k++)
            {
                for (int i = 1+k; i < a.Length; i++)
                {
                    dAd = -a[i][k]/a[k][k];
                    
                    for (int j = 0; j < a[0].Length; j++)
                    {
                        a[i][j] = a[i][j] + a[k][j] * dAd;
                        
                    }
                    b[i] += b[k]*dAd; 
                }
            }

            for (int i = b.Length-1; i >=0 ; i--)
            {
                double lS = b[i];
                for (int j = b.Length-1; j > i; j--)
                {
                    lS -= a[i][j] * r[j];
                }
                r[i] = lS / a[i][i];
            }

            return r;
        }

        private double[] QuadsMethod(int d)
        {

            double[][] sX = new double[d+1][];

            double[] sY = new double[d+1];
            

            for (int j = 0; j < d + 1; j++)
            {
                sX[j] = new double[d+1];
                for (int i = j; i < j+d+1; i++)
                {
                    foreach (var x in _xMass)
                    {
                        
                        sX[j][i-j] += Math.Pow(x, i);
                    }
                }

                for (int i = 0; i < _yMass.Length; i++)
                {
                    sY[j] += _yMass[i] * Math.Pow(_xMass[i], j);
                }
            }
            
            return GaussMethod(sX,sY);
        }
        
        private double QuadsFunc(double[] r,double x)
        {
            double q=0;

            for (int i = r.Length-1; i >=0 ; i--)
            {
                q += r[i] * Math.Pow(x, i);
            }
            return q; 
        }

        private void RecalculateWithQuadsMethod(double[] r)
        {
            _xMass = new double[_n];
            _yMass = new double[_n];
            
            XYlist = new List<XYtable>(_n);

            _xMass[0] = 0;
            _yMass[0] = QuadsFunc(r,0);
            
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));
            
            for (int i = 1; i < _n; i++)
            {
                _xMass[i] = _xMass[i - 1] + _h;
                _yMass[i] = QuadsFunc(r,_xMass[i]);
                XYlist.Add(new XYtable(_xMass[i],_yMass[i],i));
            }

        }
        
        private string _FuncToString(double[] r, object sender=null)
        {
            char[] cl = new char[10]{'⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'};
            string s = "";
            if (sender != null)
            {
                var sr = (Button) sender;
                s=sr.Content.ToString().Trim('0'); 
            }

            else s = "P(x) = ";
            for (int i = r.Length-1; i >= 0; i--)
            {
                if (i < r.Length - 1 && r[i] > 0) s += "+";
                s += Math.Round(r[i],4).ToString() + " x" + cl[i] +' ';
                if (i == r.Length / 2) s += "\n";
            }
            return s;
        }

        private double[] BirgeVietaMethod(double[] r)
        {
            double[] b = new double[r.Length];
            double[] c = new double[r.Length];
            
            List<double[]> xx = new List<double[]>(); 

            for (int i = 1; i < _n; i++)
            {
                if (_bvyMass[i - 1] > 0 && _bvyMass[i] < 0 || _bvyMass[i - 1] < 0 && _bvyMass[i] > 0)
                {
                    xx.Add(new double[]{_bvxMass[i-1], _bvxMass[i]});
                }
            }

            double[] x = new double[xx.Count];
            
            for(int i = 0; i<xx.Count;i++)
            {
                x[i] = (xx[i][0] + xx[i][1]) / 2;
            }

            for(int k=0;k<x.Length;k++)
            {
                
                for (int j = 0; j < r.Length; j++)
                {
                    
                
                b[r.Length-1] = r[r.Length-1];
                c[r.Length-1] = r[r.Length-1];
                for (int i = r.Length - 2; i >= 0; i--)
                {
                    b[i] = r[i] + x[k] * b[i + 1];
                    c[i] = b[i] + x[k] * c[i + 1];
                }

                x[k] = x[k] - (b[0] / c[1]);
                }
            }

            foreach (var xn in x)
            {
                if (0.0002 > QuadsFunc(r, xn) && QuadsFunc(r, xn) > -0.0001)
                {
                }
                else return null;
            }
            
            return  x;
        }

        private void Table_OnLoaded(object sender, RoutedEventArgs e)
        {
            if(SeriesCollection!=null) return;
            EulerMethod();
        }

        private void ButtonEuler_OnClick(object sender, RoutedEventArgs e)
        {
            EulerMethod();
        }

        private void ButtonRungeKutt_OnClick(object sender, RoutedEventArgs e)
        {
            RungeKuttMethod();
        }

        private void ButtonNewtonPolynomial_OnClick(object sender, RoutedEventArgs e)
        {
            double[][] dY;
            XYdYlist = new List<XYdYtable>(_n+1);
            
            dY = NewtonPolynomial();

            for (int i = 0; i < _n;i++)
            {
                XYdYlist.Add(new XYdYtable(_xMass[i],_yMass[i],i,dY[i]));
            }

            Table.ItemsSource = XYdYlist;
        }

        private void ButtonOk_OnClick(object sender, RoutedEventArgs e)
        {
            var h = double.Parse(HBox.Text, System.Globalization.CultureInfo.InvariantCulture);

            H = h;
        }

        private void HBox_OnPreviewTextInput(object sender, TextCompositionEventArgs e)
        {
            if (!(Char.IsDigit(e.Text, 0) || (e.Text == ".")
                  && (!HBox.Text.Contains(".")
                      && HBox.Text.Length != 0)))
            {
                e.Handled = true;
            }
        }

        private void Table_OnAutoGeneratingColumn(object sender, DataGridAutoGeneratingColumnEventArgs e)
        {
            switch (e.Column.Header.ToString())
            {
                case "I":
                    e.Column.Visibility = Visibility.Collapsed;
                    break;
                case "X":
                    e.Column.Visibility = Visibility.Collapsed;
                    break;
                case "Y":
                    e.Column.Visibility = Visibility.Collapsed;
                    break;
            }
        }

        private void Quads_OnClick(object sender, RoutedEventArgs e)
        {
            if (XYdYlist == null) return;
            _r = QuadsMethod(XYdYlist[1].dY.Length);
            
            RecalculateWithQuadsMethod(_r);

            Table1.ItemsSource = XYlist;

            QuadFunc.Content = _FuncToString(_r);
            
            
            var obsPlist = new List<ObservablePoint>();
            for(int i =0; i<_n; i++)
            {
                obsPlist.Add(new ObservablePoint(_xMass[i],_yMass[i]));
            }
            if (SeriesCollection.FirstOrDefault((item)=> item.Title=="P(x)")!=null) 
            {
                SeriesCollection.FirstOrDefault((item)=> item.Title=="P(x)").Values =
                    new ChartValues<ObservablePoint>(obsPlist);
            }
            else
            {
                SeriesCollection.Add(new LineSeries()
                {
                    Title = "P(x)",
                    Values = new ChartValues<ObservablePoint>(obsPlist),
                    Fill = Brushes.Transparent
                });
            }
            
        }

        private void Squares_OnClick(object sender, RoutedEventArgs e)
        {
            if(_r==null) return;
            double[] fr= new double[_r.Length];
            double[] fx = new double[_xMass.Length];
            double[] fy = new double[_yMass.Length];

            _r.CopyTo(fr,0);
            fr[0] -= 13;
            fr[1] -= 1.5;
            
            _xMass.CopyTo(fx,0);
            _yMass.CopyTo(fy,0);
            
            RecalculateWithQuadsMethod(fr);
            
            Table2.ItemsSource = XYlist;
            
            PxFunc.Content = _FuncToString(fr,sender);
            
            var obsPlist = new List<ObservablePoint>();
            for(int i =0; i<_n; i++)
            {
                obsPlist.Add(new ObservablePoint(_xMass[i],_yMass[i]));
            }
            
            _bvxMass = new double[_n];
            _bvyMass = new double[_n];
            
            _xMass.CopyTo(_bvxMass,0);
            _yMass.CopyTo(_bvyMass,0);
            
            _xMass = fx;
            _yMass = fy;

            if (SeriesCollection.FirstOrDefault((item)=> item.Title=="P(x)=f(x)")!=null) 
            {
                SeriesCollection.FirstOrDefault((item)=> item.Title=="P(x)=f(x)").Values =
                    new ChartValues<ObservablePoint>(obsPlist);
                
            }
            else
            {
                SeriesCollection.Add(new LineSeries()
                {
                    Title = "P(x)=f(x)",
                    Values = new ChartValues<ObservablePoint>(obsPlist),
                    Fill = Brushes.Transparent
                });
                SeriesCollection.Add(new LineSeries()
                {
                    Title = "y=0",
                    Values = new ChartValues<ObservablePoint>(new ObservablePoint[]
                    {
                        new ObservablePoint(0,0), new ObservablePoint(2,0),
                    })
                });
            }

            var y0x = BirgeVietaMethod(fr);

            if (y0x != null)
            {
                string xstring = "";

                for (int i = 0; i < y0x.Length; i++)
                {
                    xstring += "x" + i + " = " + Math.Round(y0x[i], 6) + ";  ";
                    xstring += "y" + "(x" + i + ") = " + Math.Round(QuadsFunc(fr, y0x[i]), 13) + ";  ";
                    
                    if (i == (y0x.Length-1) / 2) xstring += "\n";
                }

                Xlabel.Content = xstring;
            }

            else Xlabel.Content = "Something gone wrong";

        }
    }

    public class XYtable
    {
        public XYtable(double x, double y, int i)
        {
            this.I = i;
            this.X = x;
            this.Y = Math.Round(y,4);
        }
        
        public int I { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
    }
    public class XYdYtable: XYtable
    {
        public double[] dY;
        public XYdYtable(double x, double y, int i, double[] dy) : base(x, y, i)
        {
            this.dY  = dy;
            var l = dY.Length > 8 ? 8 : dY.Length;
            switch (l)
            {
                case 8:
                    this.dY7 = Math.Round(dY[7],4);
                    goto case 7;
                case 7:
                    this.dY6 = Math.Round(dY[6],4);
                    goto case 6;
                case 6:
                    this.dY5 = Math.Round(dY[5],4);
                    goto case 5;
                case 5:
                    this.dY4 = Math.Round(dY[4],4);
                    goto case 4;
                case 4:
                    this.dY3 = Math.Round(dY[3],4);
                    goto case 3;
                case 3:
                    this.dY2 = Math.Round(dY[2],4);
                    goto case 2;
                case 2:
                    this.dY1 = Math.Round(dY[1],4);
                    goto case 1;
                case 1:
                    this.dY0 = Math.Round(dY[0],4);
                    break;
            }
        }
        
        public double dY0{ get; set;}
        public double dY1{ get; set;}
        public double dY2{ get; set;}
        public double dY3{ get; set;}
        public double dY4{ get; set;}
        public double dY5{ get; set;}
        public double dY6{ get; set;}
        public double dY7{ get; set;}


    }
}