// Программа по выч мату.

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
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
        private int _n = 40; // number of steps
        
        private double _h = 0.05; // value of step

        private double[] _xMass { get; set; }
    
        private double[] _yMass { get; set; }
        
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

                _n = (int)((_uX - _lX) / _h);
            }
        }

        public List<XYtable> XYlist;
        public List<XYdYtable> XYdYlist;
        
        public MainWindow()
        {
            InitializeComponent();
            //EulerMethod();
            
        }

        private double _Function(double x)
        {
            double Fx = 25 * Math.Cos(7.5 * x) / (x + 2);
            return Fx;
        }

        // Euler
        private void EulerMethod()
        {
            _xMass = new double[_n+1];
            _yMass = new double[_n+1];

            _xMass[0] = 0;
            _yMass[0] = _y0;
            
            XYlist = new List<XYtable>(_n);
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));
            
            for (int i = 1; i < _n+1; i++)
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
            _xMass = new double[_n+1];
            _yMass = new double[_n+1];

            _xMass[0] = 0;
            _yMass[0] = _y0;
            
            XYlist = new List<XYtable>(_n);
            XYlist.Add(new XYtable(_xMass[0],_yMass[0],0));

            double dY = 0;

            double K1 = 0;
            double K2 = 0;
            double K3 = 0;
            double K4 = 0;
            
            for (int i = 1; i < _n+1; i++)
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

            var f = XFormatter;
            
            DataContext = this; // very important
        }
        
        private double[][] NewtonPolynomial(double[][] dy=null, int n=0)
        {
            double[][] dY = new double[_n][];
            int i = 0;
            
            if (dy == null)
            {
                dy = new double[_n+1][];
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
                for (i = 0; i < _n ; i++)
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
                }
            }

            if (n == 5) return dY;
            if (Math.Abs(dY[8][n] - dY[9][n]) > 0.05) dY = NewtonPolynomial(dY, n + 1);

            return dY;
        }

        private void GaussMethod(double[] a, double b)
        {
            
        }

        private void QuadsMethod(int d)
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
            
            // TODO: GaussMethod for matrix

            return;
        }

        private void Table_OnLoaded(object sender, RoutedEventArgs e)
        {
            EulerMethod();
        }

        private void ButtonEuler_OnClick(object sender, RoutedEventArgs e)
        {
            EulerMethod();
        }

        private void ButtonRungeKutt_OnClick(object sender, RoutedEventArgs e)
        {
            //CartesianChart1.ClearValue();
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
            QuadsMethod(XYdYlist[1].dY.Length);
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
            switch (dY.Length)
            {
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


    }
}