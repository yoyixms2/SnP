using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
//using System.Numerics;

namespace SnPReader
{
    class S2PData
    {
        public enum FrequencyUnits
        {
            GHz = 0,
            MHz = 1,
            KHz = 2,
            Hz = 3,
        };
        public enum DataFormats
        {
            MagnitudeAngle = 0,
            RealImaginary = 1,
            DecibelAngle = 2
        };
        public enum ParameterTypes
        {
            S = 0,
            Y = 1,
            Z = 2,
            H = 3,
            G = 4
        };

        public int Impedance;
        public readonly FrequencyUnits FileFrequencyUnit;
        public readonly DataFormats FileDataFormat;
        public readonly ParameterTypes ParameterType;
        public readonly Vector<float> Frequency;
        public readonly Vector<Complex32> S11;
        public readonly Vector<Complex32> S21;
        public readonly Vector<Complex32> S12;
        public readonly Vector<Complex32> S22;
        public readonly int PointNumber;
        public S2PData(string filePath)
        {
            string[] FileLines;
            try
            {
                if (File.Exists(filePath))
                {
                    float frequencyFactor =1;
                    FileLines = File.ReadAllLines(filePath);
                    foreach (string line in FileLines)
                    {
                        if (line.Substring(0, 1).ToCharArray()[0] == '#')
                        {
                            var split = line.Split(' ');
                            switch (split[1].ToUpper())
                            {
                                case "GHZ":
                                    FileFrequencyUnit = FrequencyUnits.GHz;
                                    frequencyFactor = (float)Math.Pow(10, 9);
                                    break;
                                case "MHZ":
                                    FileFrequencyUnit = FrequencyUnits.MHz;
                                    frequencyFactor = (float)Math.Pow(10, 6);
                                    break;
                                case "KHZ":
                                    FileFrequencyUnit = FrequencyUnits.KHz;
                                    frequencyFactor = (float)Math.Pow(10, 3);
                                    break;
                                case "HZ":
                                    FileFrequencyUnit = FrequencyUnits.Hz;
                                    frequencyFactor = (float)Math.Pow(10, 0);
                                    break;
                                default:
                                    FileFrequencyUnit = FrequencyUnits.GHz;
                                    frequencyFactor = (float)Math.Pow(10, 9);
                                    break;
                            }
                            switch (split[2])
                            {
                                case "S":
                                    ParameterType = ParameterTypes.S;
                                    break;
                                case "Y":
                                    ParameterType = ParameterTypes.Y;
                                    break;
                                case "Z":
                                    ParameterType = ParameterTypes.Z;
                                    break;
                                case "H":
                                    ParameterType = ParameterTypes.H;
                                    break;
                                case "G":
                                    ParameterType = ParameterTypes.G;
                                    break;
                                default:
                                    ParameterType = ParameterTypes.S;
                                    break;
                            }
                            switch (split[3])
                            {
                                case "RI":
                                    FileDataFormat = DataFormats.RealImaginary;
                                    break;
                                case "MA":
                                    FileDataFormat = DataFormats.MagnitudeAngle;
                                    break;
                                case "DB":
                                    FileDataFormat = DataFormats.DecibelAngle;
                                    break;
                                default:
                                    FileDataFormat = DataFormats.MagnitudeAngle;
                                    break;
                            }
                            if (split[4].Equals("R"))
                            {
                                int.TryParse(split[5], out Impedance);
                            }
                            break;
                        }
                    }
                    List<float> f = new List<float>();
                    List<Complex32> s11 = new List<Complex32>();
                    List<Complex32> s21 = new List<Complex32>();
                    List<Complex32> s12 = new List<Complex32>();
                    List<Complex32> s22 = new List<Complex32>();
                    Vector<double> sFirst = Vector<double>.Build.Dense(4);
                    Vector<double> sSecond = Vector<double>.Build.Dense(4);
                    Vector<double> ten = Vector<double>.Build.Dense(4);
                    ten[0] = 10; ten[1] = 10; ten[2] = 10; ten[3] = 10;
                    foreach (string line in FileLines)
                    {
                        char startChar = line.Substring(0, 1).ToCharArray()[0];
                        if (startChar != '!' && startChar != '#')
                        {
                            string[] split = line.Split(new char[] { ' ' }, System.StringSplitOptions.RemoveEmptyEntries);
                            float.TryParse(split[0], out float parsedFrequency);
                            f.Add(parsedFrequency* frequencyFactor);
                            double.TryParse(split[1], NumberStyles.Any, CultureInfo.InvariantCulture, out double s11First);
                            double.TryParse(split[2], NumberStyles.Any, CultureInfo.InvariantCulture, out double s11Second);
                            double.TryParse(split[3], NumberStyles.Any, CultureInfo.InvariantCulture, out double s21First);
                            double.TryParse(split[4], NumberStyles.Any, CultureInfo.InvariantCulture, out double s21Second);
                            double.TryParse(split[5], NumberStyles.Any, CultureInfo.InvariantCulture, out double s12First);
                            double.TryParse(split[6], NumberStyles.Any, CultureInfo.InvariantCulture, out double s12Second);
                            double.TryParse(split[7], NumberStyles.Any, CultureInfo.InvariantCulture, out double s22First);
                            double.TryParse(split[8], NumberStyles.Any, CultureInfo.InvariantCulture, out double s22Second);
                            sFirst[0] = s11First;
                            sFirst[1] = s21First;
                            sFirst[2] = s12First;
                            sFirst[3] = s22First;
                            sSecond[0] = s11Second;
                            sSecond[1] = s21Second;
                            sSecond[2] = s12Second;
                            sSecond[3] = s22Second;
                            switch (FileDataFormat)
                            {
                                case DataFormats.DecibelAngle:
                                    {
                                        sFirst = ten.PointwisePower(sFirst.Divide(20));
                                        sSecond = sSecond.Multiply(Constants.Pi).Divide(180);
                                        s11.Add(new Complex32((float)(sFirst[0] * Math.Cos(sSecond[0])), (float)(sFirst[0] * Math.Sin(sSecond[0]))));
                                        s21.Add(new Complex32((float)(sFirst[1] * Math.Cos(sSecond[1])), (float)(sFirst[1] * Math.Sin(sSecond[1]))));
                                        s12.Add(new Complex32((float)(sFirst[2] * Math.Cos(sSecond[2])), (float)(sFirst[2] * Math.Sin(sSecond[2]))));
                                        s22.Add(new Complex32((float)(sFirst[3] * Math.Cos(sSecond[3])), (float)(sFirst[3] * Math.Sin(sSecond[3]))));
                                        break;
                                    }
                                case DataFormats.MagnitudeAngle:
                                    {

                                        sSecond = sSecond.Multiply(Constants.Pi).Divide(180);
                                        s11.Add(new Complex32((float)(sFirst[0] * Math.Cos(sSecond[0])), (float)(sFirst[0] * Math.Sin(sSecond[0]))));
                                        s21.Add(new Complex32((float)(sFirst[1] * Math.Cos(sSecond[1])), (float)(sFirst[1] * Math.Sin(sSecond[1]))));
                                        s12.Add(new Complex32((float)(sFirst[2] * Math.Cos(sSecond[2])), (float)(sFirst[2] * Math.Sin(sSecond[2]))));
                                        s22.Add(new Complex32((float)(sFirst[3] * Math.Cos(sSecond[3])), (float)(sFirst[3] * Math.Sin(sSecond[3]))));

                                        break;
                                    }
                                case DataFormats.RealImaginary:
                                    {
                                        s11.Add(new Complex32((float)s11First, (float)s11Second));
                                        s21.Add(new Complex32((float)s21First, (float)s21Second));
                                        s12.Add(new Complex32((float)s12First, (float)s12Second));
                                        s22.Add(new Complex32((float)s22First, (float)s22Second));
                                        break;
                                    }
                            }


                        }

                    }
                    PointNumber = s11.Count;
                    Frequency = Vector<float>.Build.DenseOfArray(f.ToArray());
                    S11 = Vector<Complex32>.Build.DenseOfArray(s11.ToArray());
                    S21 = Vector<Complex32>.Build.DenseOfArray(s21.ToArray());
                    S12 = Vector<Complex32>.Build.DenseOfArray(s12.ToArray());
                    S22 = Vector<Complex32>.Build.DenseOfArray(s22.ToArray());

                }
            }

            catch
            {

            }
        }

}
}
