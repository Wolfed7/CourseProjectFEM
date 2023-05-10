using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProjectFEM;


public class Mesh
{
   private double _timeStart;
   private double _timeEnd;
   private int _timeSplits;
   private double _timeDischarge;
   private List<double> _timeLayers;

   private double[] _linesX;
   private double[] _linesY;
   private List<double> _meshLinesX;
   private List<double> _meshLinesY;
   private int[] _splitsX;
   private int[] _splitsY;
   private double[] _dischargeX;
   private double[] _dischargeY;
   private (int, int, int, int, int)[] _areas;
   private List<int[]> _boundaryConditions;

   private List<List<int>> _allRibs;

   // Кринж.
   private List<Point2D> _points;
   private List<int[]> _elements;
   private IDictionary<int, int> _areaNodes;
   private HashSet<int> _boundaryNodes1;
   private HashSet<int> _boundaryNodes2;
   private List<List<int>> _boundaryRibs2;
   private HashSet<int> _fictiveNodes;

   public int NodesCount => _points.Count;
   public int ElementsCount => _elements.Count;

   public ImmutableList<double> TimeLayers => _timeLayers.ToImmutableList();
   public ImmutableList<Point2D> Points => _points.ToImmutableList();
   public ImmutableList<int[]> Elements => _elements.ToImmutableList();
   public ImmutableDictionary<int, int> AreaNodes => _areaNodes.ToImmutableDictionary();
   public ImmutableHashSet<int> BoundaryNodes1 => _boundaryNodes1.ToImmutableHashSet();
   public ImmutableList<List<int>>  BoundaryRibs2 => _boundaryRibs2.ToImmutableList();
   public ImmutableHashSet<int> FictiveNodes => _fictiveNodes.ToImmutableHashSet();


   public Mesh()
   {
      _allRibs = new();
      _boundaryRibs2 = new();
      _areaNodes = new Dictionary<int, int>();
      _boundaryConditions = new();
      _boundaryNodes1 = new();
      _boundaryNodes2 = new();
      _fictiveNodes = new();

      _timeLayers = new();
      _linesX = Array.Empty<double>();
      _linesY = Array.Empty<double>();
      _meshLinesX = new();
      _meshLinesY = new();
      _points = new();
      _splitsX = Array.Empty<int>();
      _splitsY = Array.Empty<int>();
      _dischargeX = Array.Empty<double>();
      _dischargeY = Array.Empty<double>();
      _areas = Array.Empty<(int, int, int, int, int)>();
      _elements = new();
   }

   public void Input(string filepath1, string filepath2, string filepath3, string filepath4)
   {
      try
      {
         using (var sr = new StreamReader(filepath1))
         {
            sr.ReadLine();
            _linesX = sr.ReadLine().Split().Select(double.Parse).ToArray();

            sr.ReadLine();
            _linesY = sr.ReadLine().Split().Select(double.Parse).ToArray();

            sr.ReadLine();
            _areas = sr.ReadToEnd().Split("\n").Select(row => row.Split())
               .Select(value => 
               (
               int.Parse(value[0]),
               int.Parse(value[1]) - 1,
               int.Parse(value[2]) - 1,
               int.Parse(value[3]) - 1,
               int.Parse(value[4]) - 1
               )
               ).ToArray();
         }

         using (var sr = new StreamReader(filepath2))
         {
            _splitsX = new int[_linesX.Length - 1];
            _dischargeX = new double[_splitsX.Length];

            _splitsY = new int[_linesY.Length - 1];
            _dischargeY = new double[_splitsY.Length];

            var lineX = sr.ReadLine().Split();
            for (int i = 0; i < lineX.Length / 2; i++)
            {
               _splitsX[i] = int.Parse(lineX[2 * i]);
               _dischargeX[i] = double.Parse(lineX[2 * i + 1]);
            }

            var lineY = sr.ReadLine().Split();
            for (int i = 0; i < lineY.Length / 2; i++)
            {
               _splitsY[i] = int.Parse(lineY[2 * i]);
               _dischargeY[i] = double.Parse(lineY[2 * i + 1]);
            }
         }


         using (var sr = new StreamReader(filepath3))
         {
            while(!sr.EndOfStream)
               _boundaryConditions.Add(sr.ReadLine().Split().Select(int.Parse).ToArray());
         }

         using (var sr = new StreamReader(filepath4))
         {
            var data = sr.ReadLine().Split().Select(double.Parse).ToList();
            _timeStart = data[0];
            _timeEnd = data[1];

            _timeSplits = int.Parse(sr.ReadLine());
            _timeDischarge = double.Parse(sr.ReadLine());
         }
      }
      catch (Exception ex)
      {
         Console.WriteLine(ex.Message);
      }
   }

   public void Build()
   {
      // Разбиение каждой области в соответствии с её параметрами:
      // количеством разбиений;
      // коэффициенте разрядки.

      // По времени

      {
         double h;
         double sum = 0;
         double lenght = _timeEnd - _timeStart;

         for (int k = 0; k < _timeSplits; k++)
            sum += Math.Pow(_timeDischarge, k);

         h = lenght / sum;

         _timeLayers.Add(_timeStart);

         while (Math.Round(_timeLayers.Last() + h, 1) < _timeEnd)
         {
            _timeLayers.Add(_timeLayers.Last() + h);
            h *= _timeDischarge;
         }

         _timeLayers.Add(_timeEnd);
      }


      // По оси X
      for (int i = 0; i < _linesX.Length - 1; i++)
      {
         double h;
         double sum = 0;
         double lenght = _linesX[i + 1] - _linesX[i];

         for (int k = 0; k < _splitsX[i]; k++)
            sum += Math.Pow(_dischargeX[i], k);

         h = lenght / sum;

         _meshLinesX.Add(_linesX[i]);

         while (Math.Round(_meshLinesX.Last() + h, 1) < _linesX[i + 1])
         {
            _meshLinesX.Add(_meshLinesX.Last() + h);
            h *= _dischargeX[i];
         }
      }
      _meshLinesX.Add(_linesX.Last());

      // По оси Y
      for (int i = 0; i < _linesY.Length - 1; i++)
      {
         double h;
         double sum = 0;
         double lenght = _linesY[i + 1] - _linesY[i];

         for (int k = 0; k < _splitsY[i]; k++)
            sum += Math.Pow(_dischargeY[i], k);

         h = lenght / sum;

         _meshLinesY.Add(_linesY[i]);

         while (Math.Round(_meshLinesY.Last() + h, 1) < _linesY[i + 1])
         {
            _meshLinesY.Add(_meshLinesY.Last() + h);
            h *= _dischargeY[i];
         }
      }
      _meshLinesY.Add(_linesY.Last());

      // Сборка списка узлов.
      // Узлы нумеруются слева направо, снизу вверх.
      for (int i = 0; i < _meshLinesY.Count; i++)
         for (int j = 0; j < _meshLinesX.Count; j++)
            _points.Add(new(_meshLinesX[j], _meshLinesY[i]));

      // Сборка списка элементов.
      int splitsX = _splitsX.Sum();
      int splitsY = _splitsY.Sum();
      int correction = 0;
      for (int i = 0; i < splitsY; i++)
      {
         for (int j = 0; j < splitsX; j++)
         {
            _elements.Add(
            new int[4]
            {
            i * splitsX + j + correction,
            i * splitsX + j + 1 + correction,
            (i + 1) * (splitsX + 1) + j,
            (i + 1) * (splitsX + 1) + j + 1
            }
            );
         }
         correction++;
      }

      // TODO: Ща будет кринж с определением краевых нодов и фэйковых.
      for (int i = 0; i < _points.Count; i++)
      {
         for (int j = 0; j < _boundaryConditions.Count; j++)
         {
            if (_linesX[_boundaryConditions[j][2] - 1] <= _points[i].X && _points[i].X <= _linesX[_boundaryConditions[j][3] - 1])
            {
               if (_linesY[_boundaryConditions[j][4] - 1] <= _points[i].Y && _points[i].Y <= _linesY[_boundaryConditions[j][5] - 1])
               {
                  if (_boundaryConditions[j][0] == 1)
                     _boundaryNodes1.Add(i);
                  if (_boundaryConditions[j][0] == 2)
                     _boundaryNodes2.Add(i);
                  //break;
               }
            }
         }

         for (int j = 0; j < _areas.Length; j++)
         {
            if (_linesX[_areas[j].Item2] <= _points[i].X && _points[i].X <= _linesX[_areas[j].Item3])
            {
               if (_linesY[_areas[j].Item4] <= _points[i].Y && _points[i].Y <= _linesY[_areas[j].Item5])
               {
                  _areaNodes.Add(i, j);
                  break;
               }
            }
         }

         if(!_areaNodes.ContainsKey(i))
            _fictiveNodes.Add(i);
      }


      // TODO: Кринж со списком рёбер.
      _boundaryRibs2 = new List<int>[_points.Count].Select(_ => new List<int>()).ToList();
      _allRibs = new List<int>[_points.Count].Select(_ => new List<int>()).ToList();
      foreach (var element in Elements)
         foreach (var position in element)
            foreach (var node in element)
               if (position < node)
                  // Нестабильная штука.
                  if (_points[position].X == _points[node].X || _points[position].Y == _points[node].Y)
                  {
                     if(!_fictiveNodes.Contains(position) && !_fictiveNodes.Contains(node))
                        _allRibs[position].Add(node);
                     if (_boundaryNodes2.Contains(position) && _boundaryNodes2.Contains(node) && !_boundaryRibs2[position].Contains(node))
                        _boundaryRibs2[position].Add(node);
                  }

   }

   public void Output(string filepath1, string filepath2)
   {
      try
      {
         using (var sw = new StreamWriter(filepath1))
         {
            foreach (var point in _points)
            sw.WriteLine(point.ToString());
         }

         using (var sw = new StreamWriter(filepath2))
         {
            for (int i = 0; i < _allRibs.Count; i++)
               for (int j = 0; j < _allRibs[i].Count; j++)
                  sw.WriteLine($"{_points[i]} {_points[_allRibs[i][j]]}");
         }
      }
      catch (Exception ex)
      {
         Console.WriteLine(ex.Message);
      }
   }
}