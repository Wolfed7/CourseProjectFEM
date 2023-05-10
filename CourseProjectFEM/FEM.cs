using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProjectFEM;

public class FEM
{
   private SparseMatrix _globalMatrix;
   private Vector _globalVector;
   private Vector _solution;

   public Mesh Mesh;
   private Solver _solver;
   private Vector _localVector;
   private Matrix _localStiffness;
   private Matrix _localMassSigma, _localMassHee;

   private double _timeLayersNorm;
   private int _currentTimeLayerI => 2;
   private int _prevTimeLayerI => 1;
   private int _prevPrevTimeLayerI => 0;

   private double _currentTimeLayer;
   private double _prevTimeLayer;
   private double _prevPrevTimeLayer;
   private double _t, _t0, _t1;
   private Vector[] _qLayers;

   public static int NodesPerElement => 4;

   public FEM()
   {
      _qLayers = new Vector[3];
      _localStiffness = new(NodesPerElement);
      _localMassSigma = new(NodesPerElement);
      _localMassHee = new(NodesPerElement);
      _localVector = new(NodesPerElement);

      _solver = new BCG();
      Mesh = new Mesh();

      _globalMatrix = new SparseMatrix(0, 0);
      _globalVector = new(0);
      _solution = new(0);
   }

   public void SetSolver(Solver solver)
      => _solver = solver;
   public void SetMesh(Mesh mesh)
   => Mesh = mesh;


   public void Compute()
   {
      BuildPortrait();

      // Костыль - сборка двух слоёв вместо одного начального.
      for (int i = 0; i < _qLayers[0].Size; i++)
         _qLayers[0][i] = Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, Mesh.TimeLayers[0]);

      for (int i = 0; i < _qLayers[1].Size; i++)
         _qLayers[1][i] = Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, Mesh.TimeLayers[1]);

      // DEBUG: решение, которое должно быть найдено в ходе поиска третьего слоя.
      for (int i = 0; i < _qLayers[2].Size; i++)
         Console.WriteLine($"{i}:  {Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, Mesh.TimeLayers[2])}");

      for (int i = 2; i < Mesh.TimeLayers.Count; i++)
      {

         _currentTimeLayer = Mesh.TimeLayers[i];
         _prevTimeLayer = Mesh.TimeLayers[i - 1];
         _prevPrevTimeLayer = Mesh.TimeLayers[i - 2];

         _t = _currentTimeLayer - _prevPrevTimeLayer;
         _t0 = _currentTimeLayer - _prevTimeLayer;
         _t1 = _prevTimeLayer - _prevPrevTimeLayer;

         AssemblySLAE();
         AccountSecondConditions();
         AccountFirstConditions();
         ExcludeFictiveNodes();

         _solver.SetSLAE(_globalVector, _globalMatrix);
         _qLayers[2] = _solver.Solve();

         // DEBUG: погрешность на третьем слое
         //if (i == 2)
         //{
         //   for (int j = 0; j < _qLayers[2].Size; j++)
         //      Console.WriteLine($"{j}:  {Math.Abs(Parameters.U(_mesh.Points[j].X, _mesh.Points[j].Y, _mesh.TimeLayers[2]) - _qLayers[2][j]):e2}");
         //}
         _qLayers[0] = _qLayers[1];
         _qLayers[1] = _qLayers[2];

         PrintCurrentLayerErrorNorm();
      }
      //PrintTimeLayersNorm();
   }

   public void BuildPortrait()
   {
      var list = new HashSet<int>[Mesh.NodesCount].Select(_ => new HashSet<int>()).ToList();
      foreach (var element in Mesh.Elements)
         foreach (var position in element)
            foreach (var node in element)
               if (position > node)
                  list[position].Add(node);

      int offDiagonalElementsCount = list.Sum(childList => childList.Count);

      _globalMatrix = new(Mesh.NodesCount, offDiagonalElementsCount);
      _globalVector = new(Mesh.NodesCount);

      _globalMatrix._ia[0] = 0;

      for (int i = 0; i < list.Count; i++)
         _globalMatrix._ia[i + 1] = _globalMatrix._ia[i] + list[i].Count;

      int k = 0;
      foreach (var childList in list)
         foreach (var value in childList.Order())
            _globalMatrix._ja[k++] = value;

      for (int i = 0; i < _qLayers.Length; i++)
         _qLayers[i] = new Vector(Mesh.NodesCount);
   }

   private void AssemblySLAE()
   {
      _globalVector.Fill(0);
      _globalMatrix.Clear();

      for (int ielem = 0; ielem < Mesh.ElementsCount; ielem++)
      {
         AssemblyLocalSLAE(ielem);
         AddLocalMatrixToGlobal(ielem);
         AddLocalVectorToGlobal(ielem);

         _localStiffness.Clear();
         _localMassSigma.Clear();
         _localMassHee.Clear();
         _localVector.Clear();
      }

      Array.Copy(_globalMatrix._al, _globalMatrix._au, _globalMatrix._al.Length);
   }

   private void AssemblyLocalSLAE(int ielem)
   {
      double hx = Math.Abs(Mesh.Points[Mesh.Elements[ielem][3]].X - Mesh.Points[Mesh.Elements[ielem][0]].X);
      double hy = Math.Abs(Mesh.Points[Mesh.Elements[ielem][3]].Y - Mesh.Points[Mesh.Elements[ielem][0]].Y);

      double coeffG1 = hy / hx / 6;
      double[,] matrixG1 =
{
         { 2.0, -2.0, 1.0, -1.0 },
         { -2.0, 2.0, -1.0, 1.0 },
         { 1.0, -1.0, 2.0, -2.0 },
         { -1.0, 1.0, -2.0, 2.0 }
      };

      double coeffG2 = hx / hy / 6;
      double[,] matrixG2 =
{
         { 2.0, 1.0, -2.0, -1.0 },
         { 1.0, 2.0, -1.0, -2.0 },
         { -2.0, -1.0, 2.0, 1.0 },
         { -1.0, -2.0, 1.0, 2.0 }
      };

      double coeffM = hx * hy / 36;
      double[,] matrixM =
      {
         { 4.0, 2.0, 2.0, 1.0 },
         { 2.0, 4.0, 1.0, 2.0 },
         { 2.0, 1.0, 4.0, 2.0 },
         { 1.0, 2.0, 2.0, 4.0 }
      };

      var _tempMass = new Matrix(NodesPerElement);
      for (int i = 0; i < NodesPerElement; i++)
      {
         for (int j = 0; j < NodesPerElement; j++)
         {
            _localMassSigma[i, j] = Parameters.Sigma() * coeffM * matrixM[i, j];
            _localMassHee[i, j] = Parameters.Hee() * coeffM * matrixM[i, j];
            // В матрицу жёсткости запишу всю локальную А.
            _localStiffness[i, j] = Parameters.Lambda() * coeffG1 * matrixG1[i, j] + Parameters.Lambda() * coeffG2 * matrixG2[i, j]
               + 2 / _t / _t0 * _localMassHee[i, j] + (_t + _t0) / _t / _t0 * _localMassSigma[i, j];


            _tempMass[i, j] = matrixM[i, j];
         }
      }

      for (int i = 0; i < NodesPerElement; i++)
         _localVector[i] = Parameters.F(Mesh.Points[Mesh.Elements[ielem][i]].X, Mesh.Points[Mesh.Elements[ielem][i]].Y, _currentTimeLayer);

      Vector qLocalPrevPrev = new(NodesPerElement);
      Vector qLocalPrev = new(NodesPerElement);

      for (int i = 0; i < NodesPerElement; i++)
      {
         qLocalPrevPrev[i] = _qLayers[_prevPrevTimeLayerI][Mesh.Elements[ielem][i]];
         qLocalPrev[i] = _qLayers[_prevTimeLayerI][Mesh.Elements[ielem][i]];
      }

      // Вектор правой части d (тоже локальный)
      _localVector = coeffM * _tempMass * _localVector
         - 2 / _t1 / _t * _localMassHee * qLocalPrevPrev
         + 2 / _t1 / _t0 * _localMassHee * qLocalPrev
         - (_t0 / _t / _t1) * _localMassSigma * qLocalPrevPrev
         + _t / _t1 / _t0 * _localMassSigma * qLocalPrev;
   }

   private void AddLocalMatrixToGlobal(int ielem)
   {
      for (int i = 0; i < NodesPerElement; i++)
      {
         for (int j = 0; j < NodesPerElement; j++)
         {
            if (Mesh.Elements[ielem][i] == Mesh.Elements[ielem][j])
            {
               _globalMatrix._di[Mesh.Elements[ielem][i]] += _localStiffness[i, j];
               continue;
            }

            if (Mesh.Elements[ielem][i] > Mesh.Elements[ielem][j])
            {
               for (int icol = _globalMatrix._ia[Mesh.Elements[ielem][i]]; icol < _globalMatrix._ia[Mesh.Elements[ielem][i] + 1]; icol++)
               {
                  if (_globalMatrix._ja[icol] == Mesh.Elements[ielem][j])
                  {
                     _globalMatrix._al[icol] += _localStiffness[i, j];
                     break;
                  }
               }
            }
         }
      }
   }

   private void AddLocalVectorToGlobal(int ielem)
   {
      for (int i = 0; i < NodesPerElement; i++)
         _globalVector[Mesh.Elements[ielem][i]] += _localVector[i];
   }

   public void AccountSecondConditions()
   {
      for (int i = 0; i < Mesh.BoundaryRibs2.Count; i++)
      {
         for (int j = 0; j < Mesh.BoundaryRibs2[i].Count; j++)
         {
            double h = Math.Max
               (
               Math.Abs(Mesh.Points[i].X - Mesh.Points[Mesh.BoundaryRibs2[i][j]].X),
               Math.Abs(Mesh.Points[i].Y - Mesh.Points[Mesh.BoundaryRibs2[i][j]].Y)
               );

            double Theta1 = Parameters.dU_dn(Mesh.Points[i].X, Mesh.Points[i].Y, _currentTimeLayer);
            double Theta2 = Parameters.dU_dn(Mesh.Points[Mesh.BoundaryRibs2[i][j]].X, Mesh.Points[Mesh.BoundaryRibs2[i][j]].Y, _currentTimeLayer);

            _globalVector[i] += h / 6.0 * (2.0 * Theta1 + Theta2);
            _globalVector[Mesh.BoundaryRibs2[i][j]] += h / 6.0 * (Theta1 + 2.0 * Theta2);
         }
      }
   }

   public void AccountFirstConditions()
   {
      foreach (var node in Mesh.BoundaryNodes1)
      {
         int row = node;

         // На диагонали единица.
         _globalMatrix._di[row] = 1;
         
         // В векторе правой части значение краевого.
         _globalVector[row] = Parameters.U(Mesh.Points[node].X, Mesh.Points[node].Y, _currentTimeLayer);

         // Вся остальная строка 0. 
         for (int i = _globalMatrix._ia[row]; i < _globalMatrix._ia[row + 1]; i++)
            _globalMatrix._al[i] = 0;

         for (int col = row + 1; col < _globalMatrix.Size; col++)
         {
            for (int j = _globalMatrix._ia[col]; j < _globalMatrix._ia[col + 1]; j++)
            {
               if (_globalMatrix._ja[j] == row)
               {
                  _globalMatrix._au[j] = 0;
                  break;
               }
            }
         }
      }
   }

   public void ExcludeFictiveNodes()
   {
      foreach (var node in Mesh.FictiveNodes)
      {
         int row = node;

         // На диагонали единица.
         _globalMatrix._di[row] = 1;

         // В векторе правой части 0.
         _globalVector[row] = 0;

         // Вся остальная строка 0. 
         for (int i = _globalMatrix._ia[row]; i < _globalMatrix._ia[row + 1]; i++)
            _globalMatrix._al[i] = 0;

         for (int col = row + 1; col < _globalMatrix.Size; col++)
         {
            for (int j = _globalMatrix._ia[col]; j < _globalMatrix._ia[col + 1]; j++)
            {
               if (_globalMatrix._ja[j] == row)
               {
                  _globalMatrix._au[j] = 0;
                  break;
               }
            }
         }
      }
   }

   public void PrintCurrentLayerErrorNorm()
   {
      double errorNorm = 0;
      double realSolutionNorm = 0;

      // Без узлов с первыми краевыми условиями.
      for (int i = 0; i < _qLayers[2].Size; i++)
         if (!Mesh.BoundaryNodes1.Contains(i) && !Mesh.FictiveNodes.Contains(i))
         {
            errorNorm += (_qLayers[2][i] - Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, _currentTimeLayer))
               * (_qLayers[2][i] - Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, _currentTimeLayer));
            realSolutionNorm += Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, _currentTimeLayer)
               * Parameters.U(Mesh.Points[i].X, Mesh.Points[i].Y, _currentTimeLayer);
         }

      var res = Math.Sqrt(errorNorm) / Math.Sqrt(realSolutionNorm);
      _timeLayersNorm += res * res;
      Console.WriteLine($"{res:e2}");
   }

   public void PrintTimeLayersNorm()
   {
      Console.WriteLine($"{Math.Sqrt(_timeLayersNorm):e2}");
   }
}
