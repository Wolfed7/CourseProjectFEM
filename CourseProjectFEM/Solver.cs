using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace CourseProjectFEM;

public abstract class Solver
{
   protected SparseMatrix _matrix;
   protected Vector _vector;
   protected Vector _solution;

   public double Eps { get; init; }
   public int MaxIters { get; init; }
   public double SolvationTime { get; protected set; }

   public Solver(double eps = 1e-14, int maxIters = 2000)
   {
      Eps = eps;
      MaxIters = maxIters;

      _matrix = new SparseMatrix(0, 0);
      _vector = new Vector(0);
      _solution = new Vector(0);
   }

   public void SetSLAE(Vector vector, SparseMatrix matrix)
   {
      _vector = vector;
      _matrix = matrix;
   }

   public abstract Vector Solve();

   protected void DecomposeLU()
   {
      for (int i = 0; i < _matrix.Size; i++)
      {
         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
         {
            int jColumn = _matrix._ja[j];
            int jk = _matrix._ia[jColumn];
            int k = _matrix._ia[i];

            int shift = _matrix._ja[_matrix._ia[i]] - _matrix._ja[_matrix._ia[jColumn]];

            if (shift > 0)
               jk += shift;
            else
               k -= shift;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < _matrix._ia[jColumn + 1]; k++, jk++)
            {
               sumL += _matrix._al[k] * _matrix._au[jk];
               sumU += _matrix._au[k] * _matrix._al[jk];
            }

            _matrix._al[j] -= sumL;
            _matrix._au[j] -= sumU;
            _matrix._au[j] /= _matrix._di[jColumn];
         }

         double sumD = 0.0;
         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
            sumD += _matrix._al[j] * _matrix._au[j];

         _matrix._di[i] -= sumD;
      }
   }

   protected void ForwardElimination()
   {
      for (int i = 0; i < _matrix.Size; i++)
      {
         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
         {
            _solution[i] -= _matrix._al[j] * _solution[_matrix._ja[j]];
         }

         _solution[i] /= _matrix._di[i];
      }
   }

   protected void BackwardSubstitution()
   {
      for (int i = _matrix.Size - 1; i >= 0; i--)
         for (int j = _matrix._ia[i + 1] - 1; j >= _matrix._ia[i]; j--)
            _solution[_matrix._ja[j]] -= _matrix._au[j] * _solution[i];
   }

   public void PrintSolution(string format = "e14")
   {
      for (int i = 0; i < _solution.Size; i++)
         Console.WriteLine(_solution[i]);
   }
}

public class LOS : Solver
{

   public LOS(double eps = 1e-14, int maxIters = 2000)
   {
      Eps = eps;
      MaxIters = maxIters;
   }

   public override Vector Solve()
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);

      Vector r = _vector - _matrix * _solution;
      Vector z = 1 * r;
      Vector p = _matrix * z;
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double discrepancy = r * r;

      for (int i = 1; i <= MaxIters && discrepancy > Eps; i++)
      {
         alpha = (p * r) / (p * p);
         _solution += alpha * z;

         r -= alpha * p;
         tmp = _matrix * r;

         beta = -(p * tmp) / (p * p);

         z = r + beta * z;
         p = tmp + beta * p;

         discrepancy = r * r;
      }

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
   }
}

public class LU : Solver
{
   public override Vector Solve()
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);
      _matrix = _matrix.ConvertToProfile();

      Stopwatch sw = Stopwatch.StartNew();

      DecomposeLU();
      ForwardElimination();
      BackwardSubstitution();

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
   }
}

public class LOSWithLU : Solver
{
   public override Vector Solve()
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);

      SparseMatrix matrixLU = new(_matrix.Size, _matrix._ja.Length);
      SparseMatrix.Copy(_matrix, matrixLU);

      PartitialLU(matrixLU);

      Vector r = Forward(matrixLU, _vector - _matrix * _solution);
      Vector z = Backward(matrixLU, r);
      Vector p = Forward(matrixLU, _matrix * z);
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double discrepancy = r * r;


      for (int i = 1; i <= MaxIters && discrepancy > Eps; i++)
      {
         alpha = (p * r) / (p * p);
         _solution += alpha * z;

         r -= alpha * p;

         tmp = Forward(matrixLU, _matrix * Backward(matrixLU, r));
         beta = -(p * tmp) / (p * p);

         z = Backward(matrixLU, r) + beta * z;
         p = tmp + beta * p;

         discrepancy = r * r;
      }

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
   }

   protected static void PartitialLU(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            int jCol = Matrix._ja[j];
            int jk = Matrix._ia[jCol];
            int k = Matrix._ia[i];

            int sdvig = Matrix._ja[Matrix._ia[i]] - Matrix._ja[Matrix._ia[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix._ia[jCol + 1]; k++, jk++)
            {
               sumL += Matrix._al[k] * Matrix._au[jk];
               sumU += Matrix._au[k] * Matrix._al[jk];
            }

            Matrix._al[j] -= sumL;
            Matrix._au[j] -= sumU;
            Matrix._au[j] /= Matrix._di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
            sumD += Matrix._al[j] * Matrix._au[j];

         Matrix._di[i] -= sumD;
      }
   }

   protected static Vector Forward(SparseMatrix Matrix, Vector b)
   {
      var result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            result[i] -= Matrix._al[j] * result[Matrix._ja[j]];
         }

         result[i] /= Matrix._di[i];
      }

      return result;
   }

   protected static Vector Backward(SparseMatrix Matrix, Vector b)
   {
      var result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix._ia[i + 1] - 1; j >= Matrix._ia[i]; j--)
         {
            result[Matrix._ja[j]] -= Matrix._au[j] * result[i];
         }
      }

      return result;
   }
}

public class BCG : Solver
{
   public BCG(double eps = 1e-14, int maxIters = 2000)
   {
      Eps = eps;
      MaxIters = maxIters;
   }


   public override Vector Solve()
   {
      _solution = new(_vector.Size);

      Vector residual = _vector - _matrix * _solution;

      Vector p = new(residual.Size);
      Vector z = new(residual.Size);
      Vector s = new(residual.Size);

      Vector.Copy(residual, p);
      Vector.Copy(residual, z);
      Vector.Copy(residual, s);


      Stopwatch sw = Stopwatch.StartNew();

      double vecNorm = _vector.Norm();
      double discrepancy = 1;
      double prPrev = p * residual;

      for (int i = 1; i <= MaxIters && discrepancy > Eps; i++)
      {
         var Az = _matrix * z;
         double alpha = prPrev / (s * Az);

         _solution = _solution + alpha * z;
         residual = residual - alpha * Az;
         p = p - alpha * SparseMatrix.TransposedMatrixMult(_matrix, s);

         double pr = p * residual;
         double beta = pr / prPrev;
         prPrev = pr;

         z = residual + beta * z;
         s = p + beta * s;

         discrepancy = residual.Norm() / vecNorm;
      }

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
   }
}

public class BCGWithLU : Solver
{
   public override Vector Solve()
   {
      SparseMatrix reMatrixLU = new(_matrix.Size, _matrix._ja.Length);
      SparseMatrix.Copy(_matrix, reMatrixLU);

      PartitialLU(reMatrixLU);
      var reVector = Forward(reMatrixLU, _vector);


      _solution = new(_vector.Size);
      Vector residual = reVector - Forward(reMatrixLU, _matrix * Backward(reMatrixLU, _solution));
      //Vector residual = _vector - _matrix * _solution;

      Vector p = new(residual.Size);
      Vector z = new(residual.Size);
      Vector s = new(residual.Size);

      Vector.Copy(residual, p);
      Vector.Copy(residual, z);
      Vector.Copy(residual, s);

      Stopwatch sw = Stopwatch.StartNew();

      //double vecNorm = _vector.Norm();
      double vecNorm = reVector.Norm();
      double discrepancy = 1;
      double prPrev = p * residual;

      for (int i = 1; i <= MaxIters && discrepancy > Eps; i++)
      {
         var L_AU_z = Forward(reMatrixLU, _matrix * Backward(reMatrixLU, z));
         double alpha = prPrev / (s * L_AU_z);

         _solution = _solution + alpha * z;
         residual = residual - alpha * L_AU_z;

         //var L_ATU_s = SparseMatrix.TransposedMatrixMult(_matrix, s);
         //var L_ATU_s = Backward(reMatrixLU, _matrix * Forward(reMatrixLU, s));
         var L_ATU_s = Forward(reMatrixLU, SparseMatrix.TransposedMatrixMult(_matrix, Backward(reMatrixLU, s)));
         p = p - alpha * L_ATU_s;

         double pr = p * residual;
         double beta = pr / prPrev;
         prPrev = pr;

         z = residual + beta * z;
         s = p + beta * s;

         discrepancy = residual.Norm() / vecNorm;
      }

      _solution = Backward(reMatrixLU, _solution);

      sw.Stop();
      SolvationTime = sw.ElapsedMilliseconds;

      return _solution;
   }

   protected static void PartitialLU(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            int jCol = Matrix._ja[j];
            int jk = Matrix._ia[jCol];
            int k = Matrix._ia[i];

            int sdvig = Matrix._ja[Matrix._ia[i]] - Matrix._ja[Matrix._ia[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix._ia[jCol + 1]; k++, jk++)
            {
               sumL += Matrix._al[k] * Matrix._au[jk];
               sumU += Matrix._au[k] * Matrix._al[jk];
            }

            Matrix._al[j] -= sumL;
            Matrix._au[j] -= sumU;
            Matrix._au[j] /= Matrix._di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
            sumD += Matrix._al[j] * Matrix._au[j];

         Matrix._di[i] -= sumD;
      }
   }

   protected static Vector Forward(SparseMatrix Matrix, Vector b)
   {
      var result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            result[i] -= Matrix._al[j] * result[Matrix._ja[j]];
         }

         result[i] /= Matrix._di[i];
      }

      return result;
   }

   protected static Vector Backward(SparseMatrix Matrix, Vector b)
   {
      var result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix._ia[i + 1] - 1; j >= Matrix._ia[i]; j--)
         {
            result[Matrix._ja[j]] -= Matrix._au[j] * result[i];
         }
      }

      return result;
   }
}