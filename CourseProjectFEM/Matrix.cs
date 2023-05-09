using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProjectFEM;

public class Matrix
{
   private double[,] _container;
   public int Size { get; }

   public double this[int i, int j]
   {
      get => _container[i, j];
      set => _container[i, j] = value;
   }

   public Matrix(int size)
   {
      _container = new double[size, size];
      Size = size;
   }

   public void Clear()
       => Array.Clear(_container, 0, _container.Length);

   public static void Copy(Matrix source, Matrix destination)
   {
      for (int i = 0; i < destination.Size; i++)
         for (int j = 0; j < destination.Size; j++)
            destination[i, j] = source[i, j];
   }

   public static Matrix operator +(Matrix fstMatrix, Matrix sndMatrix)
   {
      Matrix resultMatrix = new(fstMatrix.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
         for (int j = 0; j < resultMatrix.Size; j++)
            resultMatrix[i, j] = fstMatrix[i, j] + sndMatrix[i, j];

      return resultMatrix;
   }

   public static Matrix operator *(double coef, Matrix Matrix)
   {
      Matrix resultMatrix = new(Matrix.Size);

      for (int i = 0; i < resultMatrix.Size; i++)
         for (int j = 0; j < resultMatrix.Size; j++)
            resultMatrix[i, j] = coef * Matrix[i, j];

      return resultMatrix;
   }
}
