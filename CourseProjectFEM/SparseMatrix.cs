using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProjectFEM;

public class SparseMatrix
{
   public int[] _ia { get; set; }
   public int[] _ja { get; set; }
   public double[] _di { get; set; }
   public double[] _al { get; set; }
   public double[] _au { get; set; }
   public int Size { get; init; }

   public SparseMatrix(int size, int elemsCount)
   {
      Size = size;
      _ia = new int[size + 1];
      _ja = new int[elemsCount];
      _al = new double[elemsCount];
      _au = new double[elemsCount];
      _di = new double[size];
   }

   public static Vector operator *(SparseMatrix matrix, Vector vector)
   {
      Vector product = new(vector.Size);

      for (int i = 0; i < vector.Size; i++)
      {
         product[i] += matrix._di[i] * vector[i];

         for (int j = matrix._ia[i]; j < matrix._ia[i + 1]; j++)
         {
            product[i] += matrix._al[j] * vector[matrix._ja[j]];
            product[matrix._ja[j]] += matrix._au[j] * vector[i];
         }
      }

      return product;
   }

   public static Vector TransposedMatrixMult(SparseMatrix matrix, Vector vector)
   {
      Vector product = new(vector.Size);

      for (int i = 0; i < vector.Size; i++)
      {
         product[i] += matrix._di[i] * vector[i];

         for (int j = matrix._ia[i]; j < matrix._ia[i + 1]; j++)
         {
            product[i] += matrix._au[j] * vector[matrix._ja[j]];
            product[matrix._ja[j]] += matrix._al[j] * vector[i];
         }
      }

      return product;
   }

   public void Clear()
   {
      for (int i = 0; i < Size; i++)
      {
         _di[i] = 0.0;

         for (int k = _ia[i]; k < _ia[i + 1]; k++)
         {
            _al[k] = 0.0;
            _au[k] = 0.0;
         }
      }
   }

   public SparseMatrix ConvertToProfile()
   {
      int sizeOffDiag = 0;
      for (int i = 0; i < Size; i++)
      {
         sizeOffDiag += i - _ja[_ia[i]];
      }

      SparseMatrix result = new(Size, sizeOffDiag);

      result._ia[0] = 0;

      for (int i = 0; i < Size; i++)
      {
         result._di[i] = _di[i];

         int rowSize = i - _ja[_ia[i]];
         result._ia[i + 1] = result._ia[i] + rowSize;

         int jPrev = _ia[i];
         for (int j = result._ia[i]; j < result._ia[i + 1]; j++)
         {
            int col = i - (result._ia[i + 1] - j);
            int colPrev = jPrev < _ia[i + 1] ? _ja[jPrev] : i;
            if (col == colPrev)
            {
               result._al[j] = _al[jPrev];
               result._au[j] = _au[jPrev++];
            }
            else
            {
               result._al[j] = 0;
               result._au[j] = 0;
            }
            result._ja[j] = col;
         }
      }

      return result;
   }

   public static void Copy(SparseMatrix source, SparseMatrix destination)
   {
      for (int i = 0; i < destination.Size + 1; i++)
      {
         destination._ia[i] = source._ia[i];
      }

      for (int i = 0; i < destination.Size; i++)
      {
         destination._di[i] = source._di[i];
      }

      for (int i = 0; i < destination._ja.Length; i++)
      {
         destination._ja[i] = source._ja[i];
         destination._al[i] = source._al[i];
         destination._au[i] = source._au[i];
      }
   }

   public void PrintDenseToConsole()
   {

   }
}