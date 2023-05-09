using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CourseProjectFEM;

public class Point2D
{
   public double X { get; init; }
   public double Y { get; init; }

   public Point2D(double x, double y) 
   {
      X = x;
      Y = y;
   }

   public override string ToString()
   {
      return $"{X:e15} {Y:e15}";
   }

   public Point2D Parse(string input)
   {
      var data = input.Split().Select(double.Parse).ToList(); 
      return new Point2D(data[0], data[1]);
   }
}
