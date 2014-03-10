/*
 * geometry.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "interval.h"

namespace plif {

/**
 * Euclidean distance between two points
 *
 */
inline
precision_type length(breakpoint a, breakpoint b)
{
	precision_type u = a.get_y()-b.get_y();
	precision_type v = a.get_x()-b.get_x();

	return sqrt(u*u+v*v);
}

/**
 * Computes the intersection of two line segments and returns the x-coordinate of the intersection point
 *
 * In case of bool left true
 	* a = -infinity, b = x co-ordinate of left most point
 	* m = left slope of 1st function, n = left limit of 1st function at b
 	* p = left slope of 2nd function, q = left limit of 2nd function at b
 * In case of bool right true
 	* a = x co-ordinate of right most point, b = infinity
 	* m = right limit of 1st function at b, n = right slope of 1st function
 	* p = right limit of 2nd function at b, q = right slope of 2nd function
 * In case of both left and right false, uses two point form of line
 	* a = x-coordinate of 1st point, b = x-coordinate of 2nd point
 	* m = y-coordinate of 1st point for 1st function, b = y-coordinate of 2nd point for 1st function
 	* p = y-coordinate of 1st point for 2nd function, q = y-coordinate of 2nd point for 2nd function
 * The intersection values used here are calculated using the line equations constructed from the given information
 * left and right cannot be simultaneously true 	.
 */
inline
precision_type intersection(precision_type a, precision_type b, precision_type m, precision_type n, precision_type p, precision_type q, bool left, bool right)
{
	precision_type ret_val;
	if(left)
	{
		if(p==m)
			ret_val = NEG_INFTY;
		else
			ret_val = b - (q-n)/(p-m);
	}
	else if(right)
	{
		if(q==n)
			ret_val = POS_INFTY;
		else
			ret_val = a + (m-p)/(q-n);
	}
	else
		ret_val = (a*(q-n) + b*(m-p))/(m-p+q-n);
	return ret_val;
}

/* Proper intersection algorithm that should be used:
/// <summary>
	/// Find the intersection point between two lines.
	/// </summary>
	/// <param name="IntersectPoint">The intersection point. A <see cref="Esteem.Geometry.PointD">PointD</see> structure.</param>
	/// <param name="L1StartPoint">The starting point of first line. A PointD structure.</param>
	/// <param name="L1EndPoint">The end point of first line. A PointD structure.</param>
	/// <param name="L2StartPoint">The starting point of second line. A PointD structure.</param>
	/// <param name="L2EndPoint">The end point of second line. A PointD structure.</param>
	/// <param name="L1IntersectPos">The intersection position at first line.</param>
	/// <param name="L2IntersectPos">The intersection position at second line.</param>
	/// <returns>Returns a boolean. True if there is intersection, false otherwise.</returns>
	/// <remarks>The formula is taken from comp.graphics.algorithms Frequently Asked Questions.</remarks>
	public static bool LineIntersect(out PointD IntersectPoint, PointD L1StartPoint, PointD L1EndPoint, PointD L2StartPoint, PointD L2EndPoint, out double L1IntersectPos, out double L2IntersectPos)
	{
		IntersectPoint = new PointD();
		double ay_cy, ax_cx, px, py;
		double dx_cx = L2EndPoint.X - L2StartPoint.X,
			dy_cy = L2EndPoint.Y - L2StartPoint.Y,
			bx_ax = L1EndPoint.X - L1StartPoint.X,
			by_ay = L1EndPoint.Y - L1StartPoint.Y;

		double de = (bx_ax) * (dy_cy) - (by_ay) * (dx_cx);
		//double tor = 1.0E-10;		//tolerance


		L1IntersectPos = 0;		 L2IntersectPos = 0;
		if (Math.Abs(de)<0.01)
			return false;
		//if (de > -tor && de < tor) return false; //line is in parallel

		ax_cx = L1StartPoint.X - L2StartPoint.X;
		ay_cy = L1StartPoint.Y - L2StartPoint.Y;
		double r = ((ay_cy) * (dx_cx) - (ax_cx) * (dy_cy)) / de;
		double s = ((ay_cy) * (bx_ax) - (ax_cx) * (by_ay)) / de;
		px = L1StartPoint.X + r * (bx_ax);
		py = L1StartPoint.Y + r * (by_ay);

		IntersectPoint.X = px;  //return the intersection point
		IntersectPoint.Y = py;	//return the intersection position
		L1IntersectPos = r;		 L2IntersectPos = s;

		return true; //indicate there is intersection
	}
	*/

}

#endif /* GEOMETRY_H_ */
