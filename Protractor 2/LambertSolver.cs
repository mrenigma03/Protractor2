using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public static class LambertSolver
{
	public const double TwoPi = 2.0d * Math.PI;
	public const double HalfPi = 0.5d * Math.PI;
	public static double InverseSquareGoldenRatio = (3.0d - Math.Sqrt(5.0d)) * 0.5d; // This is 1 - (1 / phi) == 1/phi^2
	public static double MachineEpsilon = CalculateMachineEpsilon();
	public static double SqrtMachineEpsilon = Math.Sqrt(MachineEpsilon);

	/// <summary>
	/// Find the universal time of the next lauch window from <paramref name="origin"/> to <paramref name="destination"/>.
	/// Limitations: Does not take into account ejection inclination change costs. Does not take into acccount insertion delta-v.
	/// </summary>
	/// <param name="origin">The origin body.</param>
	/// <param name="destination">The destination body. Must have the same <c>referenceBody</c> as <paramref name="origin"/>.</param>
	public static double NextLaunchWindowUT(CelestialBody origin, CelestialBody destination)
	{
		if (origin.referenceBody != destination.referenceBody)
		{
			throw new ArgumentException("Origin and destination bodies must be orbiting the same referenceBody.");
		}

		double now = Planetarium.GetUniversalTime();

		double currentPhaseAngle = CurrentPhaseAngle(origin.orbit, destination.orbit);
		double hohmannPhaseAngle = HohmannPhaseAngle(origin.orbit, destination.orbit);
		double deltaPhaseAngle = (360.0 + currentPhaseAngle - hohmannPhaseAngle) % 360.0;
		if (destination.orbit.semiMajorAxis < origin.orbit.semiMajorAxis)
		{
			deltaPhaseAngle = 360.0 - deltaPhaseAngle;
		}

		double synodicPeriod = Math.Abs(1.0 / (1.0 / origin.orbit.period - 1.0 / destination.orbit.period));
		double estimatedDeparture = now + deltaPhaseAngle / (360.0 / synodicPeriod);

		if (CoplanarOrbits(origin.orbit, destination.orbit) && origin.orbit.eccentricity == 0 && destination.orbit.eccentricity == 0)
		{ // Hohmann transfer
			return estimatedDeparture;
		}
		else
		{
			double gravParameter = origin.referenceBody.gravParameter;
			double hohmannTimeOfFlight = HohmannTimeOfFlight(origin.orbit, destination.orbit);
			double xmin = now;
			double xmax = now + synodicPeriod;
			bool longWay = false;

			Func<Vector2d, Vector2d, Range> boundsFunc = (Vector2d point, Vector2d direction) => {
				double t_xmin = 0, t_xmax = 0;
				if (direction.x != 0)
				{
					if (direction.x > 0)
					{
						t_xmin = (xmin - point.x) / direction.x;
						t_xmax = (xmax - point.x) / direction.x;
					}
					else
					{
						t_xmin = (xmax - point.x) / direction.x;
						t_xmax = (xmin - point.x) / direction.x;
					}

					if (Math.Abs(direction.y / direction.x) < 0.1)
					{ // Direction is > 90% horizontal, don't worry about the y bounds
						return new Range(t_xmin, t_xmax);
					}
				}

				double timeToOpposition = TimeAtOpposition(origin.orbit.getRelativePositionAtUT(point.x), destination.orbit, point.x + 0.5 * hohmannTimeOfFlight) - point.x;

				double ymin = longWay ? 0.95 * timeToOpposition : 0.5 * hohmannTimeOfFlight;
				double ymax = longWay ? 2.0 * hohmannTimeOfFlight : 1.05 * timeToOpposition;
				double t_ymin, t_ymax;
				if (direction.y > 0)
				{
					t_ymin = (ymin - point.y) / direction.y;
					t_ymax = (ymax - point.y) / direction.y;
				}
				else
				{
					t_ymin = (ymax - point.y) / direction.y;
					t_ymax = (ymin - point.y) / direction.y;
				}

				if (Math.Abs(direction.x / direction.y) < 0.1)
				{ // Direction is > 90% vertical, don't worry about the x bounds
					return new Range(t_ymin, t_ymax);
				}
				else
				{
					return new Range(Math.Max(t_xmin, t_ymin), Math.Min(t_xmax, t_ymax));
				}
			};

			Func<Vector2d,Vector3d> deltaVFunc = (Vector2d coords) => {
				double t1 = coords.x;
				double dt = coords.y;
				Vector3d pos1 = origin.orbit.getRelativePositionAtUT(t1);
				Vector3d pos2 = destination.orbit.getRelativePositionAtUT(t1 + dt);
				Vector3d ejectionVelocity = Solve(gravParameter, pos1, pos2, dt, longWay);
				return ejectionVelocity - origin.orbit.getOrbitalVelocityAtUT(t1);
			};

			Vector2d shortTransfer = new Vector2d(now + estimatedDeparture, 0.90 * hohmannTimeOfFlight);
			Vector3d shortDeltaVector = MinimizeDeltaV(ref shortTransfer, boundsFunc, 1e-4, deltaVFunc);

			longWay = true;
			Vector2d longTransfer = new Vector2d(shortTransfer.x, 1.10 * hohmannTimeOfFlight);
			Vector3d longDeltaVector = MinimizeDeltaV(ref longTransfer, boundsFunc, 1e-4, deltaVFunc);

			if (shortDeltaVector.sqrMagnitude <= longDeltaVector.sqrMagnitude)
			{
				return shortTransfer.x;
			}
			else
			{
				return longTransfer.x;
			}
		}
	}

	/// <summary>
	/// Determines the ejection burn parameters to take <paramref name="vessel"/> from its parent body to <paramref name="destination"/>.
	/// Limitations: Does not take into account inclination change delta-v when determining the time of the optimal transfer.
	/// </summary>
	/// <param name="ut">The universal time of the ejection burn.</param>
	/// <param name="vessel">The vessel's orbit is assumed to be (approximately) circular and equatorial.</param>
	/// <param name="destination">The destination body.</param>
	/// <param name="angleToPrograde">The angle to the vessel's parent body's prograde vector where the impulse burn should take place.</param>
	/// <param name="ejectionInclination">The inclination of the hyperbolic ejection orbit.</param>
	/// <param name="deltaV">The total delta-v of the ejection burn, including inclination change.</param>
	public static Vector3d EjectionBurn(ref double ut, Vessel vessel, CelestialBody destination)
	{
		ManeuverNode n = new ManeuverNode();
		n.
		CelestialBody origin = vessel.mainBody;
		if (origin.referenceBody != destination.referenceBody) {
			throw new ArgumentException ("Vessel must be orbiting a body with the same referenceBody as destination.");
		}

		double gravParameter = origin.referenceBody.gravParameter;
		Vector3d pos1 = origin.orbit.getRelativePositionAtUT (ut);
		Vector3d originVelocity = origin.orbit.getOrbitalVelocityAtUT (ut);
		Vector3d originNormal = origin.orbit.GetOrbitNormal ().normalized;
		double hohmannTimeOfFlight = HohmannTimeOfFlight (origin.orbit, destination.orbit);

		Vector3d deltaVector;
		if (CoplanarOrbits (origin.orbit, destination.orbit)) {
			Range bounds = new Range (0.5d * hohmannTimeOfFlight, 2.0d * hohmannTimeOfFlight);

			double timeOfFlight;
			deltaVector = MinimizeDeltaV (bounds, out timeOfFlight, 1e-4, (x) => {
				Vector3d pos2 = destination.orbit.getRelativePositionAtUT (ut + x);
				Vector3d ejectionVelocity = Solve (gravParameter, pos1, pos2, x, (x > hohmannTimeOfFlight));
				return ejectionVelocity - originVelocity; });

			Debug.Log (String.Format ("timeOfFlight, deltaVector, deltaV: {0}, {1}, {2}", timeOfFlight / 3600 / 24, deltaVector, deltaVector.magnitude));
		} else {
			double timeToOpposition = TimeAtOpposition (pos1, destination.orbit, ut + 0.5 * hohmannTimeOfFlight) - ut;
			bool longWay = false;

			Func<double,Vector3d> func = (x) => {
				Vector3d pos2 = destination.orbit.getRelativePositionAtUT (ut + x);
				Vector3d ejectionVelocity = Solve (gravParameter, pos1, pos2, x, longWay);
				return ejectionVelocity - originVelocity;
			};

			Range bounds = new Range (0.5d * hohmannTimeOfFlight, timeToOpposition);
			double shortTimeOfFlight;
			Vector3d shortDeltaVector = MinimizeDeltaV (bounds, out shortTimeOfFlight, 1e-4, func);

			longWay = true;
			bounds.lower = timeToOpposition;
			bounds.upper = 2.0d * hohmannTimeOfFlight;
			double longTimeOfFlight;
			Vector3d longDeltaVector = MinimizeDeltaV (bounds, out longTimeOfFlight, 1e-4, func);

			if (shortDeltaVector.sqrMagnitude <= longDeltaVector.sqrMagnitude) {
				deltaVector = shortDeltaVector;
				Debug.Log (String.Format ("timeOfFlight, deltaVector, deltaV: {0}, {1}, {2}", shortTimeOfFlight / 3600 / 24, deltaVector, deltaVector.magnitude));
			} else {
				deltaVector = longDeltaVector;
				Debug.Log (String.Format ("timeOfFlight, deltaVector, deltaV: {0}, {1}, {2}", longTimeOfFlight / 3600 / 24, deltaVector, deltaVector.magnitude));
			}
		}

		deltaV = deltaVector.magnitude;
		ejectionInclination = Math.Asin (Vector3d.Dot (deltaVector, originNormal) / deltaV) * 180.0d / Math.PI;
		if (vessel.orbit.eccentricity < 1.0) { // Vessel is in an elliptical orbit
			// Assume the orbit is circular and equatorial. Otherwise we'd have to do another minimization search.
			double orbitalSpeed = vessel.orbit.getOrbitalVelocityAtUT (ut).magnitude;
			double vinf = deltaV;
			deltaV = CircularToHyperbolicDeltaV (orbitalSpeed, vinf, ejectionInclination);

			double r = vessel.orbit.getRelativePositionAtUT (ut).magnitude;
			double v = orbitalSpeed + CircularToHyperbolicDeltaV (orbitalSpeed, vinf);
			double e = r * v * v / gravParameter - 1; // Eq. 4.30 simplified for a flight path angle of 0
			angleToPrograde = HyperbolicEjectionAngle (deltaVector, e, originNormal, originVelocity);
		} else {
			Vector3d deltaVInPlane = Vector3d.Exclude (originNormal, deltaVector);
			angleToPrograde = Vector3d.Angle (deltaVInPlane, originVelocity);
			if (Vector3d.Dot (Vector3d.Cross (deltaVInPlane, originVelocity), originNormal) < 0) {
				angleToPrograde = 360.0d - angleToPrograde;
			}
		}
	}

	/// <summary>
	/// Calculates the time of flight for a Hohmann transfer between <paramref name="origin"/> and <paramref name="destination"/>, assuming the orbits are circular and coplanar.
	/// </summary>
	/// <returns>The time of flight.</returns>
	/// <param name="origin">The origin orbit.</param>
	/// <param name="destination">The destination orbit.</param>
	public static double HohmannTimeOfFlight (Orbit origin, Orbit destination)
	{
		double a = (origin.semiMajorAxis + destination.semiMajorAxis) * 0.5;
		double mu = origin.referenceBody.gravParameter;
		return Math.PI * Math.Sqrt ((a * a * a) / mu);
	}

	/// <summary>
	/// Calculates the phase angle for a Hohmann transfer between <paramref name="origin"/> and <paramref name="destination"/>, assuming the orbits are circular and coplanar.
	/// </summary>
	/// <returns>The phase angle.</returns>
	/// <param name="origin">The origin orbit.</param>
	/// <param name="destination">The destination orbit.</param>
	public static double HohmannPhaseAngle (Orbit origin, Orbit destination)
	{
		return 180.0 - HohmannTimeOfFlight (origin, destination) * 360.0 / destination.period;
	}

	/// <summary>
	/// Calculates the current phase angle between <paramref name="origin"/> and <paramref name="destination"/>.
	/// </summary>
	/// <returns>The phase angle.</returns>
	/// <param name="origin">Origin.</param>
	/// <param name="destination">Destination.</param>
	public static double CurrentPhaseAngle (Orbit origin, Orbit destination)
	{
		Vector3d normal = origin.GetOrbitNormal ().normalized;
		Vector3d projected = Vector3d.Exclude (normal, destination.pos);
		double result = Vector3d.Angle (origin.pos, projected);
		if (Vector3d.Dot (Vector3d.Cross (origin.pos, projected), normal) < 0) {
			return 360.0 - result;
		} else {
			return result;
		}
	}

	/// <summary>
	/// Calculates the earliest universal time after <paramref name="after"/> when <paramref name="destination"/> will be 180 degrees from the <paramref name="origin"/> position.
	/// </summary>
	/// <returns>The universal time when <paramref name="destination"/> will be in opposition to <paramref name="origin"/>.</returns>
	/// <param name="origin">Origin position.</param>
	/// <param name="destination">Destination orbit.</param>
	/// <param name="after">Universal time after which to find the opposition.</param>
	public static double TimeAtOpposition (Vector3d origin, Orbit destination, double after = 0)
	{
		Vector3d normal = destination.GetOrbitNormal ().normalized;
		Vector3d periapsis = destination.GetEccVector ();
		Vector3d direction = Vector3d.Exclude (normal, -origin);

		double trueAnomaly = Vector3d.Angle (periapsis, direction) * Mathf.Deg2Rad;
		if (Vector3d.Dot (Vector3d.Cross (periapsis, direction), normal) < 0) {
			trueAnomaly = TwoPi - trueAnomaly;
		}

		double ut = Planetarium.GetUniversalTime () + destination.GetDTforTrueAnomaly (trueAnomaly, 0);
		while (ut <= after) {
			ut += destination.period;
		}

		return ut;
	}

	/// <summary>
	/// Solve the Lambert Problem, determining the velocity at <paramref name="pos1"/> of an orbit passing through <paramref name="pos2"/> after <paramref name="timeOfFlight"/>.
	/// </summary>
	/// <returns>The velocity vector of the identified orbit at <paramref name="pos1"/>.</returns>
	/// <param name="gravParameter">Gravitational parameter of the central body.</param>
	/// <param name="pos1">The first point the orbit passes through.</param>
	/// <param name="pos2">The second point the orbit passes through.</param>
	/// <param name="timeOfFlight">The time of flight between <paramref name="pos1"/> and <paramref name="pos2"/>.</param>
	/// <param name="longWay">If set to <c>true</c>, solve for an orbit subtending more than 180 degrees between <paramref name="pos1"/> and <paramref name="pos2"/>.</param>
	public static Vector3d Solve (double gravParameter, Vector3d pos1, Vector3d pos2, double timeOfFlight, bool longWay)
	{
		// Based on Sun, F.T. "On the Minium Time Trajectory and Multiple Solutions of Lambert's Problem"
		// AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979
		double r1 = pos1.magnitude;
		double r2 = pos2.magnitude;
		double cosAngleOfFlight = Vector3d.Dot (pos1, pos2) / (r1 * r2);
		double angleOfFlight = Math.Acos (cosAngleOfFlight);
		if (longWay) {
			angleOfFlight = TwoPi - angleOfFlight;
		}

		// Intermediate terms
		Vector3d deltaPos = pos2 - pos1;
		double c = deltaPos.magnitude;
		double m = r1 + r2 + c;
		double n = r1 + r2 - c;

		double cosHalfAngleOfFlight = Math.Cos (0.5 * angleOfFlight);
		double angleParameter = Math.Sqrt (4.0 * r1 * r2 / (m * m) * cosHalfAngleOfFlight * cosHalfAngleOfFlight);
		if (longWay) {
			angleParameter = -angleParameter;
		}

		double normalizedTime = 4.0 * timeOfFlight * Math.Sqrt (gravParameter / (m * m * m));
		double parabolicNormalizedTime = 2.0 / 3.0 * (1 - angleParameter * angleParameter * angleParameter);
		double minimumEnergyNormalizedTime = Math.Acos (angleParameter) + angleParameter * Math.Sqrt (1 - angleParameter * angleParameter);

		double x, y; // Path parameters
		Func<double,double> fy = (xn) => (angleParameter < 0) ? -Math.Sqrt (1.0 - angleParameter * angleParameter * (1.0 - xn * xn)) : Math.Sqrt (1.0 - angleParameter * angleParameter * (1.0 - xn * xn));
		if (normalizedTime == parabolicNormalizedTime) { // Parabolic solution
			x = 1.0;
			y = (angleParameter < 0) ? -1 : 1;
		} else if (normalizedTime == minimumEnergyNormalizedTime) { // Minimum energy elliptical solution
			x = 0.0;
			y = fy (x);
		} else {
			// Returns the difference between the normalized time for a path parameter of xn and normalizedTime
			Func<double,double> fdt = (xn) => {
				if (xn == 1.0) { // Parabolic
					return parabolicNormalizedTime - normalizedTime;
				} else {
					double yn = fy(xn);

					double g, h;
					if (xn > 1.0) { // Hyperbolic
						g = Math.Sqrt (xn * xn - 1.0);
						h = Math.Sqrt (yn * yn - 1.0);
						return (Acoth (yn / h) - Acoth (xn / g) + xn * g - yn * h) / (g * g * g) - normalizedTime;
					} else { // Elliptical (-1 < x < 1)
						g = Math.Sqrt (1.0 - xn * xn);
						h = Math.Sqrt (1.0 - yn * yn);
						return (Acot (xn / g) - Math.Atan(h / yn) - xn * g + yn * h) / (g * g * g) - normalizedTime;
					}
				}
			};

			// Select our bounds based on the relationship between the known normalized times and normalizedTime
			Range bounds;
			if (normalizedTime > minimumEnergyNormalizedTime) { // Elliptical high path solution
				bounds.lower = -1.0 + MachineEpsilon;
				bounds.upper = 0.0;
			} else if (normalizedTime > parabolicNormalizedTime) { // Elliptical low path solution
				bounds.lower = 0.0;
				bounds.upper = 1.0;
			} else { // Hyperbolic solution
				bounds.lower = 1.0;
				bounds.upper = 2.0;
				while (fdt(bounds.upper) > 0.0) {
					bounds.lower = bounds.upper;
					bounds.upper *= 2.0;
				}
			}

			x = FindRoot (bounds, 1e-4, fdt); // Solve for x
			y = fy (x);
		}

		double sqrtMu = Math.Sqrt (gravParameter);
		double invSqrtM = 1.0 / Math.Sqrt (m);
		double invSqrtN = 1.0 / Math.Sqrt (n);

		double vc = sqrtMu * (y * invSqrtN + x * invSqrtM);
		double vr = sqrtMu * (y * invSqrtN - x * invSqrtM);
		Vector3d ec = deltaPos * (vc / c);
		return ec + pos1 * (vr / r1);
	}

	private struct Range
	{
		public double lower, upper;

		public Range (double lwr, double upr)
		{
			lower = lwr;
			upper = upr;
		}
	}

	private static Vector3d MinimizeDeltaV (ref Vector2d p0, Func<Vector2d, Vector2d, Range> getBounds, double relativeAccuracy, Func<Vector2d, Vector3d> f)
	{
		// Uses Powell's method to find the local minimum of f(x,y) within the bounds returned by getBounds(point, direction): http://en.wikipedia.org/wiki/Powell's_method
		Queue<Vector2d> directionVectors = new Queue<Vector2d> (new Vector2d[] { new Vector2d (1, 0), new Vector2d (0, 1) });
		double sqrRelativeAccuracy = relativeAccuracy * relativeAccuracy;
		Vector3d result = new Vector3d();

		Func<Vector2d, Vector2d, Vector2d> findMinimumAlongDirection = (p, direction) => {
			double u;
			result = MinimizeDeltaV (getBounds (p, direction), out u, relativeAccuracy, (v) => {
				Vector2d point = p + v * direction;
				return f (point); });
			return p + u * direction;
		};

		for (int i = 1; i < 100; i++) {
			Vector2d pn = p0;
			foreach (Vector2d direction in directionVectors) {
				pn = findMinimumAlongDirection (pn, direction);
			}

			directionVectors.Dequeue ();
			directionVectors.Enqueue (pn - p0);

			pn = findMinimumAlongDirection (p0, directionVectors.Last ());
			double sqrDistance = (pn - p0).sqrMagnitude;
			p0 = pn;

			if (sqrDistance <= p0.sqrMagnitude * sqrRelativeAccuracy) {
				return result;
			}
		}

		throw new Exception ("LambertSolver 2D delta-v minimization failed to converge!");
	}

	private static Vector3d MinimizeDeltaV (Range bounds, out double x, double relativeAccuracy, Func<double,Vector3d> f)
	{
		// Uses Brent's method of parabolic interpolation to find a local minimum: http://linneus20.ethz.ch:8080/1_5_2.html
		x = bounds.lower + InverseSquareGoldenRatio * (bounds.upper - bounds.lower);
		double w = x;
		double v = w;
		double e = 0.0;
		Vector3d fxVector = f (x);
		double fx = fxVector.sqrMagnitude;
		double fw = fx;
		double fv = fw;
		double delta = 0;

		for (int i = 0;; i++) {
			double midpoint = 0.5d * (bounds.lower + bounds.upper);
			double tol = (SqrtMachineEpsilon + relativeAccuracy) * Math.Abs (x);
			double t2 = 2.0d * tol;

			if (Math.Abs (x - midpoint) <= t2 - 0.5d * (bounds.upper - bounds.lower)) { // Are we close enough?
				return fxVector;
			} else if (i > 100) {
				throw new Exception ("LambertSolver 1D delta-v minimization failed to converge!");
			}

			// Fit a parabola between a, x, and b 
			double p = 0, q = 0, r = 0;
			if (tol < Math.Abs (e)) {
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fv);
				p = (x - v) * q - (x - w) * r;
				q = 2.0d * (q - r);
				if (q <= 0.0) {
					q = -q;
				} else {
					p = -p;
				}
				r = e;
				e = delta;
			}

			double u;
			if (Math.Abs (p) < Math.Abs (0.5d * q * r) && p > q * (bounds.lower - x) && p < q * (bounds.upper - x)) {
				// Use parabolic interpolation for this step
				delta = p / q;
				u = x + delta;

				// We don't want to evaluate f for x within 2 * tol of a or b
				if ((u - bounds.lower) < t2 || (bounds.upper - u) < t2) {
					if (x < midpoint) {
						delta = tol;
					} else {
						delta = -tol;
					}
				}
			} else {
				// Use the golden section for this step
				if (x < midpoint) {
					e = bounds.upper - x;
				} else {
					e = bounds.lower - x;
				}
				delta = InverseSquareGoldenRatio * e;
					
				if (Math.Abs (delta) >= tol) {
					u = x + delta;
				} else if (delta > 0.0) {
					u = x + tol;
				} else {
					u = x - tol;
				}
			}

			Vector3d fuVector = f (u);
			double fu = fuVector.sqrMagnitude;

			if (fu <= fx) {
				if (u < x) {
					bounds.upper = x;
				} else {
					bounds.lower = x;
				}

				v = w;
				fv = fw;
				w = x;
				fw = fx;
				x = u;
				fxVector = fuVector;
				fx = fu;
			} else {
				if (u < x) {
					bounds.lower = u;
				} else {
					bounds.upper = u;
				}

				if (fu <= fw || w == x) {
					v = w;
					fv = fw;
					w = u;
					fw = fu;
				} else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}
	}

	private static double FindRoot (Range bounds, double relativeAccuracy, Func<double,double> f)
	{
		// Uses Brent's root finding method: http://math.fullerton.edu/mathews/n2003/BrentMethodMod.html
		double a = bounds.lower;
		double b = bounds.upper;
		double c = a;
		double fa = f (a);
		double fb = f (b);
		double fc = fa;
		double d = b - a;
		double e = d;

		if (fa >= 0) {
			return a;
		} else if (fb <= 0) {
			return b;
		}

		for (int i = 0;; i++) {
			if (Math.Abs (fc) < Math.Abs (fb)) { // If c is closer to the root than b, swap b and c
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}

			double tol = (0.5d * MachineEpsilon + relativeAccuracy) * Math.Abs (b);
			double m = 0.5d * (c - b);

			if (fb == 0 || Math.Abs (m) <= tol) {
				return b;
			} else if (i > 100) {
				throw new Exception ("LambertSolver root failed to converge!");
			}

			if (Math.Abs (e) < tol || Math.Abs (fa) <= Math.Abs (fb)) { // Use a bisection step
				d = e = m;
			} else {
				double p, q, r;
				double s = fb / fa;

				if (a == c) { // Use a linear interpolation step
					p = 2 * m * s;
					q = 1 - s;
				} else {  // Use inverse quadratic interpolation
					q = fa / fc;
					r = fb / fc;
					p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
					q = (q - 1) * (r - 1) * (s - 1);
				}

				if (p > 0) {
					q = -q;
				} else {
					p = -p;
				}

				if (2 * p < Math.Min (3 * m * q - Math.Abs (tol * q), Math.Abs (e * q))) { // Verify interpolation
					e = d;
					d = p / q;
				} else { // Fall back to bisection
					d = e = m;
				}
			}

			a = b;
			fa = fb;

			if (Math.Abs (d) > tol) {
				b += d;
			} else {
				b += (m > 0 ? tol : -tol);
			}

			fb = f (b);

			if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) { // Ensure fb and fc have different signs
				c = a;
				fc = fa;
				d = e = b - a;
			}
		}
	}

	private static double CircularToHyperbolicDeltaV (double v0, double vinf, double relativeInclination = 0.0)
	{
		double v1 = Math.Sqrt (vinf * vinf + 2 * v0 * v0);
		if (relativeInclination != 0) {
			return Math.Sqrt (v0 * v0 + v1 * v1 - 2 * v0 * v1 * Math.Cos (relativeInclination * Math.PI / 180.0d)); // Eq. 4.74
		} else {
			return v1 - v0;
		}
	}

	private static double HyperbolicEjectionAngle (Vector3d vinf, double eccentricity, Vector3d normal, Vector3d prograde)
	{
		vinf = vinf.normalized;
		prograde = prograde.normalized;

		// We have three equations of three unknowns (v.x, v.y, v.z):
		//   dot(v, vinf) = cos(eta) = -1 / e  [Eq. 4.81]
		//   norm(v) = 1  [Unit vector]
		//   dot(v, normal) = 0  [Perpendicular to normal]
		//
		// Solution is defined iff:
		//   normal.z != 0
		//   vinf.y != 0 or (vinf.z != 0 and normal.y != 0) [because we are solving for v.x first]
		//   vinf is not parallel to normal

		// Intermediate terms
		double f = vinf.y - vinf.z * normal.y / normal.z;
		double g = (vinf.z * normal.x - vinf.x * normal.z) / (vinf.y * normal.z - vinf.z * normal.y);
		double h = (normal.x + g * normal.y) / normal.z;
		double m = normal.y * normal.y + normal.z * normal.z;
		double n = eccentricity * f * normal.z * normal.z;

		// Quadratic coefficients
		double a = (1 + g * g + h * h);
		double b = -2 * (g * m + normal.x * normal.y) / n;
		double c = m / (eccentricity * f * n) - 1;

		Vector3d v;
		v.x = (-b + Math.Sqrt (b * b - 4 * a * c)) / (2 * a);
		v.y = g * v.x - 1 / (eccentricity * f);
		v.z = -(v.x * normal.x + v.y * normal.y) / normal.z;

		if (Vector3d.Dot (Vector3d.Cross (v, vinf), normal) < 0) { // Wrong orbital direction
			v.x = (-b - Math.Sqrt (b * b - 4 * a * c)) / (2 * a);
			v.y = g * v.x - 1 / (eccentricity * f);
			v.z = -(v.x * normal.x + v.y * normal.y) / normal.z;
		}

		if (Vector3d.Dot (Vector3d.Cross (v, prograde), normal) < 0) {
			return 360 - Vector3d.Angle (v, prograde);
		} else {
			return Vector3d.Angle (v, prograde);
		}
	}

	private static bool CoplanarOrbits (Orbit o1, Orbit o2)
	{
		return o1.inclination == o2.inclination && (o1.inclination == 0 || o1.LAN == o2.LAN);
	}

	private static double Acot (double x)
	{
		return HalfPi - Math.Atan (x);
	}
	
	private static double Acoth (double x)
	{
		return 0.5 * Math.Log ((x + 1) / (x - 1));
	}

	private static double CalculateMachineEpsilon ()
	{
		double result = 1.0d;

		do {
			result /= 2.0d;
		} while (1.0d + result != 1.0d);

		return result;
	}
}

