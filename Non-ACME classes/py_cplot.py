from __future__ import division
import sys
import numpy as np
from scipy import misc


with open('sederstring.eps') as f:
    sederstring = f.read()
    
width = 1.0

#Point is np.array([x,y,w])
class Curve(object):
    def __init__(self, points):
        self.degree = len(points) - 1
        self.points = points

class cubic_NURBS(object):
    def __init__(self, knots, points):
        self.degree = 3
        self.knots = knots
        self.points = points

def plotBezier(curve):
    pointlines = ["[{} {} {}]".format(*tuple(p)) for p in curve.points]
    return "[" + "".join(pointlines) + "] cplot"

def plotControlPolygon(curve):
    pointlines = ["{} {} mv ".format(*tuple(magicPoint(curve.points[0])))]
    for pt in curve.points[1:]:
        pointlines.append("{} {} ln ".format(*tuple(magicPoint(pt))))
    pointlines.append("stroke")
    return "".join(pointlines)

def subdivide(curves, params):
    """
    Subdivide ncurve at parameter value t using deCasteljua
    algorithm, assuming original domain is [0,1]. Store left portion as
    curve number nleft and right portion as nright.
    """

    ncurve, t, nleft, nright = params[0], params[1], params[2], params[3]
    pts, n = curves[ncurve].points, curves[ncurve].degree

    if pts.shape[0] == n + 1:
        P = pts.T
    else:
        P = pts

    Q = np.copy(P)
    Q0 = np.copy(P)
    
    for k in xrange(1,n+1):
        for i in xrange(n-k+1):
            Q[:,i] = (1.0-t)*Q[:,i] + t*Q[:,i+1]
        Q0[:,k] = Q[:,0]

    
    if Q.shape[1] == n+1: 
        Q = Q.T
    if Q0.shape[1] == n+1:
        Q0 = Q0.T

    curves[nleft] = Curve(Q0)
    curves[nright] = Curve(Q)
    
    return curves

def elevate(curves, params):
    """
    Degree elevate input curve n_elevate times. Store in new_curve.
    """
    ncurve, n_elevate, n_new = params[0], params[1], params[2]
    pts, n = curves[ncurve].points, curves[ncurve].degree
    P = pts

    for i in xrange(int(n_elevate)):
        alpha = 1.0*np.arange(n+2)/(n+1)
        Q = np.zeros((n+2,P.shape[1]))
        Q[0,:] = (1-alpha[0])*P[0,:]
        for j in xrange(1,n+1):
            Q[j,:] = alpha[j]*P[j-1,:] + (1-alpha[j])*P[j,:]
        Q[n+1,:] = alpha[n+1]*P[n,:]

        P = Q
        n += 1

    curves[n_new] = Curve(P)

    return curves

def magicPoint(pt):
    return pt[0]/pt[2], pt[1]/pt[2]

def readPoint(pt):
    return pt[0]*pt[2], pt[1]*pt[2], pt[2]

def getCurvature(ncurve, t, curves, pt_only=False):
    # Get P(t) and center of curvature at t for curve # ncurve
    # Also return value of curvature k


    n = curves[ncurve].degree
    pts =  curves[ncurve].points

    # Subdivide curve. Want right hand curve, which is Q
    if pts.shape[0] == n + 1:
        P = pts.T
    else:
        P = pts

    Q = np.copy(P)
    Q = np.array(Q,dtype=np.float)
    Q0 = np.copy(P)
    
    for k in xrange(1,n+1):
        for i in xrange(n-k+1):
            Q[:,i] = (1.0-t)*Q[:,i] + t*Q[:,i+1]
        Q0[:,k] = Q[:,0]

    
    if Q.shape[1] == n+1: 
        Q = Q.T # So Q is shape n+1, 3
    if Q0.shape[1] == n+1:
        Q0 = Q0.T

    if t == 1.0:
        Q = Q0[::-1,:]
    curve = (Q.T/Q[:,-1]).T[:,:-1]

    if (curve[0,:] == curve[1,:]).all():
        print "P0 = P1"
        print "t", t

    P0 = np.array(curve[0,:])
    P1 = np.array(curve[1,:])
    if curve.shape[0] > 2:
        P2 = np.array(curve[2,:])
    else:
        P2 = np.zeros(2)

    a2 = np.sum((P1-P0)**2)

    if a2 == 0:
        print "Repeated knot", P0, P1
        a2 = 1
    
    nvec = (P1 - P0)/np.sqrt(a2)
    nvec = nvec[::-1]
    # To figure out which side it goes on, take the dot product with P2-P1 
    # (Positive if you want center of curvature, negative for curvature comb)
    nvec[0] = -1*nvec[0]
    #if np.dot(nvec,P2-P1) < 0:
    #    print "need to flip vector"
    #    nvec *= -1
    #    print np.dot(nvec,P2-P1)

    #if t == 1:
    #    nvec *= -1

    if pt_only:
        return P0, -1*nvec, 0

    h = np.linalg.det(np.vstack((P1-P0,P2-P0)))/np.sqrt(a2)

    K = (n-1.0)*h/(n*a2)

    cc = P0 + nvec/K

    return P0, cc, K


def main():
    if len(sys.argv) < 3:
        
        print "Syntax:\n python py_cyplot.py testfile outputfile"
        return

    curves = {}

    NURBScurves = {}

    output = [sederstring]

    with open(sys.argv[1]) as f:
        while True:
            line = f.readline().strip().lower()[:4]

            if line == "bord":
                output.append("border stroke")

            elif line == "cplo":
                    
                output.append(plotBezier(curves[nCurve]))

            elif line == "circ":
                params = map(float, f.readline().strip().lower().split())
                output.append("{} {} {} circ".format(*params))

            elif line == "colo":
                params = map(float, f.readline().strip().lower().split())
                output.append("{} {} {} setrgbcolor".format(*params))

            elif line == "cppl":
                nCurve = int(f.readline().strip().lower())
                output.append(plotControlPolygon(curves[nCurve]))

            elif line == "exit":
                break

            elif line == "stor":
                nCurve = int(f.readline().strip().lower())
                degree = int(f.readline().strip().lower())
                pts = [readPoint(map(float, f.readline().strip().lower().split())) for v in xrange(degree+1)]
                curves[nCurve] = Curve(np.array(pts))
                print "Stored curve number ", nCurve
                print curves.keys()

            elif line == "text":
                size, x, y = map(float, f.readline().strip().lower().split())
                text = f.readline().strip()
                output.append("{} {} mv /Times-Roman findfont {} scalefont setfont ({}) show".format(x,y,size,text))

            elif line == "view":
                params = map(float, f.readline().strip().lower().split())
                output.append("{} {} {} {} viewport".format(*params))

            elif line == "wind":
                params = map(float, f.readline().strip().lower().split())
                width = min(abs(int(params[0])-int(params[2])),abs(int(params[1])-int(params[3])))
                output.append("{} {} {} {} window".format(*params))

            elif line == "widt":
                size = int(f.readline().strip().lower())
                output.append("{} setlinewidth".format(size))

            elif line == "disp":
                for i in sorted(curves):
                    c = curves[i]
                    for p in c.points:
                        pass

            elif line == "subd":
                params = map(float, f.readline().strip().lower().split())
                curves = subdivide(curves, params)

            elif line == "elev":
                params = map(float, f.readline().strip().lower().split())
                curves = elevate(curves, params)

            elif line == "ccur":
                ncurve, n = map(float, f.readline().strip().lower().split())

                for i in xrange(int(n)+1):
                    # Get P(t) and center of curvature. 
                    output.append("newpath")
                    Pt, cc, k = getCurvature(ncurve, 1.0*i/n, curves)
                    #print "{} {} moveto".format(*Pt)
                    #print "{} {} lineto".format(*cc)
                    output.append("{} {} mv".format(*Pt))
                    output.append("{} {} ln".format(*cc))
                    output.append("closepath")
                    output.append("stroke")

            elif line == "curv":
                # For curve P w/ # ncurve, print val of curvature next to P(t)
                ncurve, t = map(float, f.readline().strip().lower().split())

                Pt, cc, k = getCurvature(ncurve, t, curves)
                print Pt
                print k
                output.append("/Times-Roman findfont 12 scalefont setfont")
                output.append("0.0 0.0 0.0 setrgbcolor")
                output.append("{} {} mv".format(*(Pt)))
                #output.append("() show")
                output.append("({}) show".format(k))

            elif line == "comb":
                ncurve = map(float, f.readline().strip().lower().split())
                print "comb n:", ncurve
                # Draw a curvature comb for ncurve w/ about 100 line segments.
                # Find maximum curvature by sampling curvature at 100 pts, assign
                # maximum length for line segments to be .2*width of current window

                Pt_vals = []
                cc_vals = []
                k_vals = []
                for i in xrange(100):
                    Pt, cc, k = getCurvature(ncurve[0], i/100.0, curves)
                    Pt_vals.append(Pt)
                    cc_vals.append(cc)
                    k_vals.append(k)

                print k_vals[0], k_vals[-1]

                kmax = max(np.abs(k_vals[1:-1]))
                print "kmax", kmax
                scale = .5*width/kmax
                print "scale", scale
                for i in xrange(1,100):
                    if k_vals[i] == 0:
                        print i
                    #comb_pt = k_vals[i]*scale*(2*cc_vals[i] - Pt_vals[i])
                    comb_pt = Pt_vals[i] + k_vals[i]**2*scale*(Pt_vals[i]-cc_vals[i])
                    #comb_pt = Pt_vals[i] - scale*k_vals[i]**2*cc_vals[i]
                    output.append("newpath")
                    output.append("{} {} mv".format(*Pt_vals[i]))
                    output.append("{} {} ln".format(*comb_pt))
                    output.append("closepath")
                    output.append("stroke")

            elif line == "expl":
                print "Get explicit curve"
                ncurve, n = map(float, f.readline().strip().lower().split())
                # f0, f1, ... , fn
                print "Get fvals"
                fvals = map(float, f.readline().strip().lower().split())
                n = len(fvals)-1
                # Store in ncurve a degree n explicit Bezier curve equivalent to
                # x = t, y = f0 + f1*t + f2*t^2 + ... + fn*t^n
                
                
                yvals = fvals/misc.comb(n,np.arange(n+1))
                for i in xrange(1,len(yvals)):
                    for j in xrange(n-1,i-1,-1):
                        yvals[j] += yvals[j-1]
                pts = []
                print "Get pts"
                for i in xrange(n+1):
                    pts.append([1.0*i/n,yvals[i],1])

                print pts
                curves[ncurve] = Curve(pts)
                #print curves[ncurve]

            elif line == "offs":
                ncurve, rho = map(float, f.readline().strip().lower().split())
                # plot the offset of the specified rational Bezier curve.
                # radius > 0 => offset should be on left of tangent vector
                
                Pt_vals = []
                n_vals = []
                # Sample curve at 100 pts
                for i in xrange(101):
                    Pt, nvec, k = getCurvature(ncurve, i/100., curves, pt_only=True)
                    Pt_vals.append(Pt)
                    n_vals.append(nvec)

                i = 0
                output.append("newpath") 
                offset_pt = Pt_vals[i] + rho*n_vals[i]
                output.append("{} {} mv".format(*offset_pt))
                for i in xrange(1,100):

                    offset_pt = Pt_vals[i] + rho*n_vals[i]
                    
                    output.append("{} {} ln".format(*offset_pt))
                    output.append("closepath")
                    output.append("stroke")
                    output.append("newpath")
                    output.append("{} {} mv".format(*offset_pt))

                output.append("{} {} ln".format(*offset_pt))
                output.append("closepath")
                output.append("stroke")

            elif line == "stcn":
                # Store a cubic NURBS curve
                # n, m (curve number and number of control pts)
                # k1, k2, ... , km+2 (knots)
                # x1, y1, w1
                # ...
                # xm, ym, wm (control pts w/ weights)
                n, m = map(int, f.readline().strip().lower().split())
                knots = map(float, f.readline().strip().lower().split())
                pts = [readPoint(map(float, f.readline().strip().lower().split())) for v in xrange(m)]
                NURBScurves[n] = cubic_NURBS(np.array(knots), np.array(pts))

            elif line == "pncp":
                # Plot the control polygon for cubic NURBS n
                # n (curve number)
                nCurve = int(f.readline().strip().lower())
                output.append(plotControlPolygon(NURBScurves[nCurve]))

            elif line == "nebz":
                # For cubic NURBS curve with number n, extract the ith
                # Bezier curve and store it in Bezier curve number m
                # n, i, m
                n, i, m = map(int, f.readline().strip().lower().split())

                k = NURBScurves[n].knots
                pts = NURBScurves[n].points

                # FIX THIS
                c1 = (k[i+1]-k[i-1])/(k[i+2]-k[i-1])
                c2 = (k[i+1] - k[i])/(k[i+3] - k[i])
                c3 = (k[i+2] - k[i])/(k[i+3] - k[i])
                c4 = (k[i+2]-k[i+1])/(k[i+4]-k[i+1])

                Q1 = (1-c1)*pts[i-1,:] + c1*pts[i,:]
                Q2 = (1-c2)*pts[i,:] + c2*pts[i+1,:]
                Q3 = (1-c3)*pts[i,:] + c3*pts[i+1,:]
                Q4 = (1-c4)*pts[i+1,:] + c4*pts[i+2,:]

                q1 = (k[i+1] - k[i])/(k[i+2] - k[i])
                q4 = (k[i+2]-k[i+1])/(k[i+3]-k[i+1])

                Q1 = q1*Q1 + (1-q1)*Q2
                Q4 = q4*Q3 + (1-q4)*Q4

                curves[m] = Curve(np.vstack((Q1,Q2,Q3,Q4)))



            elif line == "ikno":
                # Insert a knot into NURBS n at knot value t. 
                # Store the result in NURBS m
                # n, t, m
                n, t, m = map(int, f.readline().strip().lower().split())

                # Update knot vector and get index of new knot
                knots = NURBScurves[n].knots
                pts = NURBScurves[n].points
                new_knots = np.zeros(len(knots)+1)
                i = 0
                while knots[i] < t:
                    new_knots[i] = knots[i]
                    i += 1
                new_knots[i] = t
                for j in xrange(i+1,len(knots)+1):
                    new_knots[j] = knots[j-1]

                # The knot is inserted between index i-1 and i
                new_pts = np.zeros((len(pts)+1,3))
                for j in xrange(i-2):
                    new_pts[j,:] = pts[j,:]

                for j in xrange(i-2,i+1):
                    # j = i-2, i-1, i
                    c = (t - knots[j-1])/(knots[j+2] - knots[j-1])
                    new_pts[j,:] = c*pts[j,:] + (1-c)*pts[j-1,:]

                for j in xrange(i+1,len(pts)+1):
                    new_pts[j,:] = pts[j-1,:]

                NURBScurves[m] = cubic_NURBS(new_knots, new_pts)
                


            else:
                raise Exception("Illegal command: {}".format(line))


        with open(sys.argv[2], 'w') as f:
            f.write("\n".join(output)+"\n")

if __name__ == "__main__":
    main()