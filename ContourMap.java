package com.sawyer.contour;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;

public class ContourMap {
    private static final Logger log = LoggerFactory.getLogger(ContourMap.class);

    private final static double CLOSED = -Double.MAX_VALUE;

    private int w;
    private int h;
    private double[] grid;

    public ContourMap(int width, int height, double[] grid) {
        this.w = width;
        this.h = height;
        this.grid = grid;
    }

    public ContourMap closed() {
        int w = this.w + 2;
        int h = this.h + 2;

        double[] grid = new double[w * h];

        for (int i = 0; i < grid.length; i++) {
            grid[i] = CLOSED;
        }

        for (int y = 0; y < this.h; y++) {
            int i = (y + 1) * w + 1;
            int j = y * this.w;

            System.arraycopy(this.grid, j, grid, i, this.w);
        }

        return new ContourMap(w, h, grid);
    }

    public Map<Double, List<List<double[]>>> contours(double[] zs) {
        long now = System.currentTimeMillis();

        List<Callable<List<List<double[]>>>> tasks = new ArrayList<>(zs.length);
        for (double z : zs) {
            tasks.add(() -> marchingSquares(this.w, this.h, z));
        }

        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        Map<Double, List<List<double[]>>> ret = new HashMap<>(zs.length);

        try {
            List<Future<List<List<double[]>>>> results = executor.invokeAll(tasks);

            executor.shutdown();

            int index = 0;
            for (double z : zs) {
                try {
                    ret.put(z, results.get(index).get());
                } catch (ExecutionException e) {
                    log.error("parse contour level[" + z + "] with error", e);
                }

                index++;
            }
        } catch (Exception e) {
            log.error("contours error", e);
        }

        log.info("contours elapsed: " + (System.currentTimeMillis() - now));

        return ret;
    }

    public List<List<double[]>> contours(double z) {
        return marchingSquares(this.w, this.h, z);
    }

    private List<List<double[]>> marchingSquares(int w, int h, double z) {
        Map<Edge, Point> edgePoint = new HashMap<>();
        Map<Point, Edge> nextEdge = new HashMap<>();

        for (int y = 0; y < h - 1; y++) {
            double up = this.at(0, y);
            double lp = this.at(0, y + 1);

            for (int x = 0; x < w - 1; x++) {
                double ul = up;
                double ur = this.at(x + 1, y);
                double ll = lp;
                double lr = this.at(x + 1, y + 1);

                up = ur;
                lp = lr;

                int squareCase = 0;
                if (ul > z) {
                    squareCase |= 1;
                }
                if (ur > z) {
                    squareCase |= 2;
                }
                if (ll > z) {
                    squareCase |= 4;
                }
                if (lr > z) {
                    squareCase |= 8;
                }

                if (squareCase == 0 || squareCase == 15) {
                    continue;
                }

                double fx = x;
                double fy = y;

                Point t = new Point(fx + fraction(ul, ur, z), fy);
                Point b = new Point(fx + fraction(ll, lr, z), fy + 1);
                Point l = new Point(fx, fy + fraction(ul, ll, z));
                Point r = new Point(fx + 1, fy + fraction(ur, lr, z));

                Edge te = new Edge(x, y, x + 1, y, y == 0);
                Edge be = new Edge(x, y + 1, x + 1, y + 1, y + 2 == h);
                Edge le = new Edge(x, y, x, y + 1, x == 0);
                Edge re = new Edge(x + 1, y, x + 1, y + 1, x + 2 == w);

                final boolean connectHigh = false;

                switch (squareCase) {
                    case 1:
                        edgePoint.put(te, t);
                        nextEdge.put(t, le);
                        edgePoint.put(le, l);
                        break;
                    case 2:
                        edgePoint.put(re, r);
                        nextEdge.put(r, te);
                        edgePoint.put(te, t);
                        break;
                    case 3:
                        edgePoint.put(re, r);
                        nextEdge.put(r, le);
                        edgePoint.put(le, l);
                        break;
                    case 4:
                        edgePoint.put(le, l);
                        nextEdge.put(l, be);
                        edgePoint.put(be, b);
                        break;
                    case 5:
                        edgePoint.put(te, t);
                        nextEdge.put(t, be);
                        edgePoint.put(be, b);
                        break;
                    case 6:
                        if (connectHigh) {
                            edgePoint.put(le, l);
                            nextEdge.put(l, te);
                            edgePoint.put(te, t);

                            edgePoint.put(re, r);
                            nextEdge.put(r, be);
                            edgePoint.put(be, b);
                        } else {
                            edgePoint.put(re, r);
                            nextEdge.put(r, te);
                            edgePoint.put(te, t);

                            edgePoint.put(le, l);
                            nextEdge.put(l, be);
                            edgePoint.put(be, b);
                        }
                        break;
                    case 7:
                        edgePoint.put(re, r);
                        nextEdge.put(r, be);
                        edgePoint.put(be, b);
                        break;
                    case 8:
                        edgePoint.put(be, b);
                        nextEdge.put(b, re);
                        edgePoint.put(re, r);
                        break;
                    case 9:
                        if (connectHigh) {
                            edgePoint.put(te, t);
                            nextEdge.put(t, re);
                            edgePoint.put(re, r);

                            edgePoint.put(be, b);
                            nextEdge.put(b, le);
                            edgePoint.put(le, l);
                        } else {
                            edgePoint.put(te, t);
                            nextEdge.put(t, le);
                            edgePoint.put(le, l);

                            edgePoint.put(be, b);
                            nextEdge.put(b, re);
                            edgePoint.put(re, r);
                        }
                        break;
                    case 10:
                        edgePoint.put(be, b);
                        nextEdge.put(b, te);
                        edgePoint.put(te, t);
                        break;
                    case 11:
                        edgePoint.put(be, b);
                        nextEdge.put(b, le);
                        edgePoint.put(le, l);
                        break;
                    case 12:
                        edgePoint.put(le, l);
                        nextEdge.put(l, re);
                        edgePoint.put(re, r);
                        break;
                    case 13:
                        edgePoint.put(te, t);
                        nextEdge.put(t, re);
                        edgePoint.put(re, r);
                        break;
                    case 14:
                        edgePoint.put(le, l);
                        nextEdge.put(l, te);
                        edgePoint.put(te, t);
                        break;
                }
            }
        }

        // pick out all boundary edgePoints
        Map<Edge, Point> boundaryEdgePoint = new HashMap<>();
        for (Map.Entry<Edge, Point> entry : edgePoint.entrySet()) {
            if (entry.getKey().boundary) {
                boundaryEdgePoint.put(entry.getKey(), entry.getValue());
            }
        }

        List<List<double[]>> contours = new ArrayList<>();

        while (edgePoint.size() > 0) {
            List<double[]> contour = new ArrayList<>();

            // find an unused edge; prefer starting at a boundary
            Edge e = null;
            if (boundaryEdgePoint.size() > 0) {
                for (Edge tmp : boundaryEdgePoint.keySet()) {
                    e = tmp;

                    if (edgePoint.get(e) != null && nextEdge.get(edgePoint.get(e)) != null) {
                        break;
                    }
                }
            } else {
                for (Edge tmp : edgePoint.keySet()) {
                    e = tmp;

                    if (edgePoint.get(e) != null && nextEdge.get(edgePoint.get(e)) != null) {
                        break;
                    }
                }
            }

            Edge e0 = e;

            // add the first point
            // (this allows closed paths to start & end at the same point)
            Point p = edgePoint.get(e);

            contour.add(new double[]{p.x, p.y});
            e = nextEdge.get(p);

            // follow points until none remain
            while (true) {
                p = edgePoint.get(e);
                if (p == null) {
                    break;
                }

                contour.add(new double[]{p.x, p.y});

                edgePoint.remove(e);
                if (e.boundary) {
                    boundaryEdgePoint.remove(e);
                }

                e = nextEdge.get(p);
            }

            // make sure the first one gets deleted in case of open paths
            edgePoint.remove(e0);
            if (e0.boundary) {
                boundaryEdgePoint.remove(e0);
            }

            // add the contour
            contours.add(contour);
        }

        return contours;
    }

    private double fraction(double z0, double z1, double z) {
        double eps = 1e-9;
        double f = 0d;

        if (z0 == CLOSED) {
            f = 0;
        } else if (z1 == CLOSED) {
            f = 1;
        } else if (z0 != z1) {
            f = (z - z0) / (z1 - z0);
        }

        f = Math.max(f, eps);
        f = Math.min(f, 1 - eps);

        return f;
    }

    private double at(int x, int y) {
        return this.grid[y * this.w + x];
    }

    private static class Point {
        double x;
        double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Point point = (Point) o;

            if (Double.compare(point.x, x) != 0) return false;
            return Double.compare(point.y, y) == 0;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            temp = Double.doubleToLongBits(x);
            result = (int) (temp ^ (temp >>> 32));
            temp = Double.doubleToLongBits(y);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    private static class Edge {
        int x0;
        int y0;
        int x1;
        int y1;
        boolean boundary;

        public Edge(int x0, int y0, int x1, int y1, boolean boundary) {
            this.x0 = x0;
            this.y0 = y0;
            this.x1 = x1;
            this.y1 = y1;
            this.boundary = boundary;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Edge edge = (Edge) o;

            if (x0 != edge.x0) return false;
            if (y0 != edge.y0) return false;
            if (x1 != edge.x1) return false;
            if (y1 != edge.y1) return false;
            return boundary == edge.boundary;
        }

        @Override
        public int hashCode() {
            int result = x0;
            result = 31 * result + y0;
            result = 31 * result + x1;
            result = 31 * result + y1;
            result = 31 * result + (boundary ? 1 : 0);
            return result;
        }
    }
}
