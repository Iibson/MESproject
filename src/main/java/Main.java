import Jama.Matrix;
import edu.princeton.cs.introcs.StdDraw;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {

    public static void main(String[] args) {
//        List<FunctionInfo> list = getAllFunctions(Integer.parseInt(args[0]));
        List<FunctionInfo> list = getAllFunctions(3);
        Matrix lMatrix = getLMatrix(list);
        Matrix res = getBMatrix(list).solve(lMatrix);
        double[][] res1 = res.getArray();
        for (double[] doubles : res1)
            for (double aDouble : doubles)
                System.out.print(" " + aDouble);
        plotResult();
    }

    private static double integrate(double a, double b, Function function) {
        return (b - a) / 2 *
                (function.functionResult((b - a) / 2 / Math.sqrt(3) + (a + b) / 2) +
                        function.functionResult((b - a) / 2 / -Math.sqrt(3) + (a + b) / 2));
    }

    private static FunctionInfo getFunctionForArea(double a, double b) {
        return new FunctionInfo(a, b, (double x) -> 2 / (b - a) * x - 2 * a / (b - a), (double x) -> 2 / (a - b) * x - 2 * b / (a - b));
    }

    private static List<FunctionInfo> getAllFunctions(int n) {
        return new ArrayList<>() {{
            double i = 0;
            while (i < 2.0) {
                add(getFunctionForArea(i, i + (4.0 / n)));
                i += 2.0 / n;
            }
            get(size() - 1).setRightFunction((double x) -> 0);
        }};
    }

    private static Matrix getLMatrix(List<FunctionInfo> functions) {
        double[][] temp = new double[functions.size()][1];
        double i = 0.0;
        int k = 0;
        while (i < 2.0) {
            int j = k;
            temp[j][0] = integrate(i, i + 2.0 / functions.size(), (double x) -> functions.get(j).getFunctionResults(x) * Math.sin(x)) +
                    integrate(i + 2.0 / functions.size(), i + (double) 4 / functions.size(), (double x) -> functions.get(j).getFunctionResults(x) * Math.sin(x));
            i += 2.0 / functions.size();
            k++;
        }
        return new Matrix(temp);
    }

    private static Matrix getBMatrix(List<FunctionInfo> functions) {
        double[][] temp = new double[functions.size()][functions.size()];
        for (int i = 0; i < functions.size(); i++) {
            for (int j = 0; j < functions.size(); j++) {
                temp[i][j] = getBFunction(functions.get(i), functions.get(j));
            }
        }
        return new Matrix(temp);
    }

    private static double getBFunction(FunctionInfo functionInfo, FunctionInfo functionInfo1) {
        double res = -functionInfo.getFunctionResults(2.0) * functionInfo1.getFunctionResults(2.0);
        double a = Math.max(functionInfo.getLeftBound(), functionInfo1.getLeftBound());
        double b = Math.min(functionInfo.getRightBound(), functionInfo1.getRightBound());
        if (a >= b)
            return res;
        if (functionInfo.equals(functionInfo1)) {
            res -= integrate(a, (b + a) / 2.0,
                    (double x) -> functionInfo.getLeftFunction().functionResult(x) *
                            functionInfo1.getLeftFunction().functionResult(x));
            res -= integrate((b + a) / 2.0, b,
                    (double x) -> functionInfo.getRightFunction().functionResult(x) *
                            functionInfo1.getRightFunction().functionResult(x));
            res += integrate(a, (b + a) / 2.0,
                    (double x) -> derivativeApproximation(functionInfo.getLeftFunction()).functionResult(x) *
                            derivativeApproximation(functionInfo1.getLeftFunction()).functionResult(x));
            res += integrate((b + a) / 2.0, b,
                    (double x) -> derivativeApproximation(functionInfo.getRightFunction()).functionResult(x) *
                            derivativeApproximation(functionInfo1.getRightFunction()).functionResult(x));
        } else {
            if (functionInfo.getRightBound() > functionInfo1.getRightBound()) {
                res -= integrate(a, b,
                        (double x) -> functionInfo.getLeftFunction().functionResult(x) *
                                functionInfo1.getRightFunction().functionResult(x));
                res += integrate(a, b,
                        (double x) -> derivativeApproximation(functionInfo.getLeftFunction()).functionResult(x) *
                                derivativeApproximation(functionInfo1.getRightFunction()).functionResult(x));
            } else {
                res -= integrate(a, b,
                        (double x) -> functionInfo.getRightFunction().functionResult(x) *
                                functionInfo1.getLeftFunction().functionResult(x));
                res += integrate(a, b,
                        (double x) -> derivativeApproximation(functionInfo.getRightFunction()).functionResult(x) *
                                derivativeApproximation(functionInfo1.getLeftFunction()).functionResult(x));
            }
        }
        return res;
    }

    private static Function derivativeApproximation(Function function) {
        return (double x) -> (function.functionResult(x + 0.001) - function.functionResult(x - 0.001)) / 0.001 / 2;
    }

    private static void plotResult() {
        // number of line segments to plot
        int n = 500;

        // the function y = sin(4x) + sin(20x), sampled at n+1 points
        // between x = 0 and x = pi
        double[] x = new double[n + 1];
        double[] y = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            x[i] = Math.PI * i / n;
            y[i] = Math.sin(4 * x[i]) + Math.sin(20 * x[i]);
        }

        // rescale the coordinate system
        StdDraw.setXscale(0, 2);
        StdDraw.setYscale(-2.0, +2.0);

        // plot the approximation to the function
        for (int i = 0; i < n; i++) {
            StdDraw.line(x[i], y[i], x[i + 1], y[i + 1]);
        }
    }
}
