import Jama.Matrix;
import edu.princeton.cs.introcs.StdDraw;

import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * Author: Jakub Libera
 * Program Takes N from console
 * -u''(x) - u(x) = sin(x)
 * u(0) = 0
 * u'(2) = u(2)
 * x e [0, 2]
 */

public class Main {

    public static void main(String[] args) {
        List<FunctionInfo> list = getAllFunctions(Integer.parseInt(new Scanner(System.in).nextLine()));
        plotResult(getBMatrix(list).solve(getLMatrix(list)), list, 2000);
    }
    
    private static double gaussianQuadrature(double a, double b, Function function) {
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
            temp[j][0] = gaussianQuadrature(i, i + 2.0 / functions.size(), (double x) -> functions.get(j).getFunctionResults(x) * Math.sin(x)) +
                    gaussianQuadrature(i + 2.0 / functions.size(), i + (double) 4 / functions.size(), (double x) -> functions.get(j).getFunctionResults(x) * Math.sin(x));
            i += 2.0 / functions.size();
            k++;
        }
        return new Matrix(temp);
    }

    private static Matrix getBMatrix(List<FunctionInfo> functions) {
        double[][] temp = new double[functions.size()][functions.size()];
        for (int i = 0; i < functions.size(); i++)
            for (int j = 0; j < functions.size(); j++)
                temp[i][j] = getBFunction(functions.get(i), functions.get(j));
        return new Matrix(temp);
    }

    private static double getBFunction(FunctionInfo functionInfo, FunctionInfo functionInfo1) {
        double res = -functionInfo.getFunctionResults(2.0) * functionInfo1.getFunctionResults(2.0);
        double a = Math.max(functionInfo.getLeftBound(), functionInfo1.getLeftBound());
        double b = Math.min(functionInfo.getRightBound(), functionInfo1.getRightBound());
        if (a >= b)
            return res;
        if (functionInfo.equals(functionInfo1)) {
            res -= gaussianQuadrature(a, (b + a) / 2.0,
                    (double x) -> functionInfo.getLeftFunction().functionResult(x) *
                            functionInfo1.getLeftFunction().functionResult(x));
            res -= gaussianQuadrature((b + a) / 2.0, b,
                    (double x) -> functionInfo.getRightFunction().functionResult(x) *
                            functionInfo1.getRightFunction().functionResult(x));
            res += gaussianQuadrature(a, (b + a) / 2.0,
                    (double x) -> derivativeApproximation(functionInfo.getLeftFunction()).functionResult(x) *
                            derivativeApproximation(functionInfo1.getLeftFunction()).functionResult(x));
            res += gaussianQuadrature((b + a) / 2.0, b,
                    (double x) -> derivativeApproximation(functionInfo.getRightFunction()).functionResult(x) *
                            derivativeApproximation(functionInfo1.getRightFunction()).functionResult(x));
        } else {
            if (functionInfo.getRightBound() > functionInfo1.getRightBound()) {
                res -= gaussianQuadrature(a, b,
                        (double x) -> functionInfo.getLeftFunction().functionResult(x) *
                                functionInfo1.getRightFunction().functionResult(x));
                res += gaussianQuadrature(a, b,
                        (double x) -> derivativeApproximation(functionInfo.getLeftFunction()).functionResult(x) *
                                derivativeApproximation(functionInfo1.getRightFunction()).functionResult(x));
            } else {
                res -= gaussianQuadrature(a, b,
                        (double x) -> functionInfo.getRightFunction().functionResult(x) *
                                functionInfo1.getLeftFunction().functionResult(x));
                res += gaussianQuadrature(a, b,
                        (double x) -> derivativeApproximation(functionInfo.getRightFunction()).functionResult(x) *
                                derivativeApproximation(functionInfo1.getLeftFunction()).functionResult(x));
            }
        }
        return res;
    }

    private static Function derivativeApproximation(Function function) {
        return (double x) -> (function.functionResult(x + 0.001) - function.functionResult(x - 0.001)) / 0.001 / 2;
    }

    private static void plotResult(Matrix matrix, List<FunctionInfo> list, int n) {
        double[] x = new double[n + 1];
        double[] y = new double[n + 1];
        for (int i = 0; i <= n; i++) {
            x[i] = 2.0 * (double) i / (double) n;
            y[i] = getYFromFunctions(matrix.getArray(), list, x[i]);
        }
        StdDraw.setXscale(0, 2);
        StdDraw.setYscale(-2.0, +2.0);
        for (int i = 0; i < n; i++)
            StdDraw.point(x[i], y[i]);
        /*
          SPODZIEWANY WYNIK
         */
//        for (int i = 0; i < n; i++)
//            StdDraw.point(x[i], (x[i] * Math.cos(x[i]) -  1.0581 * Math.sin(x[i])) / 2);
    }

    private static double getYFromFunctions(double[][] matrix, List<FunctionInfo> list, double x) {
        double res = 0;
        for (int i = 0; i < list.size(); i++)
            if (list.get(i).getLeftBound() <= x && list.get(i).getRightBound() >= x)
                res += list.get(i).getFunctionResults(x) * matrix[i][0];
        return res;
    }
}
