import Jama.Matrix;
import edu.princeton.cs.introcs.StdDraw;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Author: Jakub Libera
 * Program Takes N from console
 * -u''(x) - u(x) = sin(x)
 * u(0) = 0
 * u'(2) = u(2)
 * x e [0, 2]
 */

/**
 *
 * Program uruchamiamy z 2 argumentami
 *  1 - definiuje nasze N
 *  2 - definiuje dokładność naszego wykresu tzn ilość punktów z których zbudujemy nasz wykres
 * W przypadku nie podania argumentow program wystartuje z N = 1000
 * więcej szegółow w funkcji main
 *
 */

public class Main {

    public static void main(String[] args) {
        List<FunctionInfo> list = getAllFunctions((args.length == 1) ? Integer.parseInt(args[0]) : 1000);
        plotResult(getBMatrix(list).solve(getLMatrix(list)), list, ((args.length == 2) ? Integer.parseInt(args[1]) : 1000));
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
        StdDraw.setXscale(0, 2);
        StdDraw.setYscale(-2.0, +2.0);
        StdDraw.line(0, 2, 0, -2);
        StdDraw.line(0, 0, 2, 0);
        StdDraw.setPenColor(Color.RED);
        for (int i = 0; i <= n; i++)
            StdDraw.point(2.0 * (double) i / (double) n,
                    getYFromFunctions(matrix.getArray(), list, 2.0 * (double) i / (double) n));
        /**
          SPODZIEWANY WYNIK
         **/
        StdDraw.setPenColor(Color.GREEN);
        for (int i = 0; i < n; i++)
            StdDraw.point(2.0 * (double) i / (double) n,
                    (2.0 * (double) i / (double) n * Math.cos(2.0 * (double) i / (double) n) -  1.0581 * Math.sin(2.0 * (double) i / (double) n)) / 2);
    }

    private static double getYFromFunctions(double[][] matrix, List<FunctionInfo> list, double x) {
        double res = 0;
        for (int i = 0; i < list.size(); i++)
            if (list.get(i).getLeftBound() <= x && list.get(i).getRightBound() >= x)
                res += list.get(i).getFunctionResults(x) * matrix[i][0];
        return res;
    }
}


