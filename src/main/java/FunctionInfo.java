class FunctionInfo {
    private Function leftFunction;
    private Function rightFunction;
    private double leftBound;
    private double rightBound;

    FunctionInfo(double leftBound, double rightBound, Function leftFunction, Function rightFunction) {
        this.leftBound = leftBound;
        this.rightBound = rightBound;
        this.leftFunction = leftFunction;
        this.rightFunction = rightFunction;
    }

    double getFunctionResults(double x) {
        return (x > leftBound && x < rightBound)
                ? (x > (leftBound + rightBound) / 2
                ? rightFunction.functionResult(x)
                : leftFunction.functionResult(x))
                : 0;
    }

    double getLeftBound() {
        return leftBound;
    }

    double getRightBound() {
        return rightBound;
    }

    Function getLeftFunction() {
        return leftFunction;
    }

    Function getRightFunction() {
        return rightFunction;
    }

    void setLeftFunction(Function leftFunction) {
        this.leftFunction = leftFunction;
    }

    void setRightFunction(Function rightFunction) {
        this.rightFunction = rightFunction;
    }

}
