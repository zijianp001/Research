public class LineSearch {

    /**
     * Given a x, this method will compute the gradient of f particularly for this function
     */
    public static double[] ComGrad(double[] x) {
        double[] result = new double[2];
        result[0] = (-400) * x[0] * (x[1] - Math.pow(x[0], 2)) - 2 + 2 * x[0];
        result[1] = 200 * x[1] - 200 * (Math.pow(x[0], 2));
        return result;
    }

    /**
     * Notice this method is used to compute the next x in the iteration
     */
    public static double[] x_k_1(double[] x_k, double alpha_k, double[] p_k) {
        double[] result = new double[2];
        result[0] = x_k[0] + alpha_k * p_k[0];
        result[1] = x_k[1] + alpha_k * p_k[1];
        return result;
    }

    /**
     * Given a x, this method will compute the value of f(x) particularly for this function
     */
    public static double FunctionValue(double[] x) {
        double result = 100 * (Math.pow(x[1] - Math.pow(x[0], 2), 2)) + Math.pow(1 - x[0], 2);
        return result;
    }

    /**
     * This is used for steepest descent, and ths dot product of gradient and p_k is actually
     * the norm of gradient times cos.theta, and cos.theta is -1 in this case, this can be used
     * for any function
     */
    public static double ComputeProduct1(double[] f) {
        double result = 0.0;
        for (int i = 0; i < f.length; i++) {
            result += Math.pow(f[i], 2);
        }
        return -Math.sqrt(result);
    }


    public static double[] ComputeDirection1(double[] x) {
        double[] result = ComGrad(x);
        double[] toReturn = new double[2];
        toReturn[0] = result[0] / ComputeProduct1(result);
        toReturn[1] = result[1] / ComputeProduct1(result);
        return toReturn;
    }

    /**
     * Notice I haven't find the way to program the general hessian, i.e for x in n dimension
     */
    public static double[][] ComputeHessian(double[] x) {
        double result[][] = new double[2][2];
        result[0][0] = -400 * x[1] + 1200 * Math.pow(x[0], 2) + 2;
        result[0][1] = -400 * x[0];
        result[1][0] = -400 * x[0];
        result[1][1] = 200;
        return result;
    }

    public static double[] ComputeDirection2(double[] x) {
        double[] result = new double[2];
        double[][] hessian = ComputeHessian(x);
        double[][] copy = hessian;
        double[][] inverse = new double[2][2];
        double coefficient = 1.0 / (copy[0][0] * copy[1][1] - copy[0][1] * copy[1][0]);
        inverse[0][0] = coefficient * copy[1][1];
        inverse[1][0] = -coefficient * copy[1][0];
        inverse[0][1] = -coefficient * copy[0][1];
        inverse[1][1] = coefficient * copy[0][0];
        double[] gradient = ComGrad(x);
        result[0] = -(inverse[0][0] * gradient[0] + inverse[0][1] * gradient[1]);
        result[1] = -(inverse[1][0] * gradient[0] + inverse[1][1] * gradient[1]);
        return result;
    }

    /**
     * This can be used for any function,because gradient and direction can both be represented
     * in one dimension array, i.e a matrix with only one row/column
     */
    public static double ComputeProduct2(double[] x) {
        double[] gradient = ComGrad(x);
        double[] direction = ComputeDirection2(x);
        double result = 0.0;
        for (int i = 0; i < gradient.length; i++) {
            result += gradient[i] * direction[i];
        }
        return result;
    }

    public static void main(String args[]) {
        double alpha = 1.0;
        double a0 = 0.0;
        double[] x = {1.2, 1.2};
        double c = 0.8;
        double ro = 0.5;
        double[] p_k = ComputeDirection1(x);
        double k = ComputeProduct1(ComGrad(x));
        int i = 0;
        int j=0;
        while(FunctionValue(x_k_1(x,alpha,p_k))>(FunctionValue(x)+c*alpha*ComputeProduct1(ComGrad(x)))) {
            alpha=ro*alpha;
            j++;
        }
        a0=alpha;
        while (i < 30000000) {
            x = x_k_1(x, a0, p_k);
            p_k = ComputeDirection1(x);
            i++;
            if(Math.abs(x[0]-1)<0.05 && Math.abs(x[1]-1)<0.05) {
                break;
            }
        }
            System.out.println(ComGrad(x)[0] + " " + ComGrad(x)[1]);
            System.out.println(x[0] + " " + x[1] + " " + i);

        alpha=1.0;
        a0=1.0;
        x[0]=1.2; x[1]=1.2;
        c=0.8;
        ro=0.9;
        p_k=ComputeDirection2(x);
        j=0;
        i=0;
        while(FunctionValue(x_k_1(x,alpha,p_k))>(FunctionValue(x)+c*alpha*ComputeProduct2(ComGrad(x)))) {
            alpha=ro*alpha;
            j++;
        }
        while (i < 10000000) {
            x = x_k_1(x, a0, p_k);
            p_k = ComputeDirection2(x);
            i++;
            if(x[0]==1.0 && x[1]==1.0) {
                break;
            }
        }
        System.out.println(ComGrad(x)[0] + " " + ComGrad(x)[1]);
        System.out.println(x[0] + " " + x[1] + " " + i);



        alpha=1.0;
        a0=1.0;
        x[0]=1.2; x[1]=1.2;
        c=0.8;
        ro=0.9;
        i=0;
        double bk[][]=LineSearch.ComputeHessian(x);
        p_k=Quasi_Newton.ComputeDirection3(x,bk);
        j=0;
        while(FunctionValue(x_k_1(x,alpha,p_k))>(FunctionValue(x)+c*alpha*Quasi_Newton.ComputeProduct3(ComGrad(x),bk))) {
            alpha=ro*alpha;
            j++;
        }
        System.out.println(alpha);
        a0=alpha;
        while(i<2) {
            double[] sk=new double[2];
            double[] yk=new double[2];
            yk[0]=LineSearch.ComGrad(LineSearch.x_k_1(x,a0,p_k))[0]-LineSearch.ComGrad(x)[0];
            yk[1]=LineSearch.ComGrad(LineSearch.x_k_1(x,a0,p_k))[1]-LineSearch.ComGrad(x)[1];
            sk[0]=LineSearch.x_k_1(x,a0,p_k)[0]-x[0];
            sk[1]=LineSearch.x_k_1(x,a0,p_k)[1]-x[1];
            x=LineSearch.x_k_1(x,a0,p_k);
            bk=Quasi_Newton.ComputeB_k_1(bk,sk,yk);
            p_k=Quasi_Newton.ComputeDirection3(x,bk);
            i++;
        }
        System.out.println(ComGrad(x)[0] + " " + ComGrad(x)[1]);
        System.out.println(x[0] + " " + x[1] + " " + i);
        }


}
