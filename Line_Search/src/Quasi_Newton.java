import java.util.Arrays;

public class Quasi_Newton {

    /** This method will demonstrate the BFGS formula
     */
    public static double[][] ComputeB_k_1(double[][] bk, double[] sk, double[] yk) {

        double[][] result=new double[2][2];
        /**First we compute the second term, which is (B_k*s_k*s_k^T*B_k)/(s_k^T*B_k*s_k)
         * We may compute s_k*s_k^T first
         */
        double[][] sk_skt=new double[2][2];
        sk_skt[0][0]=sk[0]*sk[0];
        sk_skt[1][0]=sk[1]*sk[0];
        sk_skt[0][1]=sk[0]*sk[1];
        sk_skt[1][1]=sk[1]*sk[1];
        double[][] numerator=SquareMul(SquareMul(bk,sk_skt),bk);
        double[] part=new double[2];
        part[0]=sk[0]*bk[0][0]+sk[1]*bk[1][0];
        part[1]=sk[0]*bk[0][1]+sk[1]*bk[1][1];
        double denominator=part[0]*sk[0]+part[1]*sk[1];
        double[][] term2=new double[2][2];
        term2[0][0]=numerator[0][0]/denominator;
        term2[0][1]=numerator[0][1]/denominator;
        term2[1][0]=numerator[1][0]/denominator;
        term2[1][1]=numerator[1][1]/denominator;

        /** Now we need to compute the third term, which is (y_k*y_k^T)/(y_k^T*s_k)
         */
        double[][] den=new double[2][2];
        den[0][0]=yk[0]*yk[0];
        den[0][1]=yk[1]*yk[0];
        den[1][0]=yk[0]*yk[1];
        den[1][1]=yk[1]*yk[1];
        double num=yk[0]*sk[0]+yk[1]*sk[1];
        double[][] term3=new double[2][2];
        term3[0][0]=den[0][0]/num;
        term3[0][1]=den[0][1]/num;
        term3[1][0]=den[1][0]/num;
        term3[1][1]=den[1][1]/num;

        result[0][0]=bk[0][0]-term2[0][0]+term3[0][0];
        result[0][1]=bk[0][1]-term2[0][1]+term3[0][1];
        result[1][0]=bk[1][0]-term2[1][0]+term3[1][0];
        result[1][1]=bk[1][1]-term2[1][1]+term3[1][1];
        return result;
    }

    public static double[][] SquareMul(double[][] m1, double[][] m2){
        double[][] result=new double[2][2];
        result[0][0]=m1[0][0]*m2[0][0]+m1[0][1]*m2[1][0];
        result[1][0]=m1[1][0]*m2[0][0]+m1[1][1]*m2[1][0];
        result[0][1]=m1[0][0]*m2[0][1]+m1[0][1]*m2[1][1];
        result[1][1]=m1[1][0]*m2[0][1]+m1[1][1]*m2[1][1];
        return result;
    }

    public static double[] ComputeDirection3(double[] x,double[][] bk) {
        double[] result=new double[2];
        double[][] copy=bk;
        double[][] inverse=new double[2][2];
        double coefficient=1.0/(copy[0][0]*copy[1][1]-copy[0][1]*copy[1][0]);
        inverse[0][0]=coefficient*copy[1][1];
        inverse[1][0]=-coefficient*copy[1][0];
        inverse[0][1]=-coefficient*copy[0][1];
        inverse[1][1]=coefficient*copy[0][0];
        double[] gradient=LineSearch.ComGrad(x);
        result[0]=-(inverse[0][0]*gradient[0]+inverse[0][1]*gradient[1]);
        result[1]=-(inverse[1][0]*gradient[0]+inverse[1][1]*gradient[1]);
        return result;
    }
    public static double ComputeProduct3(double[] x, double[][] bk) {
        double[] gradient=LineSearch.ComGrad(x);
        double[] direction=ComputeDirection3(x,bk);
        double result=0.0;
        for(int i=0;i<gradient.length;i++){
            result+=gradient[i] * direction[i];
        }
        return result;

    }
    public static void main(String args[]) {
        double alpha=1.0;
        double a0=1.0;
        double[] x={1.2,1.2};
        double c=0.8;
        double ro=0.5;
        double bk[][]=LineSearch.ComputeHessian(x);
        double[] p_k=ComputeDirection3(x,bk);
        double k=ComputeProduct3(LineSearch.ComGrad(x),bk);
        int i=0;
        a0=1.443614920627236E-31;
        while(i<20000) {
            double[] sk=new double[2];
            double[] yk=new double[2];
            yk[0]=LineSearch.ComGrad(LineSearch.x_k_1(x,a0,p_k))[0]-LineSearch.ComGrad(x)[0];
            yk[1]=LineSearch.ComGrad(LineSearch.x_k_1(x,a0,p_k))[1]-LineSearch.ComGrad(x)[1];
            sk[0]=LineSearch.x_k_1(x,a0,p_k)[0]-x[0];
            sk[1]=LineSearch.x_k_1(x,a0,p_k)[1]-x[1];
            x=LineSearch.x_k_1(x,a0,p_k);
            alpha=ro*alpha;
            //WOLFE CONDITIONS OR GOLDSTEIN CONDITIONS HERE
            bk=ComputeB_k_1(bk,sk,yk);
            p_k=ComputeDirection3(x,bk);
            i++;
        }
        System.out.println(LineSearch.ComGrad(x)[0] +" "+LineSearch.ComGrad(x)[1]);
        System.out.println(x[0]+" "+x[1]+" "+i);
    }

}
