import java.util.Arrays;

public class Cholesky {
    public static void main(String args[]){
        double[][] a={{-3,-4,-5},{-56,-17,-3},{-1,-1,-1}};
        System.out.println(ComputeFactorization(a));
        System.out.println(ComputeMultiplier(a));
    }
    public static boolean ComputeFactorization(double[][] a){
        double [][] c=new double[a.length][a.length];
        double [][] d=new double[a.length][a.length];
        double [][] l=new double[a.length][a.length];
        for(int i=0;i<a.length;i++) {
            l[i][i]=1;
        }
        for(int j=0;j<a.length;j++) {
            double sum=0.0;
            for(int s=0;s<j;s++) {
                sum+=(d[s][s]*Math.pow(l[j][s],2.0));
            }
            c[j][j]=a[j][j]-sum;
            d[j][j]=c[j][j];
            for(int i=j+1;i<a.length;i++) {
                sum=0.0;
                for(int s=0;s<j;s++) {
                    sum+=(d[s][s]*l[i][s]*l[j][s]);
                }
                c[i][j]=a[i][j]-sum;
                l[i][j]=c[i][j]/d[j][j];
            }
        }
        System.out.println(Arrays.deepToString(l));
        System.out.println(Arrays.deepToString(d));
        for(int i=0;i<a.length;i++){
            if(d[i][i]<=0) {
                return false;
            }
        }
        return true;
    }

    public static double ComputeMultiplier(double[][] a) {
        double min=a[0][0];
        double beta=1;
        double t;
        for(int i=0;i<a.length;i++) {
            if(a[i][i]<=min){
                min=a[i][i];
            }
        }
        if(min>0) {
            t=0;
        }
        else{
            t=-min+beta;
        }
        boolean fact=false;
        while(fact==false){
            double [][] b=a;
            for(int i=0;i<a.length;i++) {
                b[i][i]+=t;
            }
            fact=ComputeFactorization(b);
            System.out.println();
            if(fact==true){
                return t;
            }
            else {
                t=Math.max(2*t,beta);
            }
        }
        return -1;
    }
}
