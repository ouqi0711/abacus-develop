#ifndef MATH_BSPLINE_H
#define MATH_BSPLINE_H

namespace ModuleBase
{


//
//DESCRIPTION:
//   A class to treat Cardinal B-spline interpolation.
//MATH:
//   Only uniform nodes are considered: xm-x[m-1]=Dx(>= 0) for control node: X={x0,x1,...,xm};
//   Any function p(x) can be written by
//   p(x)=\sum_i{ci*M_ik(x)} (k->infinity), 
//   where M_ik is the i-th k-order Cardinal B-spline base function 
//   and ci is undetermined coefficient.
//      M_i0 =  H(x-xi)-H(x-x[i+1]), H(x): step function
//                  x-xi                    x[i+k+1]-x
//      M_ik(x)=  ---------*M_i(k-1)(x)+ ----------------*M_[i+1][k-1](x)  ( xi <= x <= x[i+1] )
//                x[i+k]-xi               x[i+k+1]-x[i+1]
//   For uniform nodes: M_[i+1]k(x+Dx)=M_ik(x)
//   If we define Bk(n) stores M_ik(x+n*Dx) for x in (xi,xi+Dx):
//               x+n*Dx-xi               xi+(k-n+1)*Dx-x
//      Bk[n] = -----------*B(k-1)[n] + -----------------*B(k-1)[n-1]
//                 k*Dx                        k*Dx
//USAGE：
//   ModuleBase::Bspline bp;
//   bp.init(10,2,0.5); //Dx = 2, xi = 0.5
//   bp.getbslpine(1.1); //x = 1.1
//   cout<<bp.bspline_ele(3)<<endl; //print M_ik(xi+3*Dx+x): M_i[10](7.6) 
//AUTHOR: qianrui on 2021-9-24
//
class Bspline
{
    private:
        int norder; // the order of bezier base; norder >= 0
        double Dx; //Dx: the interval of control node
        double xi; // xi: the starting point
        double *bezier; //bezier[n] = Bk[n]

    public:
        Bspline();
        ~Bspline();

        //Init norder, Dx, xi
        void init(const int norderin, const double Dxin = 1, const double xiin = 0);

        //delete[] bezier
        void cleanp();

        //Get the result of i-th bezier base functions for different input x+xi+n*Dx.
        //x should be in [0,Dx]
        //n-th result is stored in bezier[n];
        void getbslpine(double x);
        
        //get the element of bezier
        double bezier_ele(int n);
};
}
#endif
