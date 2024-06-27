#include<iostream>
#include<algorithm> //for using max function
#include<cstdlib> //for using atof
#include<cmath> //for using fabs
#include<fstream> //for using files
#include<limits> //for using numeric_limits
#include<vector> //for using vectors
#include<sstream>
using namespace std;

int M, N, n, k, l;
double e1, s1, I1, e2, s2, I2, dx, dz, l_B, phi_D, phi_P;
double A, B, C, D;
vector<double> phi_low, phi_up;

// u runs from -M-1 to M and v runs from 0 to N-1. Thus r runs from -M to M and s runs from 0 to N

int I(int r, int s)
{
  return (r+M)*(N+1)+s;
}

//define function Pow for evaluating x^a
double Pow(double x, int a)
{
  if(a<0) return 0.;
  if(a==0) return 1.;
  if(x==0.) return 0.;
  double factor=1.;
  if((a%2!=0)&&(x<0))
    factor=-1.;
  return pow(fabs(x),((double)a))*factor;
}

//define the prefactor Q
double Q(int k)
{
  double p=-lgamma((double)(2*k+2));
  return exp(p);
}

//define the prefactor Q1
double Q1(int i, int k)
{
  double p=lgamma((double)(i+1))-lgamma((double)(k+1))-lgamma((double)(i-k+1));
  return exp(p);
}

//define the function Integral1: x^i(u+vx)^j
double Integral1(double u, double v, int i, int j)
{
  if(j>0)
    return Pow(u+v,j)/(double)(i+1)-v*((double)j/(double)(i+1))*Integral1(u,v,i+1,j-1);
  else
    return 1/(double)(i+1);
}

//define the function Integral2
double Integral2(double u, double v, int m, int n)
{
  if(m>0)
    return (1/(double)(n+1))*(Pow(1+u,m)*Pow(1+v,n+1)-Pow(u,m)*Pow(v,n+1)-m*Integral2(u,v,m-1,n+1));
  else
    return (Pow(1+v,n+1)-Pow(v,n+1))/(double)(n+1);
}

//define the function Integral
double Integral(double a, double b, double c, double d, int m, int n)
{
    if((fabs(b)<=numeric_limits<double>::epsilon())&&(fabs(d)<=numeric_limits<double>::epsilon()))
      return Pow(a,m)*Pow(c,n);
    else if((fabs(b)<=numeric_limits<double>::epsilon())&&(fabs(d)>numeric_limits<double>::epsilon()))
    {
      if(fabs(c)>numeric_limits<double>::epsilon())
      {
	double s1=0.0;
        for(int j=0; j<=n; ++j)
        {
	  s1+=Pow(1+(double)d/(double)c,j);
        }
        return Pow(a,m)*Pow(c,n)*s1/((double)(n+1));
        }
      else
      {
	return Pow(a,m)*Pow(d,n)/(double)(n+1);
      }
    }
    else if((fabs(b)>numeric_limits<double>::epsilon())&&(fabs(d)<=numeric_limits<double>::epsilon()))
    {
      if(fabs(a)>numeric_limits<double>::epsilon())
      {
	double s1=0.0;
        for(int j=0; j<=m; ++j)
        {
	  s1+=Pow(1+(double)b/(double)a,j);
        }
        return Pow(a,m)*Pow(c,n)*s1/((double)(m+1));
      }
      else
      {
	return Pow(b,m)*Pow(c,n)/(double)(m+1);
      }
    }
    else if(((fabs(b)>numeric_limits<double>::epsilon())&&(fabs(b)<1.))||((fabs(d)>numeric_limits<double>::epsilon())&&(fabs(d)<1.)))
    {
      if(m<=n)
      {
      double s1=0.;
      for(int k=0; k<=m; ++k)
      {
	s1+=Pow(a,k)*Pow(b,m-k)*Q1(m,k)*Integral1(c,d,m-k,n);
      }
      return s1;
      }
      else
      {
      double s1=0.;
      for(int k=0; k<=n; ++k)
      {
	s1+=Pow(c,k)*Pow(d,n-k)*Q1(n,k)*Integral1(a,b,n-k,m);
      }
      return s1;
      }
    }
    else
    {
      if(m<=n)
	return Pow(b,m)*Pow(d,n)*Integral2(a/b,c/d,m,n);
      else
	return Pow(b,m)*Pow(d,n)*Integral2(c/d,a/b,n,m);
    }
}

//define the function f1
double f1(int r, int s, int k, int l, double phi_b, const vector<double> &phi)
{
  double phi_r_s_mod, phi_rp1_s_mod, phi_r_sm1_mod, phi_rp1_sm1_mod;  
  if(s==0)
    return 0.;
    phi_r_s_mod=phi[I(r,s)];
    phi_r_sm1_mod=phi[I(r,s-1)];
    if(r==M)
    {
      phi_rp1_s_mod=phi_up[s];
      phi_rp1_sm1_mod=phi_up[s-1];
    }
    else
    {
      phi_rp1_s_mod=phi[I(r+1,s)];
      phi_rp1_sm1_mod=phi[I(r+1,s-1)];
    }
  if(2*k-l==0)
    return 0.;
  else if(fabs(phi_rp1_s_mod-phi_r_s_mod)<=numeric_limits<double>::epsilon())
  {
    return (2*k-l)*Pow(phi_r_s_mod-phi_b,2*k-l-1)*Integral(1.,-1.,phi_r_sm1_mod-phi_b,phi_rp1_sm1_mod-phi_r_sm1_mod,1,l);
  }
  else
  {
    return (2*k-l)*((phi_rp1_s_mod-phi_b)/(phi_rp1_s_mod-phi_r_s_mod))*Integral(phi_r_s_mod-phi_b,phi_rp1_s_mod-phi_r_s_mod,phi_r_sm1_mod-phi_b,phi_rp1_sm1_mod-phi_r_sm1_mod,2*k-l-1,l)
    -(2*k-l)*(1/(phi_rp1_s_mod-phi_r_s_mod))*Integral(phi_r_s_mod-phi_b,phi_rp1_s_mod-phi_r_s_mod,phi_r_sm1_mod-phi_b,phi_rp1_sm1_mod-phi_r_sm1_mod,2*k-l,l);
  }
}

//define the function f2
double f2(int r, int s, int k, int l, double phi_b, const vector<double> &phi)
{
  double phi_r_s_mod, phi_rm1_s_mod, phi_r_sm1_mod, phi_rm1_sm1_mod;  
  if(s==0)
    return 0.;
    phi_r_s_mod=phi[I(r,s)];
    phi_r_sm1_mod=phi[I(r,s-1)];
    if(r==-M)
    {
      phi_rm1_s_mod=phi_low[s];
      phi_rm1_sm1_mod=phi_low[s-1];
    }
    else
    {
      phi_rm1_s_mod=phi[I(r-1,s)];
      phi_rm1_sm1_mod=phi[I(r-1,s-1)];
    }
  if(2*k-l==0)
    return 0.;
  else if(fabs(phi_r_s_mod-phi_rm1_s_mod)<=numeric_limits<double>::epsilon())
  {
    return (2*k-l)*Pow(phi_rm1_s_mod-phi_b,2*k-l-1)*Integral(0.,1.,phi_rm1_sm1_mod-phi_b,phi_r_sm1_mod-phi_rm1_sm1_mod,1,l);
  }
  else
  {
    return (2*k-l)*(1/(phi_r_s_mod-phi_rm1_s_mod))*Integral(phi_rm1_s_mod-phi_b,phi_r_s_mod-phi_rm1_s_mod,phi_rm1_sm1_mod-phi_b,phi_r_sm1_mod-phi_rm1_sm1_mod,2*k-l,l)
    -(2*k-l)*((phi_rm1_s_mod-phi_b)/(phi_r_s_mod-phi_rm1_s_mod))*Integral(phi_rm1_s_mod-phi_b,phi_r_s_mod-phi_rm1_s_mod,phi_rm1_sm1_mod-phi_b,phi_r_sm1_mod-phi_rm1_sm1_mod,2*k-l-1,l);
  }
}

//define the function f3
double f3(int r, int s, int k, int l, double phi_b, const vector<double> &phi)
{
  double phi_r_s_mod, phi_rp1_s_mod, phi_r_sp1_mod, phi_rp1_sp1_mod;  
  if(s==N)
    return 0.;
    phi_r_s_mod=phi[I(r,s)];
    phi_r_sp1_mod=phi[I(r,s+1)];
    if(r==M)
    {
      phi_rp1_s_mod=phi_up[s];
      phi_rp1_sp1_mod=phi_up[s+1];
    }
    else
    {
      phi_rp1_s_mod=phi[I(r+1,s)];
      phi_rp1_sp1_mod=phi[I(r+1,s+1)];
    }
  if(l==0)
    return 0.;
  else if(fabs(phi_rp1_s_mod-phi_r_s_mod)<=numeric_limits<double>::epsilon())
  {
    return l*Pow(phi_r_s_mod-phi_b,l-1)*Integral(1.,-1.,phi_r_sp1_mod-phi_b,phi_rp1_sp1_mod-phi_r_sp1_mod,1,2*k-l);
  }
  else
  {
    return l*((phi_rp1_s_mod-phi_b)/(phi_rp1_s_mod-phi_r_s_mod))*Integral(phi_r_sp1_mod-phi_b,phi_rp1_sp1_mod-phi_r_sp1_mod,phi_r_s_mod-phi_b,phi_rp1_s_mod-phi_r_s_mod,2*k-l,l-1)
    -l*(1/(phi_rp1_s_mod-phi_r_s_mod))*Integral(phi_r_sp1_mod-phi_b,phi_rp1_sp1_mod-phi_r_sp1_mod,phi_r_s_mod-phi_b,phi_rp1_s_mod-phi_r_s_mod,2*k-l,l);
  }
}

//define the function f4
double f4(int r, int s, int k, int l, double phi_b, const vector<double> &phi)
{
  double phi_r_s_mod, phi_rm1_s_mod, phi_r_sp1_mod, phi_rm1_sp1_mod;  
  if(s==N)
    return 0.;
    phi_r_s_mod=phi[I(r,s)];
    phi_r_sp1_mod=phi[I(r,s+1)];
    if(r==-M)
    {
      phi_rm1_s_mod=phi_low[s];
      phi_rm1_sp1_mod=phi_low[s+1];
    }
    else
    {
      phi_rm1_s_mod=phi[I(r-1,s)];
      phi_rm1_sp1_mod=phi[I(r-1,s+1)];
    }
  if(l==0)
    return 0.;
  else if(fabs(phi_r_s_mod-phi_rm1_s_mod)<=numeric_limits<double>::epsilon())
  {
    return l*Pow(phi_rm1_s_mod-phi_b,l-1)*Integral(0.,1.,phi_rm1_sp1_mod-phi_b,phi_r_sp1_mod-phi_rm1_sp1_mod,1,2*k-l);
  }
  else
  {
    return l*(1/(phi_r_s_mod-phi_rm1_s_mod))*Integral(phi_rm1_sp1_mod-phi_b,phi_r_sp1_mod-phi_rm1_sp1_mod,phi_rm1_s_mod-phi_b,phi_r_s_mod-phi_rm1_s_mod,2*k-l,l)
    -l*((phi_rm1_s_mod-phi_b)/(phi_r_s_mod-phi_rm1_s_mod))*Integral(phi_rm1_sp1_mod-phi_b,phi_r_sp1_mod-phi_rm1_sp1_mod,phi_rm1_s_mod-phi_b,phi_r_s_mod-phi_rm1_s_mod,2*k-l,l-1);
  }
}

//define the function f5
double f5(int r, int s, const vector<double> &phi)
{
  if((s==0)&&(r<M)&&(r>=-M))
    return A*phi[I(r,s)]+B*phi[I(r+1,s)]+C*phi[I(r,s+1)]+D*phi[I(r+1,s+1)];
  else if((s==N)&&(r<M)&&(r>=-M))
    return A*phi[I(r,s)]+B*phi[I(r+1,s)]+C*phi[I(r,s-1)]+D*phi[I(r+1,s-1)];
  else if((s==0)&&(r==M))
    return A*phi[I(r,s)]+B*phi_up[s]+C*phi[I(r,s+1)]+D*phi_up[s+1];
  else if((s==N)&&(r==M))
    return A*phi[I(r,s)]+B*phi_up[s]+C*phi[I(r,s-1)]+D*phi_up[s-1];
  else if((s>0)&&(s<N)&&(r==M))
    return 2*A*phi[I(r,s)]+2*B*phi_up[s]+C*(phi[I(r,s-1)]+phi[I(r,s+1)])+D*(phi_up[s-1]+phi_up[s+1]);
  else if((s>0)&&(s<N)&&(r<M)&&(r>=-M))
    return 2*A*phi[I(r,s)]+2*B*phi[I(r+1,s)]+C*(phi[I(r,s-1)]+phi[I(r,s+1)])+D*(phi[I(r+1,s-1)]+phi[I(r+1,s+1)]);
  else
  {
    cerr<<"f5 impossible case"<< endl;
    exit(3);
  }
}

//define the function f6
double f6(int r, int s, const vector<double> &phi)
{
  if((s==0)&&(r>-M)&&(r<=M))
    return A*phi[I(r,s)]+B*phi[I(r-1,s)]+C*phi[I(r,s+1)]+D*phi[I(r-1,s+1)];
  else if((s==N)&&(r>-M)&&(r<=M))
    return A*phi[I(r,s)]+B*phi[I(r-1,s)]+C*phi[I(r,s-1)]+D*phi[I(r-1,s-1)];
  else if((s==0)&&(r==-M))
    return A*phi[I(r,s)]+B*phi_low[s]+C*phi[I(r,s+1)]+D*phi_low[s+1];
  else if((s==N)&&(r==-M))
    return A*phi[I(r,s)]+B*phi_low[s]+C*phi[I(r,s-1)]+D*phi_low[s-1];
  else if((s>0)&&(s<N)&&(r==-M))
    return 2*A*phi[I(r,s)]+2*B*phi_low[s]+C*(phi[I(r,s-1)]+phi[I(r,s+1)])+D*(phi_low[s-1]+phi_low[s+1]);
  else if((s>0)&&(s<N)&&(r>-M)&&(r<=M))
    return 2*A*phi[I(r,s)]+2*B*phi[I(r-1,s)]+C*(phi[I(r,s-1)]+phi[I(r,s+1)])+D*(phi[I(r-1,s-1)]+phi[I(r-1,s+1)]);
  else
  {
    cerr<<"f6 impossible case"<< endl;
    exit(4);
  }
}

// define Gradient function
double Gradient(int r, int s, const vector<double> &phi)
{
  if((s==0)||(s==N))
    return 0.;
  
  double ep, si, ic, phi_b, sum=0.0;

    if(r>=0)
    {
      ep=e1;
      si=s1;
      ic=I1;
      phi_b=0.;
    }
    else
    {
      ep=e2;
      si=s2;
      ic=I2;
      phi_b=phi_D;
    }
    
    for(int k=0; k<=n; ++k)
    {
      for(int l=0; l<=2*k; ++l)
      {
	double T1=6.022e-4*2*ic*dx*dz*l_B*Q(k)*f1(r,s,k,l,phi_b,phi);
	double T3=6.022e-4*2*ic*dx*dz*l_B*Q(k)*f3(r,s,k,l,phi_b,phi);
	sum+=T1+T3;
      }
    }
   
    double T5=ep*f5(r,s,phi)/(4.*M_PI);
    sum+=T5;

    if(r>0)
    {
      ep=e1;
      si=s1;
      ic=I1;
      phi_b=0.;
    }
    else
    {
      ep=e2;
      si=s2;
      ic=I2;
      phi_b=phi_D;
    }
    
    for(int k=0; k<=n; ++k)
    {
      for(int l=0; l<=2*k; ++l)
      {
	double T2=6.022e-4*2*ic*dx*dz*l_B*Q(k)*f2(r,s,k,l,phi_b,phi);
	double T4=6.022e-4*2*ic*dx*dz*l_B*Q(k)*f4(r,s,k,l,phi_b,phi);
	sum+=T2+T4;
      }
    }

    double T6=ep*f6(r,s,phi)/(4.*M_PI);
    sum+=T6;

    return sum;
}


int main (int args, char *arg[])
{
   if(args!=19)
   {
      cerr<<"Program1 <e1> <e2> <I1> <I2> <dx> <dz> <phi_P> <phi_D> <M> <N> <file_phi_low> <file_phi_up> <tol> <iter_max> <iter_write> <l_B> <alpha> <n>" <<endl;
      return 5;
   }
   e1 = atof(arg[1]);
   e2 = atof(arg[2]);
   I1 = atof(arg[3]);
   I2 = atof(arg[4]);
   dx = atof(arg[5]);
   dz = atof(arg[6]);
   phi_P = atof(arg[7]);
   phi_D = atof(arg[8]);
   M = atoi(arg[9]);
   N = atoi(arg[10]);
   double tol = atof(arg[13]);
   int iter_max = atoi(arg[14]);
   int iter_write = atoi(arg[15]);
   l_B = atof(arg[16]);
   double alpha = atof(arg[17]);
   n = atoi(arg[18]);
   
   A=dz/(3.*dx)+dx/(3.*dz);
   B=-dz/(3.*dx)+dx/(6.*dz);
   C=dz/(6.*dx)-dx/(3.*dz);
   D=-dz/(6.*dx)-dx/(6.*dz);

  
  phi_low.assign(N+1, 0.);
  ifstream file_phi_low(arg[11]);   // Open file "arg[11]" for reading

   // Omit the header
   while (file_phi_low.peek() == '#')
     file_phi_low.ignore(numeric_limits<int>::max(), '\n');

   // Read in two columns
   for (int s = 0; s <= N; ++s)
     {
       double tmp1, tmp2;
       file_phi_low >> tmp1 >> tmp2;
       if (!file_phi_low.good())
	 {
	   cerr << "Error in file \"" << arg[11] << "\"" << endl;
           return 6;
         }
       phi_low[s] = tmp2;
     }
   
   file_phi_low.clear();
   file_phi_low.close();
   
   
   phi_up.assign(N+1, 0.); 
   ifstream file_phi_up(arg[12]);   // Open file "arg[12]" for reading

   // Omit the header
   while (file_phi_up.peek() == '#')
     file_phi_up.ignore(numeric_limits<int>::max(), '\n');

   // Read in two columns
   for (int s = 0; s <= N; ++s)
     {
       double tmp1, tmp2;
       file_phi_up >> tmp1 >> tmp2;
       if (!file_phi_up.good())
	 {
	   cerr << "Error in file \"" << arg[12] << "\"" << endl;
           return 6;
         }
       phi_up[s] = tmp2;
     }

   file_phi_up.clear();
   file_phi_up.close();

   
  vector<double> phi((2*M+1)*(N+1),0.),phi_new((2*M+1)*(N+1),0.);
  double delta;
  for (int r=-M; r<=M; ++r)
    for (int s=0; s<=N; ++s)
    {
      phi[I(r,s)] = phi_low[s]*((double)(M-r))/((double)(2*M)) + 
                    phi_up [s]*((double)(M+r))/((double)(2*M));
    }
  for(int iter=0; iter<=iter_max; ++iter)
  {
    if((iter%iter_write==0)&&(iter>0))
    {
      stringstream name;
      name << "data_" << iter/iter_write << ".txt";
      ofstream file (name.str().c_str());
      file.precision(numeric_limits<double>::digits10);
      file << "#e1 = " << e1 << endl;
      file << "#e2 = " << e2 << endl;
      file << "#I1 = " << I1 << endl;
      file << "#I2 = " << I2 << endl;
      file << "#dx = " << dx << endl;
      file << "#dz = " << dz << endl;
      file << "#phi_P = " << phi_P << endl;
      file << "#phi_D = " << phi_D << endl;
      file << "#M = " << M << endl;
      file << "#N = " << N << endl;
      file << "#tol = " << tol << endl;
      file << "#iter_max = " << iter_max << endl;
      file << "#iter_write = " << iter_write << endl;
      file << "#iter = " << iter << endl;
      file << "#l_B = " << l_B << endl;
      file << "#alpha = " << alpha << endl;
      file << "#delta = " << delta << endl;
      file << "#n = " << n << endl;
      file << "#" << endl;
      file << "# r s phi(r,s)" << endl;

      for(int r=-M; r<=M; ++r)
        for(int s=0; s<=N; ++s)
          {
          file << r << " " << s << " " << phi[I(r,s)] << endl;
          }
      file.clear();
      file.close();
    }
    delta=0.0;
    for(int r=-M; r<=M; ++r)
      for(int s=0; s<=N; ++s)
      {
        phi_new[I(r,s)]=phi[I(r,s)]-alpha*Gradient(r,s,phi);
        delta=max(delta,fabs(phi_new[I(r,s)]-phi[I(r,s)]));
      }
   if(delta<=tol)
      break;
   phi=phi_new;
  }

  cout.precision(numeric_limits<double>::digits10);
  cout << "#e1 = " << e1 << endl;
  cout << "#e2 = " << e2 << endl;
  cout << "#I1 = " << I1 << endl;
  cout << "#I2 = " << I2 << endl;
  cout << "#dx = " << dx << endl;
  cout << "#dz = " << dz << endl;
  cout << "#phi_P = " << phi_P << endl;
  cout << "#phi_D = " << phi_D << endl;
  cout << "#M = " << M << endl;
  cout << "#N = " << N << endl;
  cout << "#tol = " << tol << endl;
  cout << "#iter_max = " << iter_max << endl;
  cout << "#iter_write = " << iter_write << endl;
  cout << "#l_B = " << l_B << endl;
  cout << "#alpha = " << alpha << endl;
  cout << "#delta = " << delta << endl;
  cout << "#n = " << n << endl;
  cout << "#" << endl;
  cout << "# r s phi(r,s)" << endl;

  for(int r=-M; r<=M; ++r)
    for(int s=0; s<=N; ++s)
      {
      cout << r << " " << s << " " << phi[I(r,s)] << endl;
      }
    
  return 0;
}


 
