#include "itensor/all.h"
#include<time.h>
//#include "itensor/mps/sites/spinhalf.h"

using namespace itensor;
using namespace std;
using std::vector;


ITensor
Szu(const SiteSet& sites, int b)
{
    ITensor Szu = sites.op("projUp",b);
    return Szu;
}

ITensor
Szd(const SiteSet& sites, int b)
{
    ITensor Szd = sites.op("projDn",b);
    return Szd;
}

/*
ITensor
Sxu(const SiteSet& sites, int b)
{
    ITensor Sxu = sites.op("projxUp",b);
    return Sxu;
}

ITensor
Sxd(const SiteSet& sites, int b)
{
    ITensor Sxd = sites.op("projxDn",b);
    return Sxd;
}
*/

/*
vector<AutoMPO>
makempoH1(SiteSet const& sites)
    {
    auto N = sites.N();
    auto mpoH1 = vector<AutoMPO>(N+1);        // N+1 ??
    for(auto b : range1(N-1))
        {
          mpoH1.at(b) += 0.5,"S+",b,"S-",b+1;    // 0.5,"S+",j,"I",j+1,"S-",j+2;
          mpoH1.at(b) += 0.5,"S-",b,"S+",b+1;
          mpoH1.at(b) +=     "Sz",b,"Sz",b+1;
        }
    return mpoH1;
    }
*/



vector<ITensor>
makeH(SiteSet const& sites)
    {
    auto N = sites.N();
    auto H = vector<ITensor>(N+1);        // N+1 ??
    for(auto b : range1(N-1))
        {
        H.at(b) = sites.op("Sz",b)*sites.op("Sz",b+1)
                   + 0.5*sites.op("S+",b)*sites.op("S-",b+1)
                   + 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        }
    return H;
    }

vector<ITensor>
makeGates(vector<ITensor> const& H,
          Real tau)
    {
    auto gates = H;
    for(auto& g : gates)
        {
        if(!g) continue;
        g = expHermitian(-tau*g);
        }
    return gates;
    }

vector<ITensor>   //swap gates
makeSwap(SiteSet const& sites)
{
  auto N = sites.N();
  auto swap = vector<ITensor>(N+1);
  for(auto b : range1(N-1))
  {
  auto sn = sites(b);      auto snp = prime(sn);
  auto sn1 = sites(b+1);   auto sn1p = prime(sn1);
  swap.at(b) = ITensor(sn,snp,sn1,sn1p);         //2x2x2x2 tensor
  swap.at(b).set(sn(1),sn1(1),snp(1),sn1p(1),1);
  swap.at(b).set(sn(1),sn1(2),snp(2),sn1p(1),1);
  swap.at(b).set(sn(2),sn1(1),snp(1),sn1p(2),1);
  swap.at(b).set(sn(2),sn1(2),snp(2),sn1p(2),1);
  }
  return swap;
}


void
doSVD(MPS & psi,
      int b,
      ITensor phi,
      Direction dir,
      Real cutoff)
    {
    auto U = psi.A(b);
    ITensor D,V;
    svd(phi,U,D,V,{"Cutoff",cutoff});
    if(dir == Fromleft)
        {
        //multiply D into V
        psi.setA(b,U);
        psi.setA(b+1,D*V);
        }
    else
        {
        //multiply D into U
        psi.setA(b,U*D);
        psi.setA(b+1,V);
        }
    }

ITensor
applyGate(ITensor phi,
          ITensor gate)
    {
    //TODO:

    phi = phi*gate;     //1. Apply gate to phi using * operator

    phi.mapprime(1,0);      //2. Restore original prime level of phi's indices
    phi /= norm(phi);      //3. Normalize phi by dividing by norm(phi)

    return phi;
    }

void
applyswap(MPS & psi,
      int b,
      ITensor swap,
      Real cutoff)
  {
    auto newtwo = psi.A(b)*psi.A(b+1)*swap;
    newtwo.mapprime(1,0);      //2. Restore original prime level of phi's indices
    auto U = psi.A(b);
    ITensor D,V;
    svd(newtwo,U,D,V,{"Cutoff",cutoff});
    psi.setA(b,U*D);
    psi.setA(b+1,V);
    //PAUSE
  }

int
main()
{
  int N = 20;
  Real tau = 0.05;
  Real cutoff = 1E-12;
  int nsweep = 20;
  //int nsweep = 1;
  auto Energy = 0.0;
  Real FiniteTE[40];
  auto Energy1 = 0.0;
  auto Energybond = 0.0;
  auto Energybond1 = 0.0;
  Real Pszu[N],Pszd[N],Psxu[N],Psxd[N],Snx[N],Snz[N];

  auto sites = SpinHalf(N);
  auto state = InitState(sites);
  auto psi = MPS(state),psi0 = MPS(state);
  clock_t start, end;

  ofstream output,outfile_E,outfile_FiniteE,outfile_C,outfile_Sn,outfile_SS;

  output.open("output.txt");
  outfile_E.open("Energy.txt"); outfile_SS.open("myfile_SS.txt");outfile_FiniteE.open("ancillat_FE.txt");
  outfile_C.open("myfile_C.txt"); outfile_Sn.open("myfile_Sn.txt");


  ITensor D,V;
  for(int n=1; n<N+1 ;n+=2)
      {
        auto U = psi0.A(n);
        auto  T = ITensor(sites(n),sites(n+1));
        T.set(sites(n)(1),sites(n+1)(2),1/sqrt(2));
        T.set(sites(n)(2),sites(n+1)(1),-1/sqrt(2));
        //(auto psin,auto psin1) = factor(T,A,B);
        svd(T,U,D,V,{"Cutoff",cutoff});
        psi0.setA(n,U);      //psi.setA(n,U*D);
        psi0.setA(n+1,D*V);
      }

  auto H = makeH(sites);
  //auto mpoH1 = makempoH1(sites);
  auto gates = makeGates(H,tau/2);
  auto swap = makeSwap(sites);

  AutoMPO ampo(sites);

  for(int j = 1; j < (N-1);j=j+2)
  //for(int j = 1; j < N; ++j)
    {
    ampo += 0.5,"S+",j,"S-",j+2;    //0.5,"S+",j,"I",j+1,"S-",j+2;
    ampo += 0.5,"S-",j,"S+",j+2;
    ampo +=     "Sz",j,"Sz",j+2;
    }
  auto effH = MPO(ampo);

  Energy = psiHphi(psi,effH,psi);
  println(Energy);

//for (int t = 0.1; t <=1 ;t+=0.1)
//for (int t = 1; t <=80 ;t+=2)
for (int t = 1; t <=40 ;t+=1)
{
  Energy = 0.0;
  psi = psi0;
  //for(auto sw : range1(nsweep))
  //for(auto sw : range1(t*10))
  for(auto sw : range1(t*2))
      {
        if (sw == nsweep)
          {
            start = clock();
            cout<<"start "<<sw<<" "<<start/CLOCKS_PER_SEC<<" clocks per second "<<CLOCKS_PER_SEC<<endl;
          }
        Energy1 = 0.0;
        //println("Starting sweep ",sw);
        for(int b=1; b<(N-1) ;b+=2)           //1,3,5,7
        {
          //auto psiold = psi;
          psi.position(b);

          applyswap(psi,b+1,swap.at(b+1),cutoff);    //swap sites b,b+1
          auto phi = psi.A(b)*psi.A(b+1);
          //println("Bond dimension between real sites = ",commonIndex(psi.A(b),psi.A(b+1)).m());
          phi = applyGate(phi,gates.at(b));
          //auto bra = prime(phi,Site);
          //Energybond1 = (bra*H.at(b)*phi).real(); //measurement of energy <S.S>
          doSVD(psi,b,phi,Fromleft,cutoff);

          applyswap(psi,b+1,swap.at(b+1),cutoff);   //swap sites b,b+1   back

          //olap = overlap(psi,psiold);
          //println("overlap2  ",olap);
          psi /= sqrt(norm(psi));
        }
        //psi /= sqrt(norm(psi));
        Energy = psiHphi(psi,effH,psi);
        //println("  Half 1: Energy = ",Energy);
        //println("  Half 1: Energysum1 = ",Energy1);
        Energy1 = 0.0;
        for(auto b = N-2; b >= 2; b-=2)            // 8,6,4,2
            {
            //auto psiold = psi;
            psi.position(b);
            applyswap(psi,b,swap.at(b),cutoff);    //swap sites b-1,b
            //Real olap = overlap(psi,psiold);

            auto phi = psi.A(b-1)*psi.A(b);

            phi = applyGate(phi,gates.at(b-1));

            doSVD(psi,b-1,phi,Fromright,cutoff);
            applyswap(psi,b,swap.at(b),cutoff);    //swap sites b-1,b

            //olap = overlap(psi,psiold);
            //println("overlap2  ",olap);
            psi /= sqrt(norm(psi));
            }

            if (sw == nsweep)
              {
                end = clock();
                cout<<"end "<<sw<<" "<<end/CLOCKS_PER_SEC<<endl;
              }
      Energy = psiHphi(psi,effH,psi);
      outfile_E<<sw<<" "<<Energy<<endl;
      //if (sw == nsweep)
      if (sw == 0)
      {
      for(int n=1;n<=N;n++)                      //do collapse
        {
          psi.position(n);                        //???
          ITensor wf = psi.A(n);
          //ITensor wf = psi.A(n)*psi.A(n+1);
          //Pszu[n]=(dag(prime(wf,Site))*Szu(sites,n)*wf).real();
          //Psxu[n]=(dag(prime(wf,Site))*Sxu(sites,n)*wf).real();
          Snx[n] = (dag(prime(wf,Site))*sites.op("Sx",n)*wf).real();
          Snz[n] = (dag(prime(wf,Site))*sites.op("Sz",n)*wf).real();
          //Pszd[n]=(dag(prime(wf,Site))*Szd(sites,n)*wf).real();
          //Psxd[n]=(dag(prime(wf,Site))*Sxd(sites,n)*wf).real();
          cout<<n<<" "<<Snx[n]<<" "<<Snz[n]<<endl;
          //outfile_Sn<<n<<" "<<Pszu[n]<<" "<<Pszd[n]<<endl;
        }
      }

      }
FiniteTE[t] = Energy;
}

Energy = psiHphi(psi,effH,psi);
println("  whole: Energy = %.10f ",Energy);
for(int t = 1; t<=40;t++)
{
  //println("temperature:%1d E: %.10f ",t*10, FiniteTE[t]);
  println("temperature:%.10f E: %.10f ",t*0.2, FiniteTE[t]);
  outfile_FiniteE<<t*0.2<<" "<<FiniteTE[t]<<endl;
}
outfile_FiniteE.close();
  return 0;
}
