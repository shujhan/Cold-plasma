#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <iomanip> 
#include <chrono>
#include <vector>
#include <omp.h>
using namespace std;
using namespace std::chrono;

vector<double> E_field(vector<double> &x_target, const vector<double> &x_source, int omega_0, double delta, vector<double> &wt);
void Euler(vector<double> &x, vector<double> &v, vector<double> &x_passive, vector<double> &v_passive, vector<double> &wt,double dt,int omega_0, double delta);
void RK4(vector<double> &x, vector<double> &v, vector<double> &x_passive, vector<double> &v_passive, vector<double> &wt,double dt,int omega_0, double delta);
void insertion(vector<double> &x, vector<double> &v, vector<double> &alpha, vector<double> &x_passive, vector<double> &v_passive, vector<double> &alpha_passive, vector<double> &wt, double d1);
void write(vector<double> &x, vector<double> &v);
// void setwt();




int main(int argc, char* argv[]) {

    double N = 200;

    // TODO: change it to double 
    const double epsilon = 0.05; 
    const double delta = 0.04;
    const double d1 = 0.05;
    const double dt = 0.04;

    const double method = 2; // 1: euler 2: RK4
    const double t_final = 20;


    // const int N = stoi(argv[1]);

    // // TODO: change it to double 
    // const double epsilon = stoi(argv[2]); 
    // const double delta = stoi(argv[3]);
    // const double d1 = stoi(argv[4]);
    // const double dt = stoi(argv[5]);

    // const double method = stoi(argv[6]); // 1: euler 2: RK4
    // const double t_final = stoi(argv[7]);


    const double omega_0 = 1.0;
    vector<double> alpha_passive(N + 1);
    vector<double> x_passive(N + 1);
    vector<double> v_passive(N + 1);

    vector<double> alpha(N);
    vector<double> x(N);
    vector<double> v(N);

    vector<double> wt(N);

    // Initialization
    for (int i = 0; i <= N; i++) {
        alpha_passive[i] = i / N;
     //   cout << alpha_passive[i] << "!";
	x_passive[i] = alpha_passive[i] + epsilon * sin(2 * M_PI * alpha_passive[i]);
        v_passive[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        alpha[i] = 0.5 * (alpha_passive[i] + alpha_passive[i + 1]);
        x[i] = alpha[i] + epsilon * sin(2 * M_PI * alpha[i]);
        v[i] = 0;
    }

    // set weights
    for (int i = 0; i < N; ++i) {
        wt[i] = alpha_passive[i + 1] - alpha_passive[i];
    }

    int Nstep = static_cast<int>(t_final / dt);


   // cout << "\n";	
   // cout << "alpha: ";
   // for (int i = 0; i < x.size(); i++) {
     //   cout << alpha[i] << " ";
  //  }
   // cout << "\n";
   // cout << "alpha_pass: ";
   // for (int i = 0; i < x.size()+1; i++) {
     //   cout << alpha_passive[i] << " ";
   // }
   // cout << "\n";






   // cout << "x: ";
   // for (int i = 0; i < x.size(); i++) {
     //   cout << x[i] << " ";
   // }
   // cout << "\n";

   // cout << "v: ";
   // for (int i = 0; i < v.size(); i++) {
     //   cout << v[i] << " ";
   // }



    auto start = high_resolution_clock::now();
    for (int step = 1; step <= Nstep; step++) {
        if (method == 1) {
            Euler(x, v, x_passive, v_passive, wt,dt,omega_0,delta);
        }

        if (method == 2) {
            RK4(x, v, x_passive, v_passive, wt,dt,omega_0,delta);
        }

        if (step == Nstep){
           write(x, v);
	   cout << x.size() << endl;
        } 
        // adaptive insertion
        insertion(x,v,alpha,x_passive,v_passive,alpha_passive,wt,d1);
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end-start);
    cout << "time used: " << duration.count() << " milliseconds"  <<  endl;


    return 0;
}

void write(vector<double> &x, vector<double> &v) {

    // Open a file for writing
    std::ofstream outputFile("vector.txt");

    // Check if the file is opened successfully
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open the file." << std::endl;
        exit(1);
    }

    // Write vector elements to the file
    outputFile << "x: ";
    for (int i = 0; i < x.size(); i++) {
        outputFile << x[i] << " ";
    }
    outputFile << "\n";

    outputFile << "v: ";
    for (int i = 0; i < v.size(); i++) {
        outputFile << v[i] << " ";
    }

    // Close the file
    outputFile.close();

    std::cout << "Vector elements written to the file successfully." << std::endl;
}


vector<double> E_field(vector<double> &x_target, const vector<double> &x_source, int omega_0, double delta, vector<double> &wt) {
    // int target_num = x_target.size();
    int source_num = x_source.size();
    int target_num = x_target.size();
    //vector<double> xsm(source_num);
    //vector<double> xtm(target_num);
    //for (int i = 0; i < source_num; ++i) {
      //  xsm[i] = fmod(x_source[i],1);
    //}
   // for (int i = 0; i < target_num; ++i) {
     //   xtm[i] = fmod(x_target[i],1);
   // }

    vector<double> E_Field(target_num);


    #pragma omp parallel for
    for (int i = 0; i < target_num; ++i) {
        for (int j = 0; j < source_num; ++j) {
            double diff = x_target[i] - x_source[j];
            while (diff < -0.5) {
                diff = diff+1;
            }
            while (diff >= 0.5) {
                diff = diff-1;
            }
          
	    
            double c_delta = sqrt(1+4*delta*delta);
            E_Field[i] -= (c_delta/2*diff/sqrt(diff*diff+delta*delta)-diff) * omega_0 * wt[j];
        }
    }

    return E_Field;

}

void Euler(vector<double> &x, vector<double> &v, vector<double> &x_passive, vector<double> &v_passive, vector<double> &wt,double dt,int omega_0, double delta) {
    int active_num = x.size();
    int passive_num = x_passive.size();
    vector<double> dx(active_num);
    vector<double> dxp(passive_num);
   for (int i=0; i < active_num; i++) {
        dx[i] = v[i];
        dxp[i] = v_passive[i];
   }
    vector<double> dv(active_num);
    vector<double> dvp(passive_num);
    dv = E_field(x,x,omega_0,delta,wt);
    dvp = E_field(x_passive,x,omega_0,delta,wt);
    for (int i = 0; i < active_num; i++) {
        x[i] += dt * dx[i];
        v[i] += dt * dv[i];
        x_passive[i] += dt * dxp[i];
        v_passive[i] += dt * dvp[i];
    }
    x_passive[passive_num-1] = x_passive[0] +1;
    v_passive[passive_num-1] = v_passive[0];
}


void RK4(vector<double> &x, vector<double> &v, vector<double> &x_passive, vector<double> &v_passive, vector<double> &wt,double dt,int omega_0, double delta) {
    int active_num = x.size();
    int passive_num = x_passive.size();

    vector<double> dx1(active_num);
    vector<double> dxp1(active_num);
    for (int i=0; i < active_num; i++) {
            dx1[i] = v[i];
            dxp1[i] = v_passive[i];
    }
    vector<double> dv1(active_num);
    vector<double> dvp1(active_num);
    vector<double> dx2(active_num);
    vector<double> dxp2(active_num);
    vector<double> dv2(active_num);
    vector<double> dvp2(active_num);
    vector<double> dx3(active_num);
    vector<double> dxp3(active_num);
    vector<double> dv3(active_num);
    vector<double> dvp3(active_num);
    vector<double> dx4(active_num);
    vector<double> dxp4(active_num);
    vector<double> dv4(active_num);
    vector<double> dvp4(active_num);

    vector<double> xold(active_num);
    vector<double> xpold(active_num);
    vector<double> vold(active_num);
    vector<double> vpold(active_num);

    dv1 = E_field(x,x,omega_0,delta,wt);
    dvp1 = E_field(x_passive,x,omega_0,delta,wt);

    for (int i = 0; i < active_num; i++) {
        xold[i] = x[i] + dt * dx1[i]/2;
        vold[i] = v[i] + dt * dv1[i]/2;
        xpold[i] = x_passive[i] + dt * dxp1[i]/2;
        vpold[i] = v_passive[i] + dt * dvp1[i]/2;
    }
    
    for (int i=0; i < active_num; i++) {
            dx2[i] = vold[i];
            dxp2[i] = vpold[i];
    }


    dv2 = E_field(xold,xold,omega_0,delta,wt);
    dvp2 = E_field(xpold,xold,omega_0,delta,wt);
    for (int i = 0; i < active_num; i++) {
        xold[i] = x[i] + dt * dx2[i]/2;
        vold[i] = v[i] + dt * dv2[i]/2;
        xpold[i] = x_passive[i] + dt * dxp2[i]/2;
        vpold[i] = v_passive[i] + dt * dvp2[i]/2;
    }


    for (int i=0; i < active_num; i++) {
            dx3[i] = vold[i];
            dxp3[i] = vpold[i];
    }


    dv3 = E_field(xold,xold,omega_0,delta,wt);
    dvp3 = E_field(xpold,xold,omega_0,delta,wt);
    for (int i = 0; i < active_num; i++) {
        xold[i] = x[i] + dt * dx3[i];
        vold[i] = v[i] + dt * dv3[i];
        xpold[i] = x_passive[i] + dt * dxp3[i];
        vpold[i] = v_passive[i] + dt * dvp3[i];
    }

    for (int i=0; i < active_num; i++) {
            dx4[i] = vold[i];
            dxp4[i] = vpold[i];
    }


    dv4 = E_field(xold,xold,omega_0,delta,wt);
    dvp4 = E_field(xpold,xold,omega_0,delta,wt);
    for (int i = 0; i < active_num; i++) {
        x[i] += dt * (dx1[i]+ 2* dx2[i] + 2*dx3[i] + dx4[i])/6;
        v[i] += dt * (dv1[i]+ 2* dv2[i] + 2*dv3[i] + dv4[i])/6;
        x_passive[i] += dt * (dxp1[i]+ 2* dxp2[i] + 2*dxp3[i] + dxp4[i])/6;
        v_passive[i] += dt * (dvp1[i]+ 2* dvp2[i] + 2*dvp3[i] + dvp4[i])/6;
    }

    x_passive[active_num] = x_passive[0] +1;
    v_passive[active_num] = v_passive[0];

}




void insertion(vector<double> &x, vector<double> &v, vector<double> &alpha, vector<double> &x_passive, vector<double> &v_passive, vector<double> &alpha_passive, vector<double> &wt, double d1) {
    int num_interval = x.size();
    int count = 1;
    int count_passive = 1;
    vector<double> xnew;
    vector<double> vnew;
    vector<double> anew;
    vector<double> xpnew;
    vector<double> vpnew;
    vector<double> apnew;

    


    for (int i = 0; i < num_interval; i++) {
        double a0 = alpha_passive[i];
        double a1 = alpha[i];
        double a2 = alpha_passive[i+1];
    
        double x0 = x_passive[i];
        double x1 = x[i];
        double x2 = x_passive[i+1];

        double v0 = v_passive[i];
        double v1 = v[i];
        double v2 = v_passive[i+1];

        //  euclidean distance in phase sapace   
        double dist = sqrt((x2-x0)*(x2-x0) + (v2-v0)*(v2-v0));
        // dist2 = sqrt((x1-(x2+x0)/2)^2 + (v1-(v2+v0)/2)^2);
        if (dist > d1) {
            // add first point to passive
            apnew.push_back(a0);
            xpnew.push_back(x0);
            vpnew.push_back(v0);

            apnew.push_back(a1);
            xpnew.push_back(x1);
            vpnew.push_back(v1);

            double middle_alpha1 = 0.5*(a1 + a0);
            double middle_alpha2 = 0.5*(a1 + a2);
            //  use newtons form interpolation
            double a = ((x2-x1)/(a2-a1) - (x1-x0)/(a1-a0))/(a2-a0);
            double b = (x1-x0)/(a1-a0);
            double c = x0;
            double insert_x1 = a*(middle_alpha1 - a0)*(middle_alpha1 - a1) + b*(middle_alpha1-a0) + c;
            double insert_x2 = a*(middle_alpha2 - a0)*(middle_alpha2 - a1) + b*(middle_alpha2-a0) + c;
            

            double a_v = ((v2-v1)/(a2-a1) - (v1-v0)/(a1-a0))/(a2-a0);
            double b_v = (v1-v0)/(a1-a0);
            double c_v = v0;
            double insert_v1 = a_v*(middle_alpha1 - a0)*(middle_alpha1 - a1) + b_v*(middle_alpha1-a0) + c_v;
            double insert_v2 = a_v*(middle_alpha2 - a0)*(middle_alpha2 - a1) + b_v*(middle_alpha2-a0) + c_v;
            
            anew.push_back(middle_alpha1);
            xnew.push_back(insert_x1);
            vnew.push_back(insert_v1);

            anew.push_back(middle_alpha2);
            xnew.push_back(insert_x2);
            vnew.push_back(insert_v2);
        }
        else {
            apnew.push_back(a0);
            xpnew.push_back(x0);
            vpnew.push_back(v0);

            anew.push_back(a1);
            xnew.push_back(x1);
            vnew.push_back(v1);
        }
    }
    // add the last point
    apnew.push_back(apnew[0]+1);
    xpnew.push_back(xpnew[0]+1);
    vpnew.push_back(vpnew[0]);

    x.clear();
    v.clear();
    alpha.clear();
    wt.clear();

    x_passive.clear();
    v_passive.clear();
    alpha_passive.clear();

    x.reserve(anew.size());
    v.reserve(anew.size());
    alpha.reserve(anew.size());
    wt.reserve(anew.size());

    x_passive.reserve(apnew.size());
    v_passive.reserve(apnew.size());
    alpha_passive.reserve(apnew.size());

    for (int i = 0; i < anew.size(); i++) {
        alpha.push_back(anew[i]);
        x.push_back(xnew[i]);
        v.push_back(vnew[i]);

        alpha_passive.push_back(apnew[i]);
        x_passive.push_back(xpnew[i]);
        v_passive.push_back(vpnew[i]);

    }


    alpha_passive.push_back(apnew[0]+1);
    x_passive.push_back(xpnew[0]+1);
    v_passive.push_back(vpnew[0]);

    for (int i = 0; i < alpha.size(); ++i) {
        wt.push_back(alpha_passive[i + 1] - alpha_passive[i]);
    }
}
