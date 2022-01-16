#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cmath>
#include <complex>
using namespace std;

double drand(){
    return rand()/(RAND_MAX - 1.0);
}

int angle_chooser(int matrix_size){
    return int(matrix_size*drand());
}

int tester(double change){
    if(drand()<=exp(-change)){
        return 1;
    }
    else {
        return 0;
    }
}


int main()
{   const   complex<double> I(0.0,1.0);
    double g = 0.2;
    double eps = 1;
    std::complex<double> S_old = {0,0}; // Action variable
    std::complex<double> S_new = {0,0};
    std::complex<double> dS; // Change in action with small change in unitary matrix
    int N = 50;
    double theta[N]; // Rows for real and imaginary parts and column for thetas.
    int THERMS = 5000;
    int SWEEPS = 10000;
    int count = 0;
    int index = 0;
    std::complex<double> temp;
    double polyakov = 0;
    int k = 0;
    double p_total = 0;
fstream polya;
polya.open("polyakov_values.txt",ios::out);
while(g<=3.2){

    for(int i = 0;i<N;i++){
        theta[i] = 0.01*drand();
        //cout<<theta[i];

    }
    for(int i =0;i<N;i++){
        S_old = S_old + -N*g*(exp(I*theta[i])+exp(-I*theta[i]))*0.5;
        for(int j = 0;j<N;j++){
            if(i!=j){
                S_old = S_old - log(sin((abs(theta[i]-theta[j])*0.5)));

            }

        }
    }

    //cout<<S_old;
    //cout<<"thermalizing\n";
    //thermalization
    while(count<THERMS){
        index = angle_chooser(N);

         for(int i =0;i<N;i++){
        if(i==index){
            temp = theta[i] + 2*(drand()-0.5)*eps ;
           S_new = S_new + -N*g*(exp(I*temp)+exp(-I*temp))*0.5;
           for(int l = 0;l<N;l++){
            if(l!=i&l!=index){
            S_new = S_new - log(sin(abs((temp-theta[l])*0.5)));}
           }
           continue;}

        S_new = S_new + -N*g*(exp(I*theta[i])+exp(-I*theta[i]))*0.5;
        for(int j = 0;j<N;j++){
            if(i!=j&j!=index){
               S_new = S_new - log(sin(abs((theta[i]-theta[j])*0.5)));
            }

        }}
        dS = real(S_new - S_old);
        k = tester(real(dS));
        if(k==1){
            theta[index] = real(temp);
            S_old = S_new;

        }
        count++;
        S_new = {0,0};
    }
    count = 0;
    cout<<"computing polyakov loop\n";
    int monte_carlo_time = 0;
    int accept = 0;

//fstream acceptance;
    //acceptance.open("Acceptance_values_g=1.txt",ios::out);

    //Sweeps
    while(count<SWEEPS){
            index = angle_chooser(N);

         for(int i =0;i<N;i++){
        if(i==index){
            temp = theta[i]+2*(drand()-0.5)*eps ;
           S_new = S_new + -N*g*(exp(I*temp)+exp(-I*temp))*0.5;
           for(int l = 0;l<N;l++){
            if(l!=i&l!=index){
            S_new = S_new - log(sin(abs((temp-theta[l])*0.5)));}
           }
           }

        S_new = S_new + -N*g*(exp(I*theta[i])+exp(-I*theta[i]))*0.5;
        for(int j = 0;j<N;j++){
            if(i!=j&j!=index){
                S_new = S_new - log(sin(abs((theta[i]-theta[j])*0.5)));;
            }

        }}
        dS = real(S_new - S_old);
        k = tester(real(dS));
        if(k==1){
            theta[index] = real(temp);
            S_old = S_new;

            accept++;
            //cout<<accept<<"\n";

        }
        temp = 0;
            for(int i = 0;i<N;i++){
                temp = temp + exp(I*theta[i]);
            }
            polyakov = polyakov + real(temp);
        count++;
        S_new = {0,0};
        monte_carlo_time++;
        //if(monte_carlo_time%100==0){
           //acceptance<<monte_carlo_time<<"\t"<<100.0*accept/monte_carlo_time<<"\n";

        //}

    }
    //acceptance.close();
    polyakov = 1.0*polyakov/(N*SWEEPS);
    polya<<g<<"\t"<<polyakov<<"\n";
    cout<<"polyakov"<<polyakov;
    polyakov = 0;
    g = g + 0.5;
    }
polya.close();
    return 0;
    }






