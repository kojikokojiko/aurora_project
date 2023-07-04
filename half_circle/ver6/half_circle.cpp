#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <iostream>
#include<fstream>
#include <numeric>
#include "MT.h"
#include <vector>
#include <tuple>
// #include <direct.h>
#include <sys/stat.h>
#include <sstream>
#include <random>
#include <iomanip> // for std::setprecision and std::fixed

using namespace std;
#define PI 3.14159265359


// 関数宣言////////////////////////////////////////////////////////


void make_mask(int mask_size, vector<int> &neighbor_row, vector<int> &neighbor_col);
double  Maxwell_velocity(double temp,double m);

void initialize_va(vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy,int NP);

void initialize_r(vector<double> &rx,vector<double> &ry, double TEMP_RHO);

void make_pairlist(vector<int> mask_col,vector<int> mask_row,
int select_gx,int select_gy,
vector<vector<int>> G_MAP,
vector<vector<int>> &pairlist);


void gmap_create4(int NP,vector<vector<int>> &g_map,
vector<vector<int>> &pairlist,
vector<double> rx,
vector<double> ry,
vector<int> mask_row,
vector<int> mask_col);


void plane_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter);

void bottom_right_corner_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter);

void upper_right_corner_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter);

void curve_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter);



void force(int NP,vector<double> &ax,vector<double> &ay,
vector<vector<int>> pairlist,
vector<double> rx,vector<double> ry,vector<double> diameter);

int BookKeeping(int NP,vector<double> Vx,vector<double> Vy);
///////////////////////////////////////////////////////////////


// int ver=3;
const double STATIC_DIA=50.0;
// const double TEMP_RHO=0.5;
const double L=STATIC_DIA*2;
const double LX=L;
const double LY=L;
const double KBT=1.0;
const double M=1.0;
const double DT=1e-5;
const double NSTEPS=1e7;
const double AFTER_POS_OUT_PERENT=0.0;
const int POS_OUT_NUM=500; 
const int POS_OUT_PERIOD=(1-AFTER_POS_OUT_PERENT)*NSTEPS/POS_OUT_NUM;

const double EPSILON=1.0;
// グリッド関連

const double LEN_BOX_X=0.4;
const double LEN_BOX_Y=LEN_BOX_X;
const int N_GX=ceil(LX/LEN_BOX_X);
const int N_GY=ceil(LY/LEN_BOX_Y);
const int N_G_ALL=N_GX*N_GY;

// マスク関連
const int MASK_SIZE=5;
const int MASK_LENGTH=2*(MASK_SIZE*(MASK_SIZE+1));
const int PAIR_LENGTH=15;

const double R_CUT=pow(2.0,1.0/6.0)*1.0;


// ABP関連
// const double V0=1.0;
// particle mobility.
const double MU=1.0;
// rotational diffusion coefficient.
const double ETA=5e-3;

const double SMALL_DIA=1.0;
int main(int argc, char *argv[]){
    const double TEMP_RHO =stod(argv[1]);
    const double V0=stod(argv[2]);
    

    // Append to the string stream
    std::ostringstream strs;
    strs << "rho_";
    strs << std::fixed << std::setprecision(2) << TEMP_RHO;
    strs << "_v0_";
    strs << std::fixed << std::setprecision(1) << V0;

    string ver=strs.str();
    // const string ver="test";
    cout<<ver<<endl;

    int reusableCount=0;

    // 乱数設定
    // random_device seed=mt19937::default_seed;
    mt19937 engine(123);
    double mean=0.0;
    // これは標準偏差
    double deviation=1.0;

    normal_distribution<> dist(mean,deviation);




    // 変数初期値
    vector<double> rx;
    vector<double> ry;
    vector<double> vx;    
    vector<double> vy;
    vector<double> WCA_x;
    vector<double> WCA_y;
    vector<double> theta;
    vector<double> diameter;
    vector<int> mask_row(MASK_LENGTH);
    vector<int> mask_col(MASK_LENGTH);

    make_mask(MASK_SIZE, mask_row, mask_col);


   

    vector<vector<int>> g_map(N_GY,vector<int>(N_GX));
    vector<vector<int>> pairlist;

    const string main_dir="./"+ver;
    const char *dir_name=main_dir.c_str();
    //mkdir(dir_name);
    mkdir(dir_name,0755);
    const string rx_file=main_dir+"/rx.dat";
    const string ry_file=main_dir+"/ry.dat";
    const string vx_file=main_dir+"/vx.dat";
    const string vy_file=main_dir+"/vy.dat";
    const string data_file=main_dir+"/data.dat";
    
    

    ofstream rxofs(rx_file);
    ofstream ryofs(ry_file);
    ofstream vxofs(vx_file);
    ofstream vyofs(vy_file);
    ofstream dataofs(data_file);
    time_t start_time,end_time;
    start_time=time(NULL);
   

    initialize_r(rx,ry,TEMP_RHO);


    int N=rx.size();

    initialize_va(WCA_x,WCA_y,vx,vy,N);


    // cout<<Nx*Ny<<endl;

    diameter.assign(N,SMALL_DIA);

    pairlist.resize(N,vector<int>(PAIR_LENGTH));
    
    
    // theta初期化
    theta.resize(N);
    for (int i=0;i<N;i++){
        theta[i]=2*PI*dist(engine);
    }


    gmap_create4(N,g_map,pairlist,rx,ry,mask_row,mask_col);



    force(N,WCA_x,WCA_y,pairlist,rx,ry,diameter);
    curve_force(WCA_x,WCA_y,rx,ry,N,diameter);
    plane_force(WCA_x,WCA_y,rx,ry,N,diameter);
    upper_right_corner_force(WCA_x,WCA_y,rx,ry,N,diameter);
    bottom_right_corner_force(WCA_x,WCA_y,rx,ry,N,diameter);




    dataofs<<"N "<<N<<endl;
    dataofs<<"LX "<<LX<<endl;
    dataofs<<"LY "<<LY<<endl;
    dataofs<<"DT "<<DT<<endl;
    dataofs<<"NSTEPS "<<NSTEPS<<endl;
    dataofs<<"POS_OUT_PERIOD "<<POS_OUT_PERIOD<<endl;
    dataofs<<"EPSILON "<<EPSILON<<endl;
    dataofs<<"R_CUT "<<R_CUT<<endl;
    dataofs<<"V0 "<<V0<<endl;
    dataofs<<"MU "<<MU<<endl;
    dataofs<<"ETA "<<ETA<<endl;
    dataofs<<"SMALL_DIA "<<SMALL_DIA<<endl;
    const double NU=(N*PI*SMALL_DIA*SMALL_DIA)/(LX*LY-PI*STATIC_DIA*STATIC_DIA/8.0);
    dataofs<<"NU "<<NU<<endl;
    const double REAL_RHO=(N)/(LX*LY-PI*STATIC_DIA*STATIC_DIA/8.0);
    dataofs<<"REAL_RHO "<<REAL_RHO<<endl;
    
    // メインループ
    for(int t=0;t<NSTEPS;t++){
        if(t%1000==0){
            cout<<t<<endl;
        }

        /////////// 速度ベルれ法１//////////////////////////////
        // 静止粒子には適用しないので1から
        
        for( int i=0;i<N;i++){
            
            rx[i]=rx[i]+vx[i]*DT;
            ry[i]=ry[i]+vy[i]*DT;
            vx[i]=V0*cos(theta[i])+MU*WCA_x[i];
            vy[i]=V0*sin(theta[i])+MU*WCA_y[i];
            // theta更新
            theta[i]=theta[i]+sqrt(2*ETA)*dist(engine)*DT;
            // 周期境界
            // 原点中心なので注意
            if(rx[i]>LX/2.0){
                rx[i]-=LX;
                if(rx[i]>LX/2.0){
                    cout<<"error"<<endl;
                }
            }
            else if(rx[i]<=-LX/2.0){
                rx[i]+=LX;
                if (rx[i]<=-LX/2.0){
                    cout<<"error"<<endl;
                }
            }

            if(ry[i]>LY/2.0){
                ry[i]-=LY;
                if(ry[i]>LY/2.0){
                    cout<<"error"<<endl;
                }
            }
            else if(ry[i]<=-LY/2.0){
                ry[i]+=LY;
                if (ry[i]<=-LY/2.0){
                    cout<<"error"<<endl;
                }
            }
        }
        

        // BookKeeping./////////////////////////////////////////
        if(reusableCount==0){
            gmap_create4(N,g_map,pairlist,rx,ry,mask_row,mask_col);
            reusableCount=BookKeeping(N,vx,vy);
            cout<<"reusableCount "<<reusableCount<<endl;


        }
        else{
            reusableCount-=1;
        }


        // for (int i = 0; i < N; ++i) {
        //     for (int j = 0; j < PAIR_LENGTH; ++j) {
        //         std::cout << pairlist[i][j] << ' ';
        //     }
        //     std::cout << '\n';
        // }

        WCA_x.assign(N,0.0);
        WCA_y.assign(N,0.0);
        force(N,WCA_x,WCA_y,pairlist,rx,ry,diameter);
        // for(int i=0;i<N;i++){
        //     cout<<WCA_x[i]<<" "<<WCA_y[i]<<endl;
        // }
        curve_force(WCA_x,WCA_y,rx,ry,N,diameter);
        plane_force(WCA_x,WCA_y,rx,ry,N,diameter);
        upper_right_corner_force(WCA_x,WCA_y,rx,ry,N,diameter);
        bottom_right_corner_force(WCA_x,WCA_y,rx,ry,N,diameter);


// end_BookKeeping./////////////////////////////////////////

        if((t>NSTEPS*AFTER_POS_OUT_PERENT)&&(t%POS_OUT_PERIOD==0)){
            for(size_t i=0;i<rx.size();i++){
                rxofs<<rx[i]<<" ";
                ryofs<<ry[i]<<" ";
                vxofs<<vx[i]<<" ";
                vyofs<<vy[i]<<" ";
            }
            rxofs<<endl;
            ryofs<<endl;
            vxofs<<endl;
            vyofs<<endl;
        }
    }


    rxofs.close();
    ryofs.close();
    vxofs.close();
    vyofs.close();
    dataofs.close();

    end_time=time(NULL);

    cout<<end_time-start_time<<endl;
    return 0;
}






int BookKeeping(int NP,vector<double> Vx,vector<double> Vy){
    double maxCandidateV;
    double maxV=0.0;
    double tLim;
    int reusableNum;
    for (int i=0;i<NP;i++){
        maxCandidateV=sqrt(Vx[i]*Vx[i]+Vy[i]*Vy[i]);
        if(maxCandidateV>maxV){
            maxV=maxCandidateV;
        }
        // R_CUTちゅうい
        tLim=(MASK_SIZE*LEN_BOX_X-R_CUT)/(2*maxV);
        reusableNum=int(tLim/(1.5*DT));
    }
    return reusableNum;
}

void force(int NP,vector<double> &ax,vector<double> &ay,
vector<vector<int>> pairlist,
vector<double> rx,vector<double> ry,vector<double> diameter)
{
    int roop_num;
    int pair_index;
    double rxij;
    double ryij;
    double r2;
    double ir,ir3;
    double ir2,ir6,ir12;
    double sigma;
    double sigma3;
    double sigma6,sigma12;
    double rc,rc2;

    double fx,fy;
    

    for (int i=0;i<NP;i++){
        roop_num=pairlist[i][PAIR_LENGTH-1];
        
        for (int j=0;j<roop_num;j++){
            pair_index=pairlist[i][j];
            rxij=rx[i]-rx[pair_index];
            
            // minimum image convention
            if (rxij>=LX/2){
                rxij=rxij-LX;
            }
            else if (rxij<-LX/2){
                rxij=rxij+LX;
            }
            else{
                rxij=rxij;
            }
            ryij=ry[i]-ry[pair_index];
            if (ryij>=LY/2){
                ryij=ryij-LY;
            }
            else if (ryij<-LY/2){
                ryij=ryij+LY;
            }
            else{
                ryij=ryij;
            }

            r2=rxij*rxij+ryij*ryij;
            ir2=1.0/r2;
            ir=sqrt(ir2);
            ir3=ir*ir*ir;
            ir6=ir2*ir2*ir2;
            ir12=ir6*ir6;

            sigma=(diameter[i]+diameter[pair_index])/2.0;
            sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
            sigma12=sigma6*sigma6;
            
            rc=pow(2.0,1.0/6.0)*sigma;
            rc2=rc*rc;

            
            if(r2>=rc2){
                fx=0.0;
                fy=0.0;
            }
            else{
                fx=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*rxij;
                fy=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*ryij;
                
            }
            // cout<<fx<<" "<<fy<<endl;
            ax[i]+=fx;
            ay[i]+=fy;
            ax[pair_index]-=fx;
            ay[pair_index]-=fy;
        }
    }
}



void curve_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter){

    double rxij;
    double ryij;
    double r2;
    double ir,ir3;
    double ir2,ir6,ir12;
    double sigma;
    double sigma3;
    double sigma6,sigma12;
    double rc,rc2;
    double fx,fy;
    
    // double rc2=(STATIC_DIA+SMALL_DIA);


    // Obstacleの中心は0,0なのでディスプレースメントは計算しない、してもいいけど

    for (int i=0;i<NP;i++){
        if (rx[i]<0.0){
            rxij=rx[i]-0.0;
            ryij=ry[i]-0.0;
            r2=rxij*rxij+ryij*ryij;
            ir2=1.0/r2;
            ir=sqrt(ir2);
            ir3=ir*ir*ir;
            ir6=ir2*ir2*ir2;
            ir12=ir6*ir6;
            sigma=(diameter[i]+STATIC_DIA)/2.0;
            sigma3=sigma*sigma*sigma;
            sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
            sigma12=sigma6*sigma6;   
            rc=pow(2.0,1.0/6.0)*sigma;
            rc2=rc*rc;
            if(r2>=rc2){
                fx=0.0;
                fy=0.0;
            }else{
                fx=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*rxij;
                fy=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*ryij;
                if (abs(fx)>100.0){
                    cout<<"Cylinder:fx"<<fx<<endl;
                    cout<<"rxij"<<rxij<<endl;
                }
                if (abs(fy)>100.0){
                    cout<<"Cylinder:fx"<<fx<<endl;
                    cout<<"rxij"<<rxij<<endl;
                }
            }
            ax[i]+=fx;
            ay[i]+=fy;
        }   
    }
}






void upper_right_corner_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter){

    double rxij;
    double ryij;
    double r2;
    double ir,ir3;
    double ir2,ir6,ir12;
    double sigma;
    double sigma3;
    double sigma6,sigma12;
    double rc,rc2;
    double fx,fy;

    double upper_right_corner_x=0.0;
    double upper_right_corner_y=STATIC_DIA/2.0;
    
    // double rc2=(STATIC_DIA+SMALL_DIA);


    // Obstacleの中心は0,0なのでディスプレースメントは計算しない、してもいいけど

    for (int i=0;i<NP;i++){
        if (rx[i]>=0.0){
            if (ry[i]>=STATIC_DIA/2.0){
                rxij=rx[i]-upper_right_corner_x;
                ryij=ry[i]-upper_right_corner_y;
                r2=rxij*rxij+ryij*ryij;
                ir2=1.0/r2;
                ir=sqrt(ir2);
                ir3=ir*ir*ir;
                ir6=ir2*ir2*ir2;
                ir12=ir6*ir6;
                sigma=diameter[i]/2.0;
                sigma3=sigma*sigma*sigma;
                sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
                sigma12=sigma6*sigma6;   
                rc=pow(2.0,1.0/6.0)*sigma;
                rc2=rc*rc;
                if(r2>=rc2){
                    fx=0.0;
                    fy=0.0;
                }else{
                    // fx=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                    // fy=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij;
                    fx=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*rxij;
                    fy=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*ryij;
                    if (abs(fx)>100.0){
                        cout<<"upper_right_corner_force:fx"<<fx<<endl;
                        cout<<"rxij"<<rxij<<endl;
                    }
                    if (abs(fy)>100.0){
                        cout<<"upper_right_corner_force:fx"<<fx<<endl;
                        cout<<"rxij"<<rxij<<endl;
                    }
                }
                ax[i]+=fx;
                ay[i]+=fy;
                
            }

        }   
    }
}





void bottom_right_corner_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter){

    double rxij;
    double ryij;
    double r2;
    double ir,ir3;
    double ir2,ir6,ir12;
    double sigma;
    double sigma3;
    double sigma6,sigma12;
    double rc,rc2;
    double fx,fy;

    double bottom_right_corner_x = 0.0;
    double bottom_right_corner_y = -STATIC_DIA/2.0;
    
    // double rc2=(STATIC_DIA+SMALL_DIA);

    // Obstacleの中心は0,0なのでディスプレースメントは計算しない、してもいいけど

    for (int i=0;i<NP;i++){
        if (rx[i]>=0.0){
            if (ry[i]<= -STATIC_DIA/2.0){
                rxij=rx[i]-bottom_right_corner_x;
                ryij=ry[i]-bottom_right_corner_y;
                r2=rxij*rxij+ryij*ryij;
                ir2=1.0/r2;
                ir=sqrt(ir2);
                ir3=ir*ir*ir;
                ir6=ir2*ir2*ir2;
                ir12=ir6*ir6;
                sigma=diameter[i]/2.0;
                sigma3=sigma*sigma*sigma;

                sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
                sigma12=sigma6*sigma6;   
                rc=pow(2.0,1.0/6.0)*sigma;
                rc2=rc*rc;
                if(r2>=rc2){
                    fx=0.0;
                    fy=0.0;
                }else{
                    fx=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*rxij;
                    fy=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*ryij;
                    if (abs(fx)>100.0){
                        cout<<"bottom_right_corner_force:fx"<<fx<<endl;
                        cout<<"rxij"<<rxij<<endl;
                    }
                    if (abs(fy)>100.0){
                        cout<<"bottom_right_corner_force:fx"<<fx<<endl;
                        cout<<"rxij"<<rxij<<endl;
                    }
                }
                ax[i]+=fx;
                ay[i]+=fy;
                
            }
        }   
    }
}




void plane_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter){

    double rxij;
    double r2;
    double ir,ir3;
    double ir2,ir6,ir12;
    double sigma;
    double sigma3;
    double sigma6,sigma12;
    double rc,rc2;
    double fx;
    
    for (int i=0;i<NP;i++){
        if (rx[i]>=0.0 ){
            if (ry[i]>-STATIC_DIA/2.0 && ry[i] < STATIC_DIA/2.0){
                rxij=rx[i]-0.0;
                r2=rxij*rxij;
                ir2=1.0/r2;
                ir=sqrt(ir2);
                ir3=ir*ir*ir;
                ir6=ir2*ir2*ir2;
                ir12=ir6*ir6;
                sigma=diameter[i]/2.0;
                sigma3=sigma*sigma*sigma;
                sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
                sigma12=sigma6*sigma6;   
                rc=pow(2.0,1.0/6.0)*sigma;
                rc2=rc*rc;
                if(r2>=rc2){
                    fx=0.0;
                }else{
                    fx=(12.0*sigma12*ir12-18*sigma6*ir6+pow(2,1.5)*3*sigma3*ir3)*ir2*rxij;
                    if (abs(fx)>100.0){
                        cout<<rx[i]<<endl;
                        cout<<ry[i]<<endl;
                        cout<<" plane:fx"<<fx<<endl;
                    }
                }

                ax[i]+=fx;

            }
        }   
    }
}

// 論文と一緒のポテンシャルにしてみたｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍｍ




void gmap_create4(int NP,vector<vector<int>> &g_map,
vector<vector<int>> &pairlist,
vector<double> rx,
vector<double> ry,
vector<int> mask_row,
vector<int> mask_col){

	int gx_map;
	int gy_map;
    int select_index;
    g_map.assign(N_GY,vector<int> (N_GX,-1));
    pairlist.assign(NP,vector<int>(PAIR_LENGTH,-1));
/////////////////////////////////////////////////
    for (int i=0;i<NP;i++){
        gx_map=int((rx[i]+LX/2.0)/LEN_BOX_X);
        gy_map=int((ry[i]+LY/2.0)/LEN_BOX_Y);
        if (gx_map>=N_GX){
            cout<<"gx_map"<<gx_map<<endl;
            cout<<"rx[i]"<<rx[i]<<endl;
            cout<<i<<endl;
        }
        if (gx_map<0){
            cout<<"gx_map"<<gx_map<<endl;
            cout<<"rx[i]"<<rx[i]<<endl;
            cout<<i<<endl;


        }
        if(gy_map>=N_GY){
            cout<<"gy_map"<<gy_map<<endl;
            cout<<"ry[i]"<<ry[i]<<endl;
            cout<<i<<endl;

        }
        if (gy_map<0){
            cout<<"gy_map"<<gy_map<<endl;
            cout<<"ry[i]"<<ry[i]<<endl;
            cout<<i<<endl;


        }
        g_map.at(gy_map).at(gx_map)=i;
    }
    for (int i=0;i<N_GY;i++){
        for(int j =0;j<N_GX;j++){
            select_index=g_map.at(i).at(j);
            if(select_index!=-1){
                make_pairlist(mask_col,mask_row,j,i,g_map,pairlist);
            }
        }
    }
}



void make_pairlist(vector<int> mask_col,vector<int> mask_row,
int select_gx,int select_gy,
vector<vector<int>> G_MAP,
vector<vector<int>> &pairlist){

    int partcle_counter=0;
    int search_gx;
    int search_gy;
    int search_index;
    int select_index=G_MAP.at(select_gy).at(select_gx);

    for(int k=0;k<MASK_LENGTH;k++){
        search_gx=select_gx+mask_col[k];
        search_gy=select_gy+mask_row[k];
        //  グリッドの周期委境界
        if(search_gx>=N_GX){
            search_gx-=N_GX;
        }
        else if (search_gx<0){
            search_gx+=N_GX;
        }
        if (search_gy>=N_GY){
            search_gy-=N_GY;
        }
        else if(search_gy<0){
            search_gy+=N_GY;
        }

        search_index=G_MAP.at(search_gy).at(search_gx);
        if(search_index!=-1){
            if (partcle_counter>=PAIR_LENGTH-1){
                cout<<"pairlist over flow"<<endl;
                exit(1);
            }
            pairlist[select_index][partcle_counter]=search_index;
            partcle_counter+=1;
        }
    }
    pairlist[select_index][PAIR_LENGTH-1]=partcle_counter;
}



void initialize_r(vector<double> &rx,vector<double> &ry,double TEMP_RHO){

    const double dx=sqrt(1.0/TEMP_RHO);
    const double dy=dx;
    double x;
    double y;
    double r2;
    const int Nx=int(LX/dx);
    const int Ny=int(LY/dy);

    const int remove_dx_num=ceil(STATIC_DIA*0.5/dx);
    for(int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){

            x=i*dx+dx/2-0.5*LX;
            y=j*dy+dy/2-0.5*LY;
            if (x<0){
                r2=x*x+y*y;
                if(r2<(0.5*STATIC_DIA+1.0)*(0.5*STATIC_DIA+1.0)){
                    continue;
                }
            }else {
                if (y<STATIC_DIA/2.0+2.0 && y>-STATIC_DIA/2.0-2.0){
                    if(x<1.0){
                        continue;
                    }

                }

            }            
            rx.push_back(x);
            ry.push_back(y);
        }
    }
}
void initialize_va(vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy,int NP){
    double vx_sum=0.0;
    double vy_sum=0.0;
    double kin0=0.0;
    int ave_vx,ave_vy;

    for(int i=0;i<NP;i++){
        vx.push_back(Maxwell_velocity(KBT,M));
        vy.push_back(Maxwell_velocity(KBT,M));
        ax.push_back(0.0);
        ay.push_back(0.0);
    }

    for (int i=0;i<vx.size();i++){
        
        vx_sum=vx_sum+vx[i];
        vy_sum=vy_sum+vy[i];

    }
    ave_vx=vx_sum/(vx.size());
    ave_vy=vy_sum/(vy.size());
    for (int i=0;i<vx.size();i++){
        vx[i]=vx[i]-ave_vx;
        // vx[i]=vx[i]+ave_flow;
        vy[i]=vy[i]-ave_vy;
    }

}

double  Maxwell_velocity(double temp,double m){
    double rand1=genrand_real1();
    double rand2=genrand_real1();
    double velocity;

    velocity=(sqrt(-2*(temp/m)*log(rand1)))*cos(2*PI*rand2);
    return velocity;
}


void make_mask(int mask_size, vector<int> &neighbor_row, vector<int> &neighbor_col) {
    //xの近接リスト作成
    int neighborRowValue = -mask_size;
    int countIndexRow = 0;
    for (int i = 0; i < mask_size + 1; i++) {
        for (int j = 0; j < mask_size; j ++) {
            neighbor_row[countIndexRow] = neighborRowValue;
            countIndexRow++;
        }
        neighborRowValue++;
    }
    for (int i = 0; i < mask_size; i++) {
        for (int j = 0; j < mask_size + 1; j ++) {
            neighbor_row[countIndexRow] = neighborRowValue;
            countIndexRow++;
        }
        neighborRowValue++;
    }
    //yの近接リスト作成
    int neighborColValue = -mask_size;
    int countIndexCol = 0;
    for (int i = 0; i < mask_size + 1; i++) {
        for (int j = 0; j < mask_size; j ++) {
            neighbor_col[countIndexCol] = neighborColValue;
            countIndexCol++;
            neighborColValue++;
        }
        neighborColValue = -mask_size;
    }
    for (int i = 0; i < mask_size; i++) {
        for (int j = 0; j < mask_size + 1; j ++) {
            neighbor_col[countIndexCol] = neighborColValue;
            countIndexCol++;
            neighborColValue++;
        }
        neighborColValue = -mask_size;
    }
}


