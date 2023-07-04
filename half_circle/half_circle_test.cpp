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

using namespace std;
#define PI 3.14159265359


// 関数宣言////////////////////////////////////////////////////////


void make_mask(int mask_size, vector<int> &neighbor_row, vector<int> &neighbor_col);
double  Maxwell_velocity(double temp,double m);

void initialize_va(vector<double> &ax,vector<double> &ay,
vector<double> &vx,vector<double> &vy,int NP,
double KBT,double M);

void initialize_r(vector<double> &rx,vector<double> &ry,
double TEMP_RHO,double LX,double LY,double STATIC_DIA);

void make_pairlist(vector<int> mask_col,vector<int> mask_row,
int select_gx,int select_gy,
vector<vector<int>> G_MAP,
vector<vector<int>> &move_PAIRLIST,
int MASK_LENGTH,int N_GX,int N_GY,int PAIR_LENGTH
);




void gmap_create4(int NP,vector<vector<int>> &g_map,
vector<vector<int>> &pairlist,
vector<double> rx,
vector<double> ry,
vector<int> mask_row,
vector<int> mask_col,
int N_GX,int N_GY,
double LEN_BOX_X,double LEN_BOX_Y,
double PAIR_LENGTH,
int MASK_LENGTH
);

void plane_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter,
double STATIC_DIA,double EPSILON
);

void curve_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter,
double STATIC_DIA,double EPSILON);

void force(int NP,vector<double> &ax,vector<double> &ay,
vector<vector<int>> move_PAIRLIST,
vector<double> rx,vector<double> ry,vector<double> diameter,
int PAIR_LENGTH,double EPSILON,double LX,double LY);

int BookKeeping(int NP,vector<double> Vx,vector<double> Vy,
int MASK_SIZE,
int LEN_BOX_X,double R_CUT,double DT);
///////////////////////////////////////////////////////////////


int main(){


    // int ver=3;
    const double STATIC_DIA=100.0;
    const double TEMP_RHO=0.8;
    const double L=STATIC_DIA*2;
    const double LX=L;
    const double LY=L;
    const double KBT=1.0;
    const double M=1.0;
    const double DT=5e-4;
    const double NSTEPS=1.5e5;
    const double POS_OUT_PERIOD=NSTEPS/10;
    const double t_end=NSTEPS*DT;
    const double h_rev=1/DT ;
    const double h2=0.5*DT*DT;


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
    const double V0=1.0;
    // particle mobility.
    const double MU=1.0;
    // rotational diffusion coefficient.
    const double ETA=1.0;

    const double SMALL_DIA=1.0;
    // const double epsilon=10.0;
    int reusableCount=0;
    // const string ver="T=1_"+to_string(rho).substr(0,4)+"_"+to_string(ave_flow).substr(0,3)+"_"+to_string(static_dia).substr(0,2);
    const string ver="test";
    cout<<ver<<endl;
    
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
    cout<<"1"<<endl;

    make_mask(MASK_LENGTH, mask_row, mask_col);
    cout<<"2"<<endl;
    
   

   
    cout<<"DFdfd"<<endl;
    static vector<vector<int>> g_map(N_GY,vector<int>(N_GX));
    cout<<"3"<<endl;
    vector<vector<int>> move_pairlist;

    cout<<"Dfd"<<endl;
    const string main_dir="./"+ver;
    const char *dir_name=main_dir.c_str();
    //mkdir(dir_name);
    mkdir(dir_name,0755);
    const string rx_file=main_dir+"/rx.dat";
    const string ry_file=main_dir+"/ry.dat";
    const string vx_file=main_dir+"/vx.dat";
    const string vy_file=main_dir+"/vy.dat";
    const string data_file=main_dir+"/lg-data.dat";
    cout<<"Dfd"<<endl;
    

    ofstream rxofs(rx_file);
    ofstream ryofs(ry_file);
    ofstream vxofs(vx_file);
    ofstream vyofs(vy_file);
    ofstream dataofs(data_file);
    time_t start_time,end_time;
    start_time=time(NULL);
   

    initialize_r(rx,ry,TEMP_RHO,LX,LY,STATIC_DIA);
    int N=rx.size();

    initialize_va(WCA_x,WCA_y,vx,vy,N,KBT,M);

    // cout<<Nx*Ny<<endl;

    diameter.assign(N,SMALL_DIA);
    // cout<<N<<endl;
    move_pairlist.resize(N,vector<int>(PAIR_LENGTH));
    
    // theta初期化
    theta.resize(N);
    for (int i=0;i<N;i++){
        theta[i]=2*PI*dist(engine);
    }

    gmap_create4(N,g_map,move_pairlist,rx,ry,mask_row,mask_col,
    N_GX,N_GY,LEN_BOX_X,LEN_BOX_Y,PAIR_LENGTH,MASK_LENGTH);



    force(N,WCA_x,WCA_y,move_pairlist,rx,ry,diameter,PAIR_LENGTH,EPSILON,LX,LY);
    curve_force(WCA_x,WCA_y,rx,ry,N,diameter,STATIC_DIA,EPSILON);
    plane_force(WCA_x,WCA_y,rx,ry,N,diameter,STATIC_DIA,EPSILON);



    // メインループ
    for(int t=0;t<NSTEPS;t++){
        cout<<t<<endl;
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
            if(rx[i]>LX){
                rx[i]-=LX;
            }
            else if(rx[i]<=0.0){
                rx[i]+=LX;
            }

            if(ry[i]>+LY){
                ry[i]-=LY;
            }
            else if(ry[i]<=0.0){
                ry[i]+=LY;
            }
        }
        
        // BookKeeping./////////////////////////////////////////
        if(reusableCount==0){
            gmap_create4(N,g_map,move_pairlist,rx,ry,mask_row,mask_col,
            N_GX,N_GY,LEN_BOX_X,LEN_BOX_Y,PAIR_LENGTH,MASK_LENGTH);
            reusableCount=BookKeeping(N,vx,vy,MASK_SIZE,LEN_BOX_X, R_CUT, DT);    
        }
        else{
            reusableCount-=1;
        }



        WCA_x.assign(N,0.0);
        WCA_y.assign(N,0.0);
        force(N,WCA_x,WCA_y,move_pairlist,rx,ry,diameter,PAIR_LENGTH,EPSILON,LX,LY);
        curve_force(WCA_x,WCA_y,rx,ry,N,diameter,STATIC_DIA,EPSILON);
        plane_force(WCA_x,WCA_y,rx,ry,N,diameter,STATIC_DIA,EPSILON);


// end_BookKeeping./////////////////////////////////////////

        if((t>NSTEPS*0.8)&&(t%200==0)){
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






int BookKeeping(int NP,vector<double> Vx,vector<double> Vy,
int MASK_SIZE,
int LEN_BOX_X,double R_CUT,double DT
){
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
vector<vector<int>> move_PAIRLIST,
vector<double> rx,vector<double> ry,vector<double> diameter,
int PAIR_LENGTH,double EPSILON,double LX,double LY
)
{
    int roop_num;
    int pair_index;
    double rxij;
    double ryij;
    double r2;
    double ir2,ir6,ir12;
    double sigma;
    double sigma6,sigma12;
    double rc,rc2;

    double fx,fy;
    

    for (int i=0;i<NP;i++){
        roop_num=move_PAIRLIST[i][PAIR_LENGTH-1];
        for (int j=0;j<roop_num;j++){
            pair_index=move_PAIRLIST[i][j];
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
                fx=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                fy=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij;
            }
            ax[i]+=fx;
            ay[i]+=fy;
            ax[pair_index]-=fx;
            ay[pair_index]-=fy;
        }
    }
}



void curve_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter,
double STATIC_DIA,double EPSILON
){

    double rxij;
    double ryij;
    double r2;
    double ir2,ir6,ir12;
    double sigma;
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
            ir6=ir2*ir2*ir2;
            ir12=ir6*ir6;
            sigma=(diameter[i]+STATIC_DIA)/2.0;
            sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
            sigma12=sigma6*sigma6;   
            rc=pow(2.0,1.0/6.0)*sigma;
            rc2=rc*rc;
            if(r2>=rc2){
                fx=0.0;
                fy=0.0;
            }else{
                fx=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                fy=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*ryij;
            }
            ax[i]+=fx;
            ay[i]+=fy;
        }   
    }
}



void plane_force(vector<double> &ax,vector<double> &ay,
vector<double> rx, vector<double> ry,int NP,vector<double> diameter,
double STATIC_DIA,double EPSILON
){

    double rxij;
    double r2;
    double ir2,ir6,ir12;
    double sigma;
    double sigma6,sigma12;
    double rc,rc2;
    double fx;
    
    for (int i=0;i<NP;i++){
        if (rx[i]>=0.0 ){
            if (ry[i]>-STATIC_DIA/2.0 && ry[i] < STATIC_DIA/2.0){
                rxij=rx[i]-0.0;
                r2=rxij*rxij;
                ir2=1.0/r2;
                ir6=ir2*ir2*ir2;
                ir12=ir6*ir6;
                sigma=diameter[i]/2.0;
                sigma6=sigma*sigma*sigma*sigma*sigma*sigma;
                sigma12=sigma6*sigma6;   
                rc=pow(2.0,1.0/6.0)*sigma;
                rc2=rc*rc;
                if(r2>=rc2){
                    fx=0.0;
                }else{
                    fx=24.0*EPSILON*(2.0*sigma12*ir12-sigma6*ir6)*ir2*rxij;
                }
                ax[i]+=fx;
            }
        }   
    }
}





void gmap_create4(int NP,vector<vector<int>> &g_map,
vector<vector<int>> &pairlist,
vector<double> rx,
vector<double> ry,
vector<int> mask_row,
vector<int> mask_col,
int N_GX,int N_GY,
double LEN_BOX_X,double LEN_BOX_Y,
double PAIR_LENGTH,
int MASK_LENGTH
){

	int gx_map;
	int gy_map;
    int select_index;
    g_map.assign(N_GY,vector<int> (N_GX,-1));
    pairlist.assign(NP,vector<int>(PAIR_LENGTH,-1));

/////////////////////////////////////////////////
    for (int i=0;i<NP;i++){
        gx_map=int(rx[i]/LEN_BOX_X);
        gy_map=int(ry[i]/LEN_BOX_Y);
        g_map.at(gy_map).at(gx_map)=i;

    }

    for (int i=0;i<N_GY;i++){
        for(int j =0;j<N_GX;j++){
            select_index=g_map.at(i).at(j);
            if(select_index!=-1){
                make_pairlist(mask_col,mask_row,j,i,g_map,pairlist,MASK_LENGTH,N_GX,N_GY,PAIR_LENGTH);
            }
        }
    }

}



void make_pairlist(vector<int> mask_col,vector<int> mask_row,
int select_gx,int select_gy,
vector<vector<int>> G_MAP,
vector<vector<int>> &move_PAIRLIST,
int MASK_LENGTH,int N_GX,int N_GY,int PAIR_LENGTH
){

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
            move_PAIRLIST[select_index][partcle_counter]=search_index;
            partcle_counter+=1;
        }
    }
    move_PAIRLIST[select_index][PAIR_LENGTH-1]=partcle_counter;
}



void initialize_r(vector<double> &rx,vector<double> &ry,
double TEMP_RHO,double LX,double LY,double STATIC_DIA
){

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
                if(r2<(0.5*STATIC_DIA+0.5)*(0.5*STATIC_DIA+0.5)){
                    continue;
                }
            }else {
                if (y<STATIC_DIA/2.0 && y>-STATIC_DIA/2.0){
                    if(x<0.5){
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
vector<double> &vx,vector<double> &vy,int NP,
double KBT,double M
){
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


