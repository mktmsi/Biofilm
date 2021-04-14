#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <GLUT/glut.h>
#include <time.h>
#include "Monitor.hpp"
//#include <opencv/cv.h>
#include "Vector.h"
#include <string.h>
//#include "DefineCell.h"
#include <unistd.h>
#include <sys/stat.h>
#include <png.h>
//#include <windows.h>

#define dt (0.001) //0.001              //タイムステップ
double dx = 0.1, dxdx = dx * dx; //タイムステップ
#define MAX_STEP (5000000)       //最大タイムステップ数50000000
#define PI (3.14159265359)
#define N (101) //一辺のマス数,N*N＝セル数
//#define N0 (100)

//ディスプレイ
#define csize (2.5)
double circlesize;
int zoom =5; //数値大→遠くなる N=100で5 N=1000で28
#define interval (7.5)
double xcenter=10.0; //大→左へずれる（フィールドが） N=1000で100 N=100で10
double ycenter=7.0;  //大→下へずれる
//int cell_max=N*N;

Monitor monitor;

//要素
double u_next[N], v_next[N], z_next[N], memoryF[N], w_next[N];
double u_old[N], v_old[N], z_old[N], w_old[N];
double LapV[N], LapU[N], LapZ[N], LapW[N];
int ActiveFlag[N];
/*先生に送ったやつ
//初期値
const double u0 = (30.0) * 10.0, v0 = 100.00, z0 = 0.00, w0 = 0.0; //100.000*0.30000:anti      //u0 = 1.000, v0 = 1.00, z0 = 0.00;
double threshV = 50.00;                                      //v0 / 2.0;
const double beta = (21.00000000) * 10.0, alpha = beta + 00.0000;  //閾値 beta = 21.00000000, alpha = beta - 0.0000;
int N_initTadius_v = 10;
double rawLapU, rawU_sum, rawLapU_max;
//初期値
const double u_in = u0, w_in = 0.00;                                                   //   u_in, uin_haji = 0.100, w_in = 00.00;
const double k_1 = 1.0, k_2 = 19.0, kh = 100.0, k_vz = 5.5; //double k_1 = 3000.0* (1.0 / 3.0) * 1.01, k_2 = 1.90, k_vz = 50.5000, kh = 3000.0;
const double kz = 0.5, ku = 0.001, kv = 0.0001, kw = 0.001;                            //double kz = 0.020, ku = 0.001, kv = 0.0001, kw = 0.001;
const double kf = 0.01;                                                                //, kf_u = (kf * 0.1);                                                                              //    kf = 0.0005, kf_u = kf * 0.1;/
const double Du = 2.33, Dz = 0.23, Dv = 0.000007, Dw = 3.00 / 3.0;
const double k3 = 0.0, k4 = 0.0;
*/

const double u0 = (30.0) * 10.0, v0 = 100.00, z0 = 0.00, w0 = 0.0; //100.000*0.30000:anti      //u0 = 1.000, v0 = 1.00, z0 = 0.00;
double threshV = 50.00;                                            //v0 / 2.0;

int N_initTadius_v = N * 0.10;
double rawLapU, rawU_sum, rawLapU_max;

const double ku = 0.001;
const double kv = 0.0001; //0.01
const double kz = 0.5;
const double kf = 0.01;
const double k_vz = 5.5; //5.5
const double kh = 100.0; //300.0
const double k_1 = 1.0;
const double k_2 = 19.0;
const double Du = 2.33;     //0.25
const double Dv = 0.000007; //0.01
const double Dz = 0.23;     //0.3
const double u_in = 300.0;  //100.0
const double alpha = 210.0; //2.1
const double beta = 210.0;

int cellcounts = 0;
const double w_in = 0.00, k3 = 0.0, k4 = 0.0, Dw = 0.0;

//画面書き出し用の関数
//double aw = 0.1, au = (0.1) , av = 0.1, az = 0.50, a_interval = 0.2500; //N=1000
double au = (0.02), av = 0.001, az = 3.0, aw = 10.0, a_interval = 0.2500; //N=100
//double au = 2.0, av = 1.60, az = 10.0, a_interval = 0.2500; //original/
int DisplayInterval = 1000;
int DataSaveInterval = DisplayInterval;
int tsJump = 00000;

double SumVLeft = 0.0, SumVRight = 0.0;

int flag_center = 1;
int flag_sigm = 0;
int flag_monitor = 1, flag = 0, flag_variable = 0;
int flag_drug = 0;
int save_flag = 1; //0:画像保存しない，1:画像保存する
int Flaaaag = 0, FlagDataSave = 1;

int count1 = 0, count2 = 0, count3 = 0;
int ParameterSetCount = 0, EndSetCount = 0;
//int save_csv = 0;q
int SignalingFlag[N] = {};
int Signaling = 0;

void capture(int *);

double t;
int winid;
int ts = 0; //現在のタイムステップ
int initial = 0;
FILE *fpValues, *fpV, *fpU, *fp, *fpZ;
FILE *fpSumV;
char filename[256] = "", strFolderParaset[150] = "", strParaSetTxt[256] = "";
char strFolderMovie[256] = "", strMovieFiles[256] = "", filenameSumV[256] = "";
char filenameV[256] = "", filenameU[256] = "", filenameZ[256] = "";
void mouse(int, int, int, int);
void resize(int, int);
void FolderAndFileInit();
void init();
void FieldCycle();
void keyboard(unsigned char, int, int);
void idle();
void display();

void FolderAndFileInit()
{
  char tmp[256] = "";
  char *p;
  int i = 0;
  p = getcwd(tmp, sizeof(tmp) - 1);
  if (p == NULL)
  {
    printf("ERROR,getcwd() ret=NULL\n");
    perror("getcwd\n");
    exit(EXIT_FAILURE);
  }
  /*
  sprintf(strFolderMovie, "./MovieDir/ParameterSet%d",ParameterSetCount);
  if (mkdir(strFolderMovie, S_IRWXU) == 0)
  { //フォルダ作成
    printf("フォルダ作成に成功しました。\n");
  }
  else
  {
    printf("フォルダ作成に失敗しました。\n");
    printf("%s\n", strFolderParaset);
    exit(0);
  }
  sprintf(strMovieFiles, "%s/temp_",strFolderMovie);
  monitor.SetMovieName(strMovieFiles);
*/
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  if (mkdir(strFolderParaset, S_IRWXU) == 0)
  { //フォルダ作成
    printf("フォルダ作成に成功しました。\n");
  }
  else
  {
    //printf("フォルダ作成に失敗しました。\n");
    do
    {
      ParameterSetCount++;
      sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
    } while (mkdir(strFolderParaset, S_IRWXU) != 0);
    printf("フォルダ作成に成功しました。\n");
    //printf("%s\n", strFolderParaset);
    //exit(0);
  }
  //sprintf(filenamev, "Values.csv");
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(strParaSetTxt, "%s/ParaSet.txt", strFolderParaset);
  fpValues = fopen(strParaSetTxt, "w");
  fprintf(fpValues, "Du=%lf,Dv=%lf,Dz=%lf,u0=%lf,v0=%lf,z0=%lf,beta=%lf,alpha=%lf,NinitTadius_v=%d,u_in=%lf,k_1=%lf,k_2=%lf,k_vz=%lf,kh=%lf,kz=%lf,ku=%lf,kv=%lf,kf=%lf\n", Du, Dv, Dz, u0, v0, z0, beta, alpha, N_initTadius_v, u_in, k_1, k_2, k_vz, kh, kz, ku, kv, kf);
  fclose(fpValues);
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(filename, "%s/ParameterSet%d.csv", strFolderParaset, ParameterSetCount);
  fp = fopen(filename, "w");
  if (fpValues == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filename=%s\n", filename);
    exit(0);
  }
  fprintf(fp, "ts,u,v,z,signal\n");
  fprintf(fp, "%d,%f,%f,%f,%d\n", ts, u_next[N / 2], v_next[N / 2], z_next[N / 2], SignalingFlag[N / 2]);
  fclose(fp);
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(filenameV, "%s/V_para.csv", strFolderParaset);
  fpV = fopen(filenameV, "w");
  if (fpV == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filenameV=%s\n", filenameV);
    exit(0);
  }
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(filenameU, "%s/U_para.csv", strFolderParaset);
  fpU = fopen(filenameU, "w");
  if (fpU == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filenameU=%s\n", filenameU);
    exit(0);
  }
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(filenameZ, "%s/Z_para.csv", strFolderParaset);
  fpZ = fopen(filenameZ, "w");
  if (fpV == NULL)
  {
    printf("ファイル作成失敗fpZ\n");
    printf("filenameZ=%s\n", filenameZ);
    exit(0);
  }

  fprintf(fpV, "%d,", ts);
  fprintf(fpU, "%d,", ts);
  fprintf(fpZ, "%d,", ts);
  for (i = 0; i < N; i++)
  {
    fprintf(fpV, "%lf,", v_next[i]);
    fprintf(fpU, "%lf,", u_next[i]);
    fprintf(fpZ, "%lf,", z_next[i]);
  }
  fprintf(fpV, "\n");
  fclose(fpV);
  fprintf(fpU, "\n");
  fclose(fpU);
  fprintf(fpZ, "\n");
  fclose(fpZ);
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(filenameSumV, "%s/SumV.csv", strFolderParaset);
  fpSumV = fopen(filenameSumV, "w");
  if (fpSumV == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filename=%s\n", filenameSumV);
    exit(0);
  }
  //fprintf(fpCellcount, "ts,left,right\n");
  //fprintf(fpSumV, "%d,%lf,%lf\n", ts, SumVLeft,SumVRight);
  fclose(fpSumV);
  SumVLeft = 0.0, SumVRight = 0.0;
}

void init()
{


  
  if (N / 100 == 1)
  {
    dx = 0.1;
    xcenter=10.0;
    zoom=5;
au = (0.02), av = 0.001, az = 3.0, aw = 10.0, a_interval = 0.2500;
    dxdx = dx * dx;
  }
  else if (N / 100 == 10)
  {
    zoom=28;
    xcenter=100.0;
    
    aw = 0.1, au = (0.1) , av = 0.1, az = 0.50, a_interval = 0.2500; //N=1000
    dx = 0.01;
    dxdx = dx * dx;
  }
  monitor.SetZoom(zoom);
  monitor.SetCenter(xcenter, ycenter);

  /*
  if (dt > (dxdx / (2 * Dv)))
  {
    printf("error! Not stable> Dv\n");
    exit(1);
  }
  else if (dt > (dxdx / (2 * Du)))
  {
    printf("error! Not stable> Du\n");
    exit(1);
  }
  else if (dt > (dxdx / (2 * Dz)))
  {
    printf("error! Not stable> Dz\n");
    exit(1);
  }
  else if (dt > (dxdx / (2 * Dw)))
  {
    printf("error! Not stable> Dw\n");
    exit(1);
  }
  */
  //変数の無次元化

  int x, i, cell_count;
  //glCaptureの設定
  // glCapture.setWriteFile("output.avi");
  /*
  if (Flaaaag == 0)
  {
    printf("Start ParameterSetcount?\n");
    scanf("%d", &ParameterSetCount);
    printf("End ParameterSetcount?\n");
    scanf("%d", &EndSetCount);
    count1 = ParameterSetCount % 3;
    count2 = ParameterSetCount / 3;
    count3 = ParameterSetCount / 9;
  }
  
  if (FlagDataSave == 1&&Flaaaag)
  {
    if (count1 % 2 == 0 && count1 != 0)
    {
      count1 = 0;
      count2++;
    }
    else 
    {
      count1++;
    }

    if (count2 % 3 == 0 && count2 != 0)
    {
      count2 = 0;
      count3++;
    }
  }
    Flaaaag = 1;
    printf("ParameterSetCount=%d\n", ParameterSetCount);
    printf("count1=%d count2=%d count3=%d\n", count1, count2,count3);
  //k2 = (double)(count1)*1.0;
  //ks = 0.001 * pow(10.0, (double)count2);
  Du = 0.0001 * pow(10.0, (double)count1);
  Dv = 0.000001 * pow(10.0, (double)count2);
  Dz = 0.01 * pow(10.0, (double)count3);
*/
  ts = 0;
  t = 0.0;
  //DataCount = 0;

//乱数初期化
srand((unsigned int)time(NULL));
  //initialize Field

  for (x = 0; x < N; x++)
  {
    z_next[x] = 0.00000000+(rand()*(1.0+1.0)/(1.0+RAND_MAX));
    z_old[x] = z_next[x];
    u_next[x] = u0+(rand()*(1.0+1.0)/(1.0+RAND_MAX));;
    u_old[x] = u_next[x];
    w_next[x] = 0.0;
    w_old[x] = w_next[x];

    SignalingFlag[x] = 0;
    if ((0 <= x && x <= N_initTadius_v) || (N - 1 - N_initTadius_v <= x && x <= N - 1))
    { //初期バイオフィルムの範囲
      //  u(x,y)=(ku0+1.500)+2.0*((ku0+1.500)*0.01)*(double)rand()/RAND_MAX-((ku0+1.500)*0.01);
      v_next[x] = v0+(rand()*(1.0+1.0)/(1.0+RAND_MAX));
      v_old[x] = v_next[x];
      if (x < N / 2)
      {
        //cellcountLeft++;
        SumVLeft += v_next[x];
      }
      else
      {
        //cellcountRight++;
        SumVRight += v_next[x];
      }
    }
    else
    {
      v_next[x] = 0.0;
      v_old[x] = v_next[x];
    }
  }
  if (FlagDataSave == 1)
  {
    FolderAndFileInit();
  }
}

void FieldCycle()
{
  int i = 0;
  int x;
  double u_sum = 0.0, v_sum = 0.0, z_sum = 0.0, w_sum = 0.0;
  double f, h, ax;
  //zmin = z_next[o];
  for (i = 0; i < DisplayInterval; i++)
  {
    SumVRight = 0.0, SumVLeft = 0.0;

    for (x = 0; x < N; x++)
    {
      u_sum = 0.0, v_sum = 0.0, z_sum = 0.0, w_sum = 0.0;
      //拡散項を求める
      if (x != 0)
      {
        v_sum += v_old[x - 1];
        z_sum += z_old[x - 1];
        u_sum += u_old[x - 1];
      }
      else
      {
        v_sum += v_old[x];
        z_sum += z_old[x];
        u_sum += u_old[x];
      }

      if (x != N - 1)
      {
        v_sum += v_old[x + 1];
        z_sum += z_old[x + 1];
        u_sum += u_old[x + 1];
      }
      else
      {
        u_sum += u_old[x];
        v_sum += v_old[x];
        z_sum += z_old[x];
      }
      LapU[x] = 2.0 * (u_sum / 2.0 - u_old[x]) / dxdx;
      LapV[x] = 2.0 * (v_sum / 2.0 - v_old[x]) / dxdx;
      LapZ[x] = 2.0 * (z_sum / 2.0 - z_old[x]) / dxdx;

      //方程式を解く
      cellcounts = 0;

      ax = k_1 * (u_old[x] - k_2 * z_old[x] - beta);

      f = v_old[x] * (1.0 / (1.0 + exp(-ax))); //膜電位・食欲

      SignalingFlag[x] = 0;
      h = (1.0 - 1.0 / (1.0 + exp(-kh * (u_old[x] - alpha)))); //シグナル
      if (h > 0.00001)
      {
        SignalingFlag[x] = 1;
      }
      memoryF[x] = f;
      z_next[x] = z_old[x] + (-kz * z_old[x] + Dz * LapZ[x] + k_vz * h * v_old[x]) * dt;

      v_next[x] = v_old[x] + (-kv * v_old[x] + Dv * LapV[x] + kf * f) * dt;

      u_next[x] = u_old[x] + (-ku * u_old[x] + Du * LapU[x] - kf * f) * dt;

      if (x == (N) / 2) //int同士は切り捨て
      {
        u_next[x] += u_in * dt;
      }

      if (x < N / 2)
      {
        SumVLeft += v_next[x];
      }
      else
      {
        SumVRight += v_next[x];
      }
    }
    for (x = 0; x < N; x++)
    {
      u_old[x] = u_next[x];
      v_old[x] = v_next[x];
      z_old[x] = z_next[x];
    }
    ts++;
    t = ts * dt;
  }
}

void keyboard(unsigned char key, int x, int y)
{
  double tmpcin;
  switch (key)
  {
  case 'i':
    if (flag_variable == 0)
    {
      au = au + a_interval;
    }
    else if (flag_variable == 1)
    {
      av = av + a_interval;
    }
    else
    {
      az = az + a_interval;
    }
    break;
  case 'd':
    if (flag_variable == 0)
    {
      au = au - a_interval;
    }
    else if (flag_variable == 1)
    {
      av = av - a_interval;
    }
    else
    {
      az = az - a_interval;
    }
    break;
  case 'f':
    flag_drug = 1 - flag_drug;
    break;
  case 'q':
    exit(0);
    break;
  case '\033': /* '\033' = ESC */
    exit(0);
    break;
  case 'h':
    monitor.SetCenter(1.0, 0);
    break;
  case 'l':
    monitor.SetCenter(-1.0, 0);
    break;
  case 'j':
    monitor.SetCenter(0, 10.0);
    break;
  case 'k':
    monitor.SetCenter(0, -10.0);
    break;
  case 'c':
    flag_variable += 1;
    if (flag_variable > 2)
    {
      flag_variable = 0;
    }
    break;
  case 'w':

    break;
  case 'z':
    monitor.SetZoom(1.1 / 1.0);
    break;
  case 'x':
    monitor.SetZoom(1.0 / 1.1);
    break;
  }
}

void idle(void)
{
  glutSetWindow(winid);
  glutPostRedisplay();
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  int i, j;
  //double r,g,b,y=0.95,x=0.68,width=500.0,gap=560.0,y_vBottom=1890.0,y_vUp=y_vBottom+width,y_zUp=y_vBottom-gap,y_zBottom=y_zUp-width,y_uUp=y_zBottom-gap,y_uBottom=y_uUp-width;
  double r, g, b, y = 0.85, x = 0.68, width = 500.0, gap = 700.0, y_vBottom = 1690.0, y_wBottom = y_vBottom - gap - width, y_uBottom = y_wBottom - gap - width;
  double f;
  static double difference;
  static double oldu, oldv;
  //ディスプレイ表示用文字列
  char str[20];

  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  FieldCycle(); //シミュレーションを回す

  //yup=(interval*(N-1)+interval/2-N)/2.0;
  if (ts == 0 || flag_monitor == 1 && ts >= tsJump)
  {
    //printf("v=%.2f u=%.4f z=%.5f w=%.5f\n", v_next[0], u_next[0], z_next[0], w_next[0]);

    for (j = 0; j <= N - 2; j++)
    {
      monitor.SetAllColor(0.0, 0.0, 0.0);
      monitor.DrawLine(0 * 0.2, au * alpha, N * 0.2, au * alpha, 0.10);

      monitor.SetAllColor(0.0, 1.0, 0.0);
      if (alpha > u_next[j])
      {
        monitor.SetAllColor(1.0, 0.0, 0.0);
      }
      monitor.DrawLine(j * 0.2, au * u_next[j], (j + 1) * 0.2, au * u_next[j + 1], 2.0);
      if (memoryF[j] < 0.0001)
      {
        monitor.SetAllColor(1.0, 0.60, 0.0);
      }
      else
      {
        monitor.SetAllColor(0, 0, 1.0);
      }
      monitor.DrawLine(j * 0.2, av * v_next[j], (j + 1) * 0.2, av * v_next[j + 1], 2.0);
      monitor.SetAllColor(1.0, 0.0, 0.0);
      monitor.DrawLine(j * 0.2, az * z_next[j], (j + 1) * 0.2, az * z_next[j + 1], 2.0);
      monitor.SetAllColor(0.0, 0.0, 0.0);
      monitor.DrawLine(j * 0.2, aw * w_next[j], (j + 1) * 0.2, aw * w_next[j + 1], 2.0);
    }
    monitor.SetAllColor(0.0, 0.0, 0.0);
    //monitor.DrawLine(0, 0.0, 0.2 * N, 0.0, 2.0);
    /*
    sprintf(strv, "V");
    sprintf(stru, "U");
    sprintf(strz, "Z");
    if (save_csv == 1)
    {
      fprintf(fpV, "\n");
      fprintf(fpU, "\n");
      fprintf(fpZ, "\n");
    }
    */
    monitor.SetAllColor(0.0, 0.0, 0.0);

    /*
    sprintf(stru0,"u0=%.3f",u0);
    sprintf(strv0,"v0=%.3f",v0);
    sprintf(strz0,"z0=%.3f",z0);
    sprintf(struin,"uin=%.3f",u_in);
    sprintf(strkf,"kf=%.3f",kf);
    sprintf(strkh,"kh=%.3f",kh);
    sprintf(strkvz,"k_vz=%.3f",k_vz);
    sprintf(strkz,"kz=%.3f",kz);
    sprintf(strku,"ku=%.3f",ku);
    sprintf(strkv,"kv=%.3f",kv);
    sprintf(strk2,"k2=%.3f",k_2);
    sprintf(strk3,"k3=%.3f",k_3);
    sprintf(strDu,"Du=%.3f",Du);
    sprintf(strDz,"Dz=%.3f",Dz);
    sprintf(strDv,"Dv=%.3f",Dv);
    //最大値
    sprintf(strumax,"umax=%f",umax);
    sprintf(strvmax,"vmax=%f",vmax);
    sprintf(strzmax,"zmax=%f",zmax);
    sprintf(strumin,"umin=%f",umin);
    sprintf(strvmin,"vmin=%f",vmin);
    sprintf(strzmin,"zmin=%f",zmin);
    sprintf(strLapUmax,"LapUmax=%f",LapUmax);
    sprintf(strLapVmax,"LapVmax=%f",LapVmax);
    sprintf(strLapZmax,"LapZmax=%f",LapZmax);

    sprintf(strau,"au=%f",au);
    sprintf(strav,"av=%f",av);
    sprintf(straz,"az=%f",az);
*/
    /*monitor.String(x-1.4,y-0.5,strv);
    monitor.String(x-1.4,y-0.9,strz);
    monitor.String(x-1.4,y-1.4,stru);*/
    //sprintf(str_threshV, "thresV=%f", threshV);
    sprintf(str, "ts=%d", ts);
    monitor.String(x - 0.5, y, str);
    y = y - 0.2;
    sprintf(str, "N=%d", N);
    monitor.String(x - 0.5, y, str);
    y = y - 0.1;
    //monitor.String(x,y-0.25,str_threshV);y=y-0.05;
    /*if (save_csv == 1)
    {
      fprintf(fp1, "%d,%d\n", ts, cellcounts);
    }*/
    /*
    if(flag_monitor==0){
            monitor.String(x,y,stru);y=y-0.05;
        }else if(flag_monitor==1){
            monitor.String(x,y,strv);y=y-0.05;
        }else{
            monitor.String(x,y,strz);y=y-0.05;
        }
    monitor.String(x,y,stru0);y=y-0.05;
    monitor.String(x,y,strv0);y=y-0.05;
    monitor.String(x,y,strz0);y=y-0.05;
    monitor.String(x,y,struin);y=y-0.05;
    monitor.String(x,y,strkf);y=y-0.05;
    monitor.String(x,y,strkh);y=y-0.05;
    monitor.String(x,y,strkvz);y=y-0.05;
    monitor.String(x,y,strkz);y=y-0.05;
    monitor.String(x,y,strku);y=y-0.05;
    monitor.String(x,y,strkv);y=y-0.05;
    monitor.String(x,y,strk2);y=y-0.05;
    monitor.String(x,y,strk3);y=y-0.05;
    //monitor.String(x,y,strk4);y=y-0.05;
    monitor.String(x,y,strDu);y=y-0.05;
    monitor.String(x,y,strDz);y=y-0.05;
    monitor.String(x,y,strDv);y=y-0.05;

    y=y-0.05;
    monitor.String(x,y,strumax);y=y-0.05;
    monitor.String(x,y,strumin);y=y-0.05;
    monitor.String(x,y,strvmax);y=y-0.05;
    monitor.String(x,y,strvmin);y=y-0.05;
    monitor.String(x,y,strzmax);y=y-0.05;
    monitor.String(x,y,strzmin);y=y-0.05;
    monitor.String(x,y,strLapUmax);y=y-0.05;
    monitor.String(x,y,strLapVmax);y=y-0.05;
    monitor.String(x,y,strLapZmax);y=y-0.05;
    y=y-0.05;
    monitor.String(x,y,strau);y=y-0.05;
    monitor.String(x,y,strav);y=y-0.05;
    monitor.String(x,y,straz);y=y-0.05;
    */
    if (FlagDataSave == 1)
    {
      sprintf(filename, "%s/ParameterSet%d.csv", strFolderParaset, ParameterSetCount);
      fpValues = fopen(filename, "a");
      if (fpValues == NULL)
      {
        printf("ファイル作成失敗fpV\n");
        printf("filename=%s\n", filename);
        exit(0);
      }
      fprintf(fpValues, "%d,%f,%f,%f,%d\n", ts, u_next[N / 2], v_next[N / 2], z_next[N / 2], SignalingFlag[N / 2]);
      fclose(fpValues);
      ////////////SumV用
      fpSumV = fopen(filenameSumV, "a");
      if (fpSumV == NULL)
      {
        printf("ファイル作成失敗fpV\n");
        printf("filename=%s\n", filenameSumV);
        exit(0);
      }
      //fprintf(fpCellcount, "ts,left,right\n");
      fprintf(fpSumV, "%d,%lf,%lf\n", ts, SumVLeft, SumVRight);
      fclose(fpSumV);
      SumVLeft = 0.0, SumVRight = 0.0;
      fpV = fopen(filenameV, "a");
      fpU = fopen(filenameU, "a");
      fpZ = fopen(filenameZ, "a");
      fprintf(fpV, "%d,", ts);
      fprintf(fpU, "%d,", ts);
      fprintf(fpZ, "%d,", ts);
      for (i = 0; i < N; i++)
      {
        fprintf(fpV, "%lf,", v_next[i]);
        fprintf(fpU, "%lf,", u_next[i]);
        fprintf(fpZ, "%lf,", z_next[i]);
      }
      fprintf(fpV, "\n");
      fclose(fpV);
      fprintf(fpU, "\n");
      fclose(fpU);
      fprintf(fpZ, "\n");
      fclose(fpZ);
      //
    }
    capture(&ts);
    //monitor.SavePPMData();
    glFlush();
    glutSwapBuffers();

    if (ts > MAX_STEP)
    {
      if (ParameterSetCount == EndSetCount)
      {
        fclose(fpValues);
        exit(0);
      }
      ParameterSetCount++;
      init(); //初期化＆再スタート
    }
  }
  /*
  if (ts % (DisplayInterval * 10) == 0)
  {
    printf("ts=%d\n", ts);
  }*/
  //ステップごとに画像を保存
  //if (save_flag == 1 && ts % DisplayInterval == 0 || ts == 0)
  // {

  //}

  /*if (end_flag == 1 || ts > MAX_STEP)
  {
    for (i = 0; i < 3; i++)
    {
      monitor.SavePPMData();
    }
    fclose(fp1);
    fclose(fpU);
    fclose(fpV);
    //fclose(fp2);
    exit(0);
  }*/
}

void mouse(int button, int state, int x, int y)
{
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
    }
    else
    {
      glutIdleFunc(idle);
    }
    break;

  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "middle: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "middle: off" << std::endl;
    }
    break;

  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "right: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "right: off" << std::endl;
    }
    break;
  }
}

void resize(int w, int h)
{
  glViewport(0, 0, w, h);

  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-w / 0.0, w / 20.0, -h / 20.0, h / 20.0, -3.0, 3.0);
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}

void OpenGL_init(int *argcp, char **argv)
{
  init(); //初期条件を設定
  glutInit(argcp, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(monitor.GetWindowSize(Monitor::X), monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition(10, 100);
  winid = glutCreateWindow("simulation");
  glutDisplayFunc(display);
  glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
}

void monitor_init()
{
  //   monitor.SetMode( 0 );
  monitor.SetWindowSize(400, 200);
  monitor.SetMovieMode(1);
  //monitor.SetMovieName("./MovieDir/temp_");
  //   monitor.SetGridMode( 0 );
  //   monitor.SetGridWidth( 2.0 );
  //monitor.SetZoom(zoom);
  //monitor.SetCenter(xcenter, ycenter);
  // for( int i = 0 ; i < 22 ; i++ ) monitor.SetZoom( 1.1 );
  // for( int i = 0 ; i < 7; i++ ) monitor.SetCenter( 0 , 1.0 );
  // for( int i = 0 ; i <10 ; i++ ) monitor.SetCenter( 1.0 , 0.0 );
}

void capture(int *pts)
{
  char filepath[100]; //= "./MovieDir/output.png";
  sprintf(filepath, "./MovieDir/%d.png", *pts);
  png_bytep raw1D;
  png_bytepp raw2D;
  int i;
  int width = glutGet(GLUT_WINDOW_WIDTH);
  int height = glutGet(GLUT_WINDOW_HEIGHT);

  // 構造体確保
  FILE *fp = fopen(filepath, "wb");
  png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  png_infop ip = png_create_info_struct(pp);
  // 書き込み準備
  png_init_io(pp, fp);
  png_set_IHDR(pp, ip, width, height,
               8,                   // 8bit以外にするなら変える
               PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // ピクセル領域確保
  raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
  raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
  for (i = 0; i < height; i++)
    raw2D[i] = &raw1D[i * png_get_rowbytes(pp, ip)];
  // 画像のキャプチャ
  glReadBuffer(GL_FRONT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
  glReadPixels(0, 0, width, height,
               GL_RGBA,          // RGBA以外にするなら変える
               GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
               (void *)raw1D);
  // 上下反転
  for (i = 0; i < height / 2; i++)
  {
    png_bytep swp = raw2D[i];
    raw2D[i] = raw2D[height - i - 1];
    raw2D[height - i - 1] = swp;
  }
  // 書き込み
  png_write_info(pp, ip);
  png_write_image(pp, raw2D);
  png_write_end(pp, ip);
  // 開放
  png_destroy_write_struct(&pp, &ip);
  fclose(fp);
  free(raw1D);
  free(raw2D);

  //printf("write out screen capture to '%s'\n", filepath);
}

int main(int argc, char *argv[])
{
  monitor_init();
  std::cout << "monitor init OK" << std::endl;

  OpenGL_init(&argc, argv);
  std::cout << "OpenGL init OK" << std::endl;

  // glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  glutMainLoop(); //無限ループ

  return 0;
}
