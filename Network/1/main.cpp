#include <iostream>
#include <stdio.h>
//#include <math.h>
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

#define dt (0.001)               //0.001              //タイムステップ
double dx = 0.1, dxdx = dx * dx; //タイムステップ
#define MAX_STEP (50000000)      //50000000   //最大タイムステップ数50000000
#define PI (3.14159265359)
#define N (3) //一辺のマス数,N*N＝セル数
//#define N0 (100)

//ディスプレイ
#define csize (2.5)
double circlesize;
int zoom = 2; //数値大→遠くなる N=100で5 N=1000で28
#define interval (7.5)
#define rij(i, j) rij[N * (i) + (j)]
#define dVij(i, j) dVij[N * (i) + (j)]
#define dUij(i, j) dUij[N * (i) + (j)]
#define dZij(i, j) dZij[N * (i) + (j)]
double xcenter = 10.0; //大→左へずれる（フィールドが） N=1000で100 N=100で10
double ycenter = 7.0;  //大→下へずれる
//int cell_max=N*N;

Monitor monitor;

//要素
double u_next[N], v_next[N], z_next[N], memoryF[N], w_next[N];
double u_old[N], v_old[N], z_old[N], w_old[N];
double LapV[N], LapU[N], LapZ[N], LapW[N];
double rij[N * N] = {}, dZij[N * N] = {}, dUij[N * N] = {}, dVij[N * N] = {};
int ActiveFlag[N];

const double v0 = 100.00, z0 = 0.00, w0 = 0.0; //100.000*0.30000:anti      //u0 = 1.000, v0 = 1.00, z0 = 0.00;
 double u0 = 2800.0;
double threshV = 50.00; //v0 / 2.0;

const int N_initTadius_v = 2;
double rawLapU, rawU_sum, rawLapU_max;

const double ku = 0.001*1.0;
const double kv = 0.0001; //0.01
const double kz = 0.5;
const double kf = 0.01*100.0;
const double k_vz = 5.5; //5.5
const double kh = 100.0; //300.0
const double k_1 = 1.0;
const double k_2 = 19.0;
 double Du = 20.33;     //0.25
const double Dv = 0.000001;      //vは拡散しない//0.000007 * 5.0; //0.01
const double Dz = 0.023;     //0.3
 double u_in = u0;  //100.0
const double alpha = 210.0; //2.1
const double beta = 210.0;
 double rijInit=1.0;

int cellcounts = 0;
const double w_in = 0.00, k3 = 0.0, k4 = 0.0, Dw = 0.0;

//画面書き出し用の関数
//double aw = 0.1, au = (0.1) , av = 0.1, az = 0.50, a_interval = 0.2500; //N=1000
double au = (0.01), av = 0.001, az = 0.50*10.0, aw = 10.0, a_interval = 0.2500; //N=100
//double au = 2.0, av = 1.60, az = 10.0, a_interval = 0.2500; //original/
int DisplayInterval = 100; //1000
int DataSaveInterval = DisplayInterval;
int tsJump = 00000;

int flag_center = 1;
int flag_sigm = 0;
int flag_monitor = 1, flag = 0, flag_variable = 0;
int flag_drug = 0;
int save_flag = 1; //0:画像保存しない，1:画像保存する
int Flaaaag = 0, FlagDataSave = 1;

int count1 = 0, count2 = 0, count3 = 0;
int ParameterSetCount = 0, EndSetCount = 35;
//int save_csv = 0;q
int SignalingFlag[N] = {};
int Signaling = 0;

void capture(int *);

double t;
int winid;
int ts = 0; //現在のタイムステップ
int initial = 0;
FILE *fpValues, *fpV, *fpU, *fp, *fpZ;

char filename[256] = "", strFolderParaset[150] = "", strParaSetTxt[256] = "";
char strFolderMovie[256] = "", strMovieFiles[256] = "";
char filenameV[256] = "", filenameU[256] = "", filenameZ[256] = "";
void mouse(int, int, int, int);
void resize(int, int);
void FolderAndFileInit();
void init();
void FieldCycle();
void keyboard(unsigned char, int, int);
void idle();
void display();
void CsvRead();

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
  }
  //sprintf(filenamev, "Values.csv");
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(strParaSetTxt, "%s/ParaSet.txt", strFolderParaset);
  fpValues = fopen(strParaSetTxt, "w");
  fprintf(fpValues, "Du=%lf,Dv=%lf,Dz=%lf,u0=%lf,v0=%lf,z0=%lf,beta=%lf,alpha=%lf,NinitTadius_v=%d,u_in=%lf,k_1=%lf,k_2=%lf,k_vz=%lf,kh=%lf,kz=%lf,ku=%lf,kv=%lf,kf=%lf,rij=%lf,dx=%lf\n", Du, Dv, Dz, u0, v0, z0, beta, alpha, N_initTadius_v, u_in, k_1, k_2, k_vz, kh, kz, ku, kv, kf,rijInit,dx);
  fclose(fpValues);

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
}

void CsvRead()
{
  char tmp[256];
  char *p, *p2, *r;
  char **endptr = NULL;
  int fileNum = 0, i = 0;
  char strPos[128];

  //現在のディレクトリを取得
  p = getcwd(tmp, sizeof(tmp) - 1);
  printf("tmp=%s\n", tmp);
  strPos[0] = '/';
  strcpy(strPos + 1, tmp);

  //フォルダの番号を探す
  char str[] = "";
  int sizeNum = 0, sizeCharas = 0;
  r = strtok(strPos, "/"); //一文字目を切り出します。
  sizeCharas = strlen(p);
  strcat(str, "/");
  strcat(str, r);
  while (r = strtok(NULL, "/"))
  { //次の文字を切り出します。無ければNULLを返します。
    fileNum = (int)(strtol(r, endptr, 10));
    sizeNum = strlen(r);
    printf("%s\n", r);
    strcat(str, "/");
    strcat(str, r);
  }
  printf("fileNum=%d\n", fileNum);
  printf("%s\n", str);
  printf("sizeCharas=%d\n", sizeCharas);
  printf("sizeNum=%d\n", sizeNum);

  //parameters.csvのディレクトリを抜き出す
  char strDir[sizeCharas - sizeNum + 1];
  strncpy(strDir, tmp, sizeCharas - sizeNum);
  sizeCharas - sizeNum + 1;
  strDir[sizeCharas - sizeNum] = '\0';

  //parameterset.csvをオープンする
  strcat(strDir, "parameters.csv"); 
  printf("strDir=%s\n", strDir);
  char *fname = strDir; 
  FILE *fp;
  fp = fopen(fname, "r");
  if (fp == NULL)
  {
    printf("%sファイルが開けません\n", fname);
    exit(0);
  }

  //csvファイルからパラメータ値を代入する
  char temp;
  char StrC1[200], StrC2[200]; //, StrC3[200], StrC4[200], StrC5[20];
  printf("Check1\n");
  fscanf(fp, "%[^,],%s", StrC1, StrC2);
  printf("Check2\n");
  for (i = 0; i < fileNum; i++)
  { //csvファイルの1行目の数値はなぜか読み取ってくれない
    fscanf(fp, "%d,%lf", &ts, &Du);
  }
}

void init()
{

  int x, y, cell_count = 0;
  //ディスプレイ調整
  if (N / 100 == 1 && Flaaaag == 0)
  {
    dx = 0.1;
    xcenter = 10.0;
    zoom = 5;
    au = (0.02), av = 0.001, az = 3.0, aw = 10.0, a_interval = 0.2500;
    dxdx = dx * dx;
    monitor.SetZoom(zoom);
    monitor.SetCenter(xcenter, ycenter);
  }
  else if (N / 100 == 10 && Flaaaag == 0)
  {
    zoom = 28;
    xcenter = 100.0;
    aw = 0.1, au = (0.1), av = 0.1, az = 0.50, a_interval = 0.2500; //N=1000
    //dx = 0.01;
    dxdx = dx * dx;
    monitor.SetZoom(zoom);
    ycenter = 50.0;
    monitor.SetCenter(xcenter, ycenter);
  }

  //CsvRead();

  if (FlagDataSave == 1 && Flaaaag)
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
    if (count2 % 4 == 0 && count2 != 0)
    {
      count2 = 0;
      count3++;
    }
  }
  Flaaaag = 1;

  printf("\n\nParameterSetCount=%d\n", ParameterSetCount);
  printf("count1=%d count2=%d count3=%d\n", count1, count2, count3);

  double rijArray[]={1.0,10.0,100.0};
  double u0Array[] = {0.5,5.0,50.0,500.0};
  double DuArray[]={0.2,2.0,20.0,200.0};
  /*
  rijInit=rijArray[count1];
  u0=u0Array[count2];
  u_in=u0;
  Du=DuArray[count3];
  */
  
  printf("\n");
  ts = 0;
  t = 0.0;
  //DataCount = 0;

  //乱数初期化
  srand((unsigned int)time(NULL));
  //initialize Field

  for (x = 0; x < N; x++)
  {
    z_next[x] = 0.00000000;
    z_old[x] = z_next[x];
    u_next[x] = u0 + (rand() * ( 1.0) / (1.0 + RAND_MAX));
    u_old[x] = u_next[x];
    w_next[x] = 0.0;
    w_old[x] = w_next[x];

    SignalingFlag[x] = 0;
    //Coupling
    if ((x!=N/2))
    { //初期バイオフィルムの範囲
      //  u(x,y)=(ku0+1.500)+2.0*((ku0+1.500)*0.01)*(double)rand()/RAND_MAX-((ku0+1.500)*0.01);
      v_next[x] = v0 + (rand() * ( 1.0) / (1.0 + RAND_MAX));
      v_old[x] = v_next[x];
    }
    else
    {
      v_next[x] = 0.0;
      v_old[x] = v_next[x];
    }
  }

  for (x = 0; x < N; x++)
  {
    for (y = 0; y < N; y++)
    {
      if (x != y)
      {
        rij(x, y) = rijInit;
        dZij(x, y) = z_old[y] - z_old[x];
        dUij(x, y) = u_old[y] - u_old[x];
        dVij(x, y) = v_old[y] - v_old[x];
      }
      else
      {
        rij(x, y) = 0.0;
        dZij(x, y) = 0.0;
        dUij(x, y) = 0.0;
        dVij(x, y) = 0.0;
      }
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
  int x, y = 0;
  double f, h, ax;
  //zmin = z_next[o];
  for (i = 0; i < DisplayInterval; i++)
  {
    for (x = 0; x < N; x++)
    {
      LapU[x] = 0.0;
      LapV[x] = 0.0;
      LapZ[x] = 0.0;
      for (y = 0; y < N; y++)
      {
        if (rij(x, y) != 0.0)
        {
          LapU[x] += dUij(x, y) / rij(x, y);
          LapV[x] += dVij(x, y) / rij(x, y);
          LapZ[x] += dZij(x, y) / rij(x, y);
        }
      }
    }

    //方程式を解く
    for (x = 0; x < N; x++)
    {
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

      if (x == N/2) //int同士は切り捨て
      {
        u_next[x] += u_in * dt;
      }
    }

    for (x = 0; x < N; x++)
    {
      u_old[x] = u_next[x];
      v_old[x] = v_next[x];
      z_old[x] = z_next[x];
    }

    for (x = 0; x < N; x++)
    {
      for (y = 0; y < N; y++)
      {
        if (x != y)
        {
          dZij(x, y) = z_old[y] - z_old[x];
          dUij(x, y) = u_old[y] - u_old[x];
          dVij(x, y) = v_old[y] - v_old[x];
        }
        else
        {
          dZij(x, y) = 0.0;
          dUij(x, y) = 0.0;
          dVij(x, y) = 0.0;
        }
      }
    }

    ts++;
    t = ts * dt;
  }
  //printf("u=%lf   v=%lf   z=%lf\n",u_next[N/2],v_next[N/2],z_next[N/2]);
}

void keyboard(unsigned char key, int x, int y)
{
  int xx = 0;
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
    for (xx = 0; xx < N / 2; xx++)
    {
      z_next[xx] += 5.0;
      z_old[xx] = z_next[xx];
    }
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
  double base = 10.0 * (double)(j - N / 2), aaa = 3.0, bbb = 0.20;
  //yup=(interval*(N-1)+interval/2-N)/2.0;
  if (ts == 0 || flag_monitor == 1 && ts >= tsJump)
  {
    monitor.DrawLine(aaa * (double)(0 - N / 2), 0.0, aaa * (double)(2 - N / 2), 0.0, 0.1);

    for (j = 0; j < N; j++)
    {
      //ノードを描写
      monitor.SetAllColor(0.0, 0.0, 0.0);
      monitor.DrawCircle(aaa * (double)(j - N / 2), 0.0, 1.0);
      monitor.SetAllColor(1.0, 1.0, 1.0);
      monitor.DrawCircle(aaa * (double)(j - N / 2), 0.0, 0.90);
      //変数値を描写
      //v
      monitor.SetAllColor(0.0, 0.0, 1.0);
      monitor.DrawLine(aaa * (double)(j - N / 2) - bbb, 0.0, aaa * (double)(j - N / 2) - bbb, av * v_next[j], 5.0);
      //u
      monitor.SetAllColor(0.0, 1.0, 0.0);
      monitor.DrawLine(aaa * (double)(j - N / 2) , 0.0, aaa * (double)(j - N / 2) , au * u_next[j], 5.0);
      //z
      monitor.SetAllColor(1.0, 0.0, 0.0);
      monitor.DrawLine( aaa * (double)(j - N / 2)+bbb , 0.0 , aaa * (double)(j - N / 2)+bbb , az*z_next[j], 5.0 );
    }

    //printf("v0=%lf  u0=%lf\n", v_old[0],u_old[0]);
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
    y = y - 0.2;
    //sprintf(str, "Distance=%d", Distance);
    //monitor.String(x - 0.5, y, str);
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
    glFlush();
    glutSwapBuffers();

    if (ts > MAX_STEP) //|| Distance <= 2
    {
      if (ParameterSetCount >= EndSetCount)
      {
        fclose(fpValues);
        exit(0);
      }
      ParameterSetCount++;
      init(); //初期化＆再スタート
    }
  }

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
  monitor.SetZoom(zoom);
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
