#include<bits/stdc++.h>
#include<fstream>
using namespace std;
typedef long long ll;
#define For(i,n,k) for(ll i=(n);i<(k);i++)
#define PI 3.1415926535897932385
#define N 3 
string gomi;//ゴミ捨て場
string filename;//ファイル名
string inname,dumpname,dataname,outname;//ファイル名＋拡張子
ll nowstep,atomnum,bondnum,Hbondnum;//ステップ数・原子数・結合数・水素結合数
ll dumpnum,runnum,calcnum;//dump出力間隔・総ステップ数・dump出力回数
ifstream fin,datain;
ofstream fout;
vector<int> donor={11,5,9},acceptor={4,10,1,7,13,14};//CPならこの組を入れる
vector<vector<int>> topology;//隣接リスト
vector<pair<int,int>> Hbondpair;//水素結合をしてる組
double maxlength=4;//水素結合の基準長さ
double maxcostheta=0;//水素結合の基準角度
double cell[3];//セルサイズ
double a[3][3],U[3][3];//対角化に使う
class atom{
    public:
        double coord[3];
        int imagine[3];
        int id;
        int type;
};
double length2(atom a,atom b){//ABの距離の二乗を返す　可換
    double len=0;
    For(c,0,3){
        len+=(a.coord[c]-b.coord[c])*(a.coord[c]-b.coord[c]);
    }
    return len;
}
double cosangle(atom a,atom b,atom c){//角ABCのコサインを返す　可換じゃないので注意
    vector<double> ba(3),bc(3);
    for(int i=0;i<3;i++){
        ba.at(i)=a.coord[i]-b.coord[i];
        bc.at(i)=c.coord[i]-b.coord[i];
    }
    double costheta=0;
    for(int i=0;i<3;i++){
        costheta+=ba.at(i)*bc.at(i);
    }
    double len=sqrt(length2(a,b)*length2(b,c));
    costheta/=len;
    return costheta;        
}
vector<atom> atomlist;
void diagonalization(double a[N][N],double U[N][N])
{
        int i,j,l,m,p,q,count;
        double max,theta;
        double u[N][N],oldU[N][N],newa[N][N];
/*計算結果としてだされた直行行列を
格納するための配列を単位行列に初期化しておく。*/
        for(p=0;p<N;p++) {
                for(q=0;q<N;q++) {
                        if(p==q) U[p][q]=1.0;
                        else U[p][q]=0.0;
                }
        }
 
        for(count=0;count<=10000;count++) {
/*配列olduは新たな対角化計算を行う前に
かけてきた直行行列を保持する。*/
                for(p=0;p<N;p++) {
                        for(q=0;q<N;q++) {
                                oldU[p][q]=U[p][q];
/*対角化するときの個々の直交行列を入れる配列uの初期化。(単位行列にする。)*/
                                if(p==q) u[p][q]=1.0;
                                else u[p][q]=0.0;
                        }
                }
/*非対角要素の中から絶対値の最大のものを見つける*/
                max=0.0;
                for(p=0;p<N;p++) {
                        for(q=0;q<N;q++) {
                                if(max<fabs(a[p][q]) && p!=q) {
                                        max=fabs(a[p][q]);
/*その最大のものの成分の行と列にあたる数を記憶しておく。*/
                                        i=p;
                                        j=q;
                                }
                        }
                }
/*先ほど選んだ最大のものが指定の値より小さければ対角化終了*/
                if(fabs(a[i][j]) < 1.0e-12) {
                        break;
                }
/*条件によってシータの値を決める*/
                if(a[i][i]==a[j][j]) theta=PI/4.0;
                else theta=atan(-2*a[i][j]/(a[i][i]-a[j][j]))/2.0;
/*ここでこのときに実対称行列にかける個々の直行行列uが決まるが
特にここでの計算の意味はない。(する必要はない。)*/
                u[i][j]=sin(theta);
                u[j][i]=-sin(theta);
                u[i][i]=u[j][j]=cos(theta);
/*ここでいままで実対称行列にかけてきた直行行列を配列Uに入れる。*/
                for(p=0;p<N;p++) {
                        U[p][i]=oldU[p][i]*cos(theta)-oldU[p][j]*sin(theta);
                        U[p][j]=oldU[p][i]*sin(theta)+oldU[p][j]*cos(theta);
                }
//対角化計算によってでた新たな実対称行列の成分を配列newaに入れる。
                newa[i][i]=a[i][i]*cos(theta)*cos(theta)
                        +a[j][j]*sin(theta)*sin(theta)
                        -2.0*a[i][j]*sin(theta)*cos(theta);
                newa[j][j]=a[i][i]*sin(theta)*sin(theta)
                        +a[j][j]*cos(theta)*cos(theta)
                        +2.0*a[i][j]*sin(theta)*cos(theta);
                newa[i][j]=newa[j][i]=0.0;
                for(l=0;l<N;l++) {
                        if(l!=i && l!=j) {
                                newa[i][l]=a[i][l]*cos(theta)
                                    -a[j][l]*sin(theta);
                                newa[l][i]=newa[i][l];
                                newa[j][l]=a[i][l]*sin(theta)
                                    +a[j][l]*cos(theta);
                                newa[l][j]=newa[j][l];
                        }
                }
                for(l=0;l<N;l++) {
                        for(m=0;m<N;m++) {
                            if(l!=i && l!=j && m!=i && m!=j) {
                                    newa[l][m]=a[l][m];
                            }
                        }
                }
//次の対角化計算を行う行列の成分を配列aへ上書きする。
                for(p=0;p<N;p++) {
                        for(q=0;q<N;q++) {
                                a[p][q]=newa[p][q];
                        }
                }
 
                if(count==10000) {
                        printf("対角化するために作業を繰り返す必要がまだあり.\n");
                }
        }
}
vector<vector<int>> readData(){//dataファイルから原子タイプと結合を読み込んで、ドナーアクセプターの原子タイプを選ばせる。戻り値は隣接リスト表現のトポロジー
    donor={};acceptor={};
    gomi="";
    while(gomi!="Masses"){
        datain>>gomi;
        if(gomi=="Amber_Cornell_ext.params+b)"){
            datain>>gomi;
            atomnum=stoi(gomi);
            datain>>gomi;
            datain>>gomi;
            bondnum=stoi(gomi);
            break;
        }
    }
    while(gomi!="Masses"){
        datain>>gomi;
    }
    cout<<endl;
    string s;datain>>s;
    while(s!="Atoms"){
        For(i,0,4){
            cout<<s<<" ";
            datain>>s;
        }
        cout<<endl;
    }
    cout<<endl;
    int num;
    cout<<"////How many donors?////"<<endl;
    cin>>num;
    cout<<"////numbers?////"<<endl;
    For(i,0,num){
        int d;cin>>d;
        donor.push_back(d);
    }
    cout<<endl;
    cout<<"////How many acceptor?////"<<endl;
    cin>>num;
    cout<<"////numbers?////"<<endl;
    For(i,0,num){
        int d;cin>>d;
        acceptor.push_back(d);
    }
    while(gomi!="Bonds") datain>>gomi;
    vector<vector<int>> top(atomnum+1);
    For(i,0,bondnum){
        int p,q;
        datain>>gomi>>gomi;
        datain>>p>>q;
        datain>>gomi>>gomi>>gomi>>gomi;  
        top.at(p).push_back(q);
        top.at(q).push_back(p);      
    }
    return top;   
}
vector<atom> readDump(){//1ステップ分の原子の情報を読み込む
    For(i,0,2)fin>>gomi;//最初の方はゴミが多い
    if(gomi==""){
        vector<atom> x(1);
        x.at(0).id=-1;
        return x;
    }
    fin>>nowstep;
    fout<<nowstep<<" ";
    For(i,0,4)fin>>gomi;
    fin>>atomnum;
    For(i,0,6)fin>>gomi;
    For(i,0,3){
        double lo,hi;
        fin>>lo>>hi;
        cell[i]=hi-lo;
    }
    For(i,0,10)fin>>gomi;
    vector<atom> atomlist(atomnum+1);
    For(i,1,atomnum+1){
        fin>>atomlist.at(i).id
        >>atomlist.at(i).type
        >>atomlist.at(i).coord[0]
        >>atomlist.at(i).coord[1]
        >>atomlist.at(i).coord[2]
        >>atomlist.at(i).imagine[0]
        >>atomlist.at(i).imagine[1]
        >>atomlist.at(i).imagine[2];
        For(j,0,3){
            atomlist.at(i).coord[j]+=cell[j]*atomlist.at(i).imagine[j];
        }
    }
    return atomlist;
}
vector<pair<int,int>> searchHbond(){//Hbondのペアを書き出す　今は原子タイプ・距離・角度・bondなし
    vector<pair<int,int>> Hbondpair;
    For(i,0,atomnum){
        for(auto Donor:donor){
            if(atomlist.at(i).type==Donor){//iはドナー
                For(j,0,atomnum){
                    if(topology.at(i).at(0)==j) continue;//結合で繋がっていたらだめ(ドナーは水素==結合先は一つのハズなのでこう書く)
                    for(auto Acceptor:acceptor){
                        if(atomlist.at(j).type==Acceptor){//jはアクセプター
                            if(length2(atomlist.at((topology.at(i).at(0))),atomlist.at(j))<maxlength*maxlength){//距離の平方で判断 
                                if(cosangle(atomlist.at((topology.at(i).at(0))),atomlist.at(i),atomlist.at(j))<maxcostheta){//角度がなんとか度以上なら
                                    Hbondpair.push_back(make_pair(i,j));
                                }
                            }
                            break;
                        }
                    }             
                }
                break;
            }
        }
    }
    Hbondnum=Hbondpair.size();
    return Hbondpair;
}
void calcOP(vector<pair<int,int>> HbondPair){
    Hbondnum=HbondPair.size();
    vector<vector<double>> HbondVector(Hbondnum,vector<double> (3));
    For(i,0,Hbondnum){
        For(j,0,3){
            HbondVector.at(i).at(j)=atomlist.at(HbondPair.at(i).first).coord[j]-atomlist.at(HbondPair.at(i).second).coord[j];
        }
        double length=sqrt(HbondVector.at(i).at(0)*HbondVector.at(i).at(0)+HbondVector.at(i).at(1)*HbondVector.at(i).at(1)+HbondVector.at(i).at(2)*HbondVector.at(i).at(2));
        For(j,0,3){
            HbondVector.at(i).at(j)/=length;
        }
    }
    double a[3][3]={0};
    For(i,0,Hbondnum){
        For(j1,0,3){
            For(j2,0,3){
                a[j1][j2]+=HbondVector.at(i).at(j1)*HbondVector.at(i).at(j2);
            }
        }
    }
    For(j1,0,3){
        For(j2,0,3){
            a[j1][j2]*=3;
            a[j1][j2]/=2;
            a[j1][j2]/=Hbondnum;
            if(j1==j2) a[j1][j2]-=0.5;
        }
    }
    double U[N][N];//固有行列
    diagonalization(a,U);
    double x=-1;
    int t=0;
    For(j1,0,3) {
        if(x<3*a[j1][j1]){
            t=j1;
            x=3*a[j1][j1];
        }
    }
    fout<<" OP= "<<x<<" director= "<<U[t][0]<<" "<<U[t][1]<<" "<<U[t][2]<<endl;
}
int main(){
    cout<<"////filename is////"<<endl;
    cin>>filename;
    inname=filename+".in",dumpname=filename+".dump",dataname=filename+".data",outname=filename+"Hbond.log";
    datain.open(inname);
    while(true){
        datain>>gomi;
        if(gomi=="dump"){
            datain>>gomi;
            if(gomi=="1"){
                datain>>gomi>>gomi>>dumpnum;
            }
        }
        else if(gomi=="run"){
            datain>>runnum;
            break;
        }
    }
    calcnum=runnum/dumpnum;
    datain.close();
    cout<<calcnum<<endl;
    datain.open(dataname);
    topology=readData();
    datain.close();
    fout.open(outname);
    fin.open(dumpname);
    atomlist=readDump();
    vector<pair<int,int>> HbondPair=searchHbond();
    fout<<"timestep xx  yy  zz  Hbondnum = "<<Hbondnum<<endl;
    calcOP(HbondPair);
    For(i,0,calcnum){
        atomlist=readDump();
        calcOP(HbondPair);    
    }
    cout<<"////finish////"<<endl;
    return 0;
}
