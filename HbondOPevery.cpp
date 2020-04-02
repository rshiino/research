#include<bits/stdc++.h>
#include<fstream>
using namespace std;
#define PI 3.1415926535897932385
#define MAX_badget 250
#define MAX_atoms 100000
class atom{
    public:
        double coord[3];
        int imagine[3];
        int type;
};
string filename,inname,dumpname,dataname,outname;
string line;
int dumpnum,runnum;
int nowstep,atomnum,bondnum,Hbondnum,typenum;
ifstream in,data,dump;
ofstream out;
vector<vector<vector<vector<int>>>> dbadget(MAX_badget,vector<vector<vector<int>>>(MAX_badget,vector<vector<int>>(MAX_badget)));
vector<vector<vector<vector<int>>>> abadget(MAX_badget,vector<vector<vector<int>>>(MAX_badget,vector<vector<int>>(MAX_badget)));
vector<vector<int>> topology(MAX_atoms);//隣接リスト
vector<string> atomline;//原子の情報
vector<pair<int,int>> Hbondpair;//水素結合をしてる組
double maxlength=4.0;//水素結合の基準長さ
double maxcostheta=0.0;//水素結合の基準角度
double cell[3][2];//xlo,xhi,ylo,yhi,zlo,zhi
int badgetnum[3];//バゲットの区切り数　毎回変える
vector<atom> atomlist(MAX_atoms);
set<int> donor,acceptor;

bool checkFileExistence(const string& str){//ファイルの存在確認
    ifstream ifs(str);
    return ifs.is_open();
}

bool checkline(string line,string keyword){//先頭がkeywordかどうかの確認
    if(line.size()<keyword.size()) return false;
    for(int i=0;i<keyword.size();i++) if(line[i]!=keyword[i]) return false;
    return true;
}

double length(int a,int b){//a-bの距離
    double coord[2][3];
    double len=0;
    for(int i=0;i<3;i++){
        coord[0][i]=atomlist[a].coord[i]+(cell[i][1]-cell[i][0])*atomlist[a].imagine[i];
        coord[1][i]=atomlist[b].coord[i]+(cell[i][1]-cell[i][0])*atomlist[b].imagine[i];
        len+=(coord[0][i]-coord[1][i])*(coord[0][i]-coord[1][i]);
    }
    len=sqrt(len);
    return len;
}

double ISlength(int a,int b){
    return length(a,b)<maxlength;
}

bool ISangle(int a,int b,int c){//a-b-cの角度
    double coord[3][3];
    double vec[2][3];
    double ip=0;
    for(int i=0;i<3;i++){
        coord[0][i]=atomlist[a].coord[i]+(cell[i][1]-cell[i][0])*atomlist[a].imagine[i];
        coord[1][i]=atomlist[b].coord[i]+(cell[i][1]-cell[i][0])*atomlist[b].imagine[i];
        coord[2][i]=atomlist[c].coord[i]+(cell[i][1]-cell[i][0])*atomlist[c].imagine[i];
        vec[0][i]=coord[1][i]-coord[0][i];
        vec[1][i]=coord[1][i]-coord[2][i];
        ip+=vec[0][i]*vec[1][i];
    }
    double cos=ip/length(a,b)/length(b,c);
    return cos<maxcostheta;
}

double func(vector<double> coeff,double x){
    return coeff[0]*x*x*x+coeff[1]*x*x+coeff[2]*x+coeff[3];
}
double bisec(vector<double> coeff){
    double left=-1,right=1;
    bool lplus=(func(coeff,left)>0);
    bool rplus=(func(coeff,right)>0);
    while(right-left>1e-6){
        double center=(right+left)/2;
        bool cplus=(func(coeff,center)>0);
        if(lplus==cplus) left=center;
        else right=center;
    }
    return (left+right)/2;
}

vector<double> calceigen(vector<vector<double>> a){
    vector<double> coeff(4,0);
    coeff[0]=1;
    coeff[1]-=(a[0][0]+a[1][1]+a[2][2]);
    coeff[2]+=(a[0][0]*a[1][1]+a[1][1]*a[2][2]+a[2][2]*a[0][0]);
    coeff[2]+=(a[0][1]*a[1][0]+a[1][2]*a[2][1]+a[2][0]*a[0][2]);
    coeff[3]-=(a[0][0]*a[1][1]*a[2][2]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]);
    coeff[3]-=(a[0][2]*a[1][1]*a[2][0]+a[0][0]*a[1][2]*a[2][1]+a[0][1]*a[1][0]*a[2][2]);
    double lambda=bisec(coeff);
}

void row_reduction(vector<vector<double>> &inv_a,vector<vector<double>> a){//掃き出し法
    double buf; //一時的なデータを蓄える
    int n=3;
    //単位行列を作る
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
        inv_a[i][j]=((i==j)?1.0:0.0);
        }
    }
    for(int i=0;i<n;i++){
        buf=1/a[i][i];
        for(int j=0;j<n;j++){
            a[i][j]*=buf;
            inv_a[i][j]*=buf;
        }
        for(int j=0;j<n;j++){
            if(i!=j){
                buf=a[j][i];
                for(int k=0;k<n;k++){
                    a[j][k]-=a[i][k]*buf;
                    inv_a[j][k]-=inv_a[i][k]*buf;
                }
            }
        }   
    }
}

void input(){//ファイル名を入力する
    while(true){
        try{
            cout<<"////enter filename////\n";//filename(拡張子抜き)を聞く
            cin.clear();
            cin.ignore(1024,'\n');
            cin.seekg(0,ios::beg);
            cin>>filename;
            inname=filename+"lmp_set";//winmostarで書くと.inではなく.lmp_setで出る
            dumpname=filename+".dump";
            dataname=filename+".data";
            outname=filename+"hbond.out";
            if(!checkFileExistence(inname)){
                cout<<"////inputfile is not exist////\n";
                throw exception();
            }
            else if(!checkFileExistence(dumpname)){
                cout<<"////dumpfile is not exist////\n";
                throw exception();
            }
            else if(!checkFileExistence(dataname)){
                cout<<"////datafile is not exist////\n";
                throw exception();
            }
        }
        catch(...){
            cout<<"////illegal input////\n";
            continue;
        }
        break;
    }
}

void readin(){//inファイルを読み込む
    while(true){
        getline(in,line);
        if(checkline(line,"dump")) {
            cout<<line<<endl;
            continue;
        }
        else if(checkline(line,"thermo")) {
            cout<<line<<endl;
            continue;
        }
        else if(checkline(line,"run")) {
            cout<<line<<endl;
            break;
        }
    }
    cout<<endl;
    cout<<"////dump interval////"<<endl;
    cin>>dumpnum;
    cout<<"////running steps////"<<endl;
    cin>>runnum;
}

void readdata(){//dataファイルを読み込む
    data>>atomnum>>line>>bondnum;
    while(true){
        getline(data,line);
        if(line=="") break;
    }
    data>>typenum;
    while(true){
        getline(data,line);
        if(line=="Masses") break;
    }
    getline(data,line);
    for(int i=0;i<typenum;i++){
        getline(data,line);
        atomline.push_back(line);
    }
    while(true){
        getline(data,line);
        if(line=="Bonds") break;
    }
    for(int i=0;i<bondnum;i++){
        int a,b,c,d;
        cin>>a>>b>>c>>d;
        topology[c].push_back(d);
        topology[d].push_back(c);
    }
}

void defCluster(){//対象の原子組の条件を受け取る
    int d;string s;
    donor={};
    acceptor={};
    for(auto s:atomline) cout<<s<<endl;
    cout<<"////which type donor is?////"<<endl;
    while(true){
        try{
            cin>>d;
            if(d>=typenum||d<1){
                cout<<"illegal number"<<endl;
                throw exception();
            } 
            donor.insert(d);
        }
        catch(...){
            cout<<"finished? put y or n"<<endl;
            cin>>s;
            if(s=="y") break;
            else continue;
        }
    }
    cout<<endl;
    cout<<"////which type acceptor is?////"<<endl;
    while(true){
        try{
            cin>>d;
            if(d>=typenum||d<1){
                cout<<"illegal number"<<endl;
                throw exception();
            } 
            acceptor.insert(d);
        }
        catch(...){
            cout<<"finished? put y or n"<<endl;
            cin>>s;
            if(s=="y") break;
            else continue;
        }
    }
    cout<<"////how long Hbond border is?////"<<endl; 
    cin>>maxlength;
    cout<<"////how small Hbond angle border is?////"<<endl;
    cin>>maxcostheta;
    maxcostheta=cos(maxcostheta*PI/180.0);
}

void putToBadget(int number,int DorA){//原子をbadgetに入れる
    double atomcoord[3];
    int badgetcoord[3];
    if(DorA==0){//Donor
        for(int i=0;i<3;i++){
            atomcoord[i]=atomlist[number].coord[i]-cell[i][0];
            badgetcoord[i]=min(badgetnum[i]-1,int(ceil(atomcoord[i]/maxlength)));
        }   
        dbadget[badgetcoord[0]][badgetcoord[1]][badgetcoord[2]].push_back(number);
    }
    else if(DorA==1){//Acceptor
        for(int i=0;i<3;i++){
            atomcoord[i]=atomlist[number].coord[i]-cell[i][0];
            badgetcoord[i]=min(badgetnum[i]-1,int(ceil(atomcoord[i]/maxlength)));
        }   
        abadget[badgetcoord[0]][badgetcoord[1]][badgetcoord[2]].push_back(number);        
    }
}

void readdump(){//dumpファイルを読み込む
    while(true){
        getline(dump,line);
        if(line=="ITEM: BOX BOUNDS pp pp pp") break;
    }
    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++) cin>>cell[i][j];
        badgetnum[i]=ceil((cell[i][1]-cell[i][0])/maxlength);
    }
    getline(dump,line);
    for(int i=0;i<atomnum;i++){
        int number;cin>>number;
        cin>>atomlist[number].type
        >>atomlist[number].coord[0]
        >>atomlist[number].coord[1]
        >>atomlist[number].coord[2]
        >>atomlist[number].imagine[0]
        >>atomlist[number].imagine[1]
        >>atomlist[number].imagine[2];
        if(acceptor.find(atomlist[number].type)!=acceptor.end()){//acceptorだったらbadgetに入れる
            putToBadget(number,0);
        }
        else if(donor.find(atomlist[number].type)!=donor.end()){//donorだったら(ry
            putToBadget(number,1);
        }
    }
}

bool checkHbond(int anumber,int dnumber){
    if(topology[dnumber][0]==anumber) return false;//繋がっているとアウト
    else if(ISlength(dnumber,anumber)) return false;//結合距離
    else if(ISangle(topology[dnumber][0],dnumber,anumber)) return false;//角度
    else return true;
}

void searchHbond(){//Hbondを列挙する
    for(int i=0;i<badgetnum[0];i++){
        for(int j=0;j<badgetnum[1];j++){
            for(int k=0;k<badgetnum[2];k++){
                for(auto anumber:abadget[i][j][k]){
                    for(int i2=i-1;i2<i+2;i2++){
                        for(int j2=j-1;j2<j+2;j2++){
                            for(int k2=k-1;k2<k+2;k2++){
                                for(auto dnumber:dbadget[(i2+badgetnum[0])%badgetnum[0]][(j2+badgetnum[1])%badgetnum[1]][(k2+badgetnum[2])%badgetnum[2]]){
                                    if(checkHbond(anumber,dnumber)){
                                        Hbondpair.push_back(make_pair(anumber,dnumber));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void calcOP(){//OPを計算する
    Hbondnum=Hbondpair.size();
    vector<vector<double>> a(3,vector<double> (3,0));
    vector<vector<double>> HbondVector(Hbondnum,vector<double> (3));
    for(int i=0;i<Hbondnum;i++){
        pair<int,int> p=Hbondpair[i];
        for(int j=0;j<3;j++){
            HbondVector[i][j]=atomlist[p.first].coord[j]+(cell[j][1]-cell[j][0])*atomlist[p.first].imagine[j]
                             -atomlist[p.second].coord[j]+(cell[j][1]-cell[j][0])*atomlist[p.second].imagine[j];
        }
        double len=length(p.first,p.second);
        for(int j=0;j<3;j++){
            HbondVector[i][j]/=len;
        }   
    }
    for(int i=0;i<Hbondnum;i++){
        for(int j1=0;j1<3;j1++){
            for(int j2=0;j2<3;j2++){
                a[j1][j2]+=HbondVector[i][j1]*HbondVector[i][j2];
            }
        }
    }
    for(int j1=0;j1<3;j1++){
        for(int j2=0;j2<3;j2++){
            a[j1][j2]*=3;
            a[j1][j2]/=2;
            a[j1][j2]/=Hbondnum;
            if(j1==j2) a[j1][j2]-=0.5;        
        }
    }
    vector<double> eigenvalue=calceigen(a);
    vector<vector<double>> inv_a(3,vector<double> (3));//逆行列
    row_reduction(inv_a,a);
}

void output(){//結果を出力する

}

int main(){
    input();
    in.open(inname);
    readin();//inファイルを読み込む　済
    in.close();
    data.open(dataname);
    readdata();//dataファイルを読み込む　済
    data.close();
    defCluster();//Hbondの定義をする　済
    dump.open(dumpname);
    out.open(outname);
    cout<<endl;
    cout<<"//////////////////////////"<<endl;
    cout<<"////calculation start/////"<<endl;
    cout<<"//////////////////////////"<<endl;
    cout<<endl;
    for(nowstep=0;nowstep<=dumpnum;nowstep+=runnum){
        for(int i=0;i<3;i++) badgetnum[i]=-1;
        Hbondnum={};
        abadget={};
        dbadget={};
        readdump();//dumpファイルを読み込む+badgetに放り込む　済
        searchHbond();//各donorについて隣接badgetまでHbond探索をする　済
        calcOP();//HbondのOPを計算する 途中
        output();//OPと数を出力する
        cout<<nowstep<<" steps finished"<<endl;//計算状況を出力する
    }
    cout<<endl;
    cout<<"//////////////////////////"<<endl;
    cout<<"////all steps finished////"<<endl;
    cout<<"//////////////////////////"<<endl;
}
