#include <iostream>
#include <math.h>

using namespace std;

const int SIZE = 7; //określenie rozmiaru macierzy
void sub2(long double t1[], long double t2[], long double *t3);
void wypisz(long double a[SIZE], long double b[SIZE-1], long double c[SIZE-1]);

int main(int argc, const char * argv[]) {
    
    long double diag[SIZE]; //wektor z 4, diagonalny
    long double upper[SIZE-1]; //wektor z 1, powyżejdiagonali
    long double lower[SIZE-1]; //wektor z 1, poniżej
    long double x[SIZE]; //wektor z wynikami końcowymi x1 ... x7
    long double z[SIZE]; //wektor z wynikami z1...z7
    long double q[SIZE]; //wektor z wynikami q1...q7
    long double c[SIZE]; //wektor c - zawiera liczby od 1 do 7
    long double u[SIZE]; //wektor u = [1,0,0,0,0,0,1]^T
    
    long double C,S,pom; //zmienne pomocnicze
    
    for(int i=0;i<SIZE;i++){ //wypełnianie tablic A, C danymi
        
        diag[i] = 4;
        if(i<SIZE-1){
            upper[i] = 1;
            lower[i]=1;}
        
        c[i] = i+1;
        u[i] = 0;
        if(i==SIZE-1||i==0) u[i]=1;
        
    }
    
    //Od macierzy A1 odejmuję u*v^T ==>> uzyskam macierz A (znikną jedynki w rogach)
    diag[0] = 3;
    diag[SIZE-1]=3;
    
    cout<<"Macierz A: "<<endl;
    wypisz(diag,upper,lower);
    

    /*
     
     Do rozwiązania problemu użyję algorytmu Shermana-Morrisona
     
     Aby rozwiązać układ równań wykorzystam obroty givensa oraz metodę back substitution.
     Program skończy pracę w czasie liniowym zamiast O(N^3)
     
     Na początek od macierzy A1 odejmę macierz u*v^T ==>> uzyskam macierz A
     (u mnie złożona z 3 wektorów reprezentujących diagonale: lower, diag, upper)
     
     Następnie obliczam wektor Z = A^(-1)*b
     Oraz wektor Q = A^(-1)*U
     
     Obydwa wyliczam poprzez rozwiązanie równiania Ax = b
     
     Korzystam z Obrotów Givensa by wykorzystać strukturę macierzy.
     
     Ten problem pozwala użyć zmodyfikowanego algorytmu back substition
     który wykonuje się w czasie O(N) zamist O(N^2)
     
     Ostatecznie wyliczam rozwiązania x1 .. xn korzystając z
     wektorów Q oraz Z
     
     Wektor x =  Z - ( (v^T)*Z / 1+(v^T)*q )*q
     
     */

    
    for(int i=0; i<SIZE-1;i++){ //wykonujemy n-1 obrotów (w tym przypadku 6)
        
        C = diag[i]/sqrt(diag[i]*diag[i]+lower[i]*lower[i]); //wyliczamy c,s potrzebne do wykonania obrotu
        S = lower[i]/sqrt(diag[i]*diag[i]+lower[i]*lower[i]);
        
        cout<<"Wartości C, S"<<endl;
        cout<<"C: "<<C<<" S: "<<S<<endl<<endl;
        
        //mnożenie części macierzy C i U przez G - macierz Givensa
        pom = c[i];
        c[i] = C * c[i] + S*c[i+1];
        c[i+1]= (-S) * pom + C*c[i+1];
        
        pom = u[i];
        u[i] = C * u[i] + S*u[i+1];
        u[i+1]= (-S) * pom + C*u[i+1];
        
        //Mnożymy części macierzy A (u nas macierz a to 3 wektory) przez macierz Givensa
        if(i==SIZE-2)
        {
            pom = diag[i];
            diag[i] = C*diag[i] + S*lower[i];
            pom = upper[i];
            upper[i] = C*upper[i] + S*diag[i+1];
            diag[i+1] = (-S)*pom + C*diag[i+1];
        }
        else
        {
            pom = diag[i];
            diag[i]= C*diag[i] + S*lower[i];
            lower[i]= (-S)*pom + C*lower[i];
            pom=upper[i];
            upper[i] = C * upper[i] + S*diag[i+1];
            diag[i+1] = (-S)*pom + C*diag[i+1];
            pom=lower[i];
            lower[i] =  S*upper[i+1];
            upper[i+1] = C*upper[i+1];
            
        }
        //dzięki temu że mnożymy tylko część elementów (zera są pomijane), złożoność jest wyrażnie mniejsza
        
        //wypisuję macierze po każdym obrocie
        cout<<endl<<"Obrót: "<<i+1<<endl;
        cout<<"Macierz A: "<<endl;
        wypisz(diag,upper,lower);
        cout<<endl<<"Macierz C: "<<endl;
        for(int i=0;i<SIZE;i++) cout<<c[i]<<endl;
        cout<<endl<<"Wektor U: "<<endl;
        for(int i=0;i<SIZE;i++) cout<<u[i]<<endl;
    }
    
    cout<<endl<<"Po obrotach Givensa: "<<endl<<endl;
    
    //back substitution
    for(int i=SIZE-1;i>=0;i--){
        if(i==SIZE-1) {
            z[i] = c[i]/diag[i]; //dla najniższego wiersza
            q[i] = u[i]/diag[i];
        }
        else if(i==SIZE-2) {
            z[i]=(c[i]-upper[i]*z[i+1])/diag[i]; //dla przedostatniego wiersza
            q[i]=(u[i]-upper[i]*q[i+1])/diag[i];
        }
        else {
            z[i]=(c[i]-upper[i]*z[i+1]-lower[i]*z[i+2])/diag[i];
            q[i]=(u[i]-upper[i]*q[i+1]-lower[i]*q[i+2])/diag[i];
        }
        /* dla pozostałych
         warto zwróćić uwagę, że maksymalnie 3 elementy w wierszu są niezerowe dlatego
         rozwiązanie jest w czasie liniowym O(N) zamiast O(N^2)
         */
    }
    //wypisanie wektora Z i Q
    cout<<endl<<"Wektor Z"<<endl;
    for(int i=0;i<SIZE;i++) cout<<"z["<<i+1<<"] = "<<z[i]<<endl;
    cout<<endl<<"Wektor Q"<<endl;
    for(int i=0;i<SIZE;i++) cout<<"q["<<i+1<<"] = "<<q[i]<<endl;
    
    //Wektor x =  Z - ( (v^T)*Z / 1+(v^T)*q )*q
    //Warto zauważyć, że (v^T)*Z / 1+(v^T)*q jest liczbą.
    
    pom = (z[0]+z[SIZE-1])/(1+q[0]+q[SIZE-1]);
    cout<<endl<<"(v^T)*Z / 1+(v^T)*q = "<<pom<<endl;
    for(int i=0;i<SIZE;i++) q[i]*=pom; //mnożymy wektor Q przez wyliczone pom
    sub2(z,q,x); //X = Z - Q
    cout<<endl<<"Wektor X"<<endl;
    for(int i=0;i<SIZE;i++) cout<<"x["<<i+1<<"] = "<<x[i]<<endl;
    
    return 0;
    
}

void sub2(long double t1[], long double t2[], long double *t3){ //odejmowanie wektorów
    for(int i=0; i<SIZE;i++){ //odejmujemy każdy element od każdego
        t3[i]=(t1[i] - t2[i]);
    }
}

//odejmuje macierz U*V^T
void sub(long double tab[SIZE][SIZE]){
    tab[SIZE-1][0] -=1;
    tab[SIZE-1][SIZE-1] -=1;
    tab[0][0] -=1;
    tab[0][SIZE-1] -=1;
}

void wypisz(long double a[SIZE], long double b[SIZE-1], long double c[SIZE-1]){
    for(int i=0;i<SIZE;i++){
        for(int j=0;j<SIZE;j++)
        {
            if(i==j) cout<<a[i]<<"\t\t";
            else if(i==j-1) cout<<b[i]<<"\t\t";
            else if(i==j-2) cout<<c[i]<<"\t\t";
            else cout<<0<<"\t\t";
        }
        cout<<endl;
    }
}

/*
Results: 

Vector Z
 z[1] = 0.2281
 z[2] = 0.315699
 z[3] = 0.509103
 z[4] = 0.647887
 z[5] = 0.899347
 z[6] = 0.754723
 z[7] = 2.08176
 
Vector Q
 q[1] = 0.366197
 q[2] = -0.0985915
 q[3] = 0.028169
 q[4] = -0.0140845
 q[5] = 0.028169
 q[6] = -0.0985915
 q[7] = 0.366197
 
 (v^T)*Z / 1+(v^T)*q = 1.33333
 
Vector X
 x[1] = -0.260163
 x[2] = 0.447154
 x[3] = 0.471545
 x[4] = 0.666667
 x[5] = 0.861789
 x[6] = 0.886179
 x[7] = 1.5935

*/