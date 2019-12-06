#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdint>
#include <cstring>
#include "bigint.h"
#include "sha256.h"

struct sig
{
    BigInteger gamma,delta;
};

struct data
{
    BigInteger p,q,alpha,a,beta;
};

struct data generateData(const BigInteger&p,const BigInteger&q,const BigInteger&g);
struct sig sign(const char* s ,const BigInteger&k,const data&d);
BigInteger hash2BigInteger(const uint8_t* hash);
bool verify(const char*s,struct sig& si,const data&d);
unsigned swap_uint32(const unsigned &x);

int main()
{
    const BigInteger p(
                "1707678874597591438281790056162780764022"
                "6368296588893844575071209893079489801905"
                "3350746183899942190230348202497274526845"
                "0971426943787028521830242884057779297127"
                "0246842811869895607758402548961537325231"
                "5365141656478879448986019442115302380086"
                "3554978522744824210361162491424596543829"
                "6675597253892346698962539341558898420773"
                "6649575603690287581963025244894396712853"
                "7079518200995492394755758439652774419449"
                "8301169922818875796512663823496778703396"
                "1590238350847934997270962962664692175276"
                "9644322862898792365394664775663353564080"
                "4959907295932525165462231825906598427004"
                "4513630454576620978846504513416166782476"
                "11110805615245953");
    const BigInteger q(
                "8586727146255868338239290877850142140482"
                "1495432764570749371157194193678385377");
    const BigInteger g(
                "1489479783570884039322743825678085797303"
                "9787202586450772583528386292462324776546"
                "9915618585578050942445721284027237277971"
                "5103464991929255083417569033947596080262"
                "5113594958834589938987183744188611464432"
                "9324517954036780169639474138503476539874"
                "0062868429040839175767996081826767720952"
                "7400153134357104789680912688993778967215"
                "9489244263852388187278986533578513121404"
                "6749752670937002093723163337784552790647"
                "1224864732671884047236017222906428109837"
                "0118363345608737958610883275984475440284"
                "8506372150061249575700714051186016460386"
                "7182769088925648376250536084051656930390"
                "8514556665764958000126588537154781939041"
                "71187728662780119");
    const char s[]="SchoolofDataandComputerScience,Sunyat-senUniversity";

    clock_t t=clock();
    struct data d=generateData(p,q,g);
    BigInteger k=BigInteger::ZERO;
    bool ans=0;
    while(1)
    {
        k=BigInteger::randomNumber(q.bitLength());
        if(k.compareTo(d.q)<0)
        {
            struct sig si=sign(s,k,d);
            if(si.delta.isZero()||si.gamma.isZero())
                continue;
            printf("delta:\t");
            si.delta.show();
            printf("\ngamma:\t");
            si.gamma.show();
            ans=verify(s,si,d);
            break;
        }
    }
    t=clock()-t;
    printf("\nVerification:\t%s\n",ans?"true":"false");
    printf("\nTime elpased:\t%lfms\n\n",1000.0*t/CLOCKS_PER_SEC);
    getchar();
    return 0;
}

struct sig sign(const char* s ,const BigInteger&k,const data&d)
{
    struct sig ans{BigInteger::ZERO,BigInteger::ZERO};
    ans.gamma=d.alpha.modPow(k,d.p).mod(d.q);
    uint8_t hash[32];
    sha256(s,strlen(s),hash);
    BigInteger sha=hash2BigInteger(hash);
    BigInteger k_inv=k.modInverse(d.q);
    ans.delta=sha.add(d.a.multiply(ans.gamma)).multiply(k_inv);
    return ans;
}

struct data generateData(const BigInteger&p,const BigInteger&q,const BigInteger&g)
{
    struct data d={p,q,g,BigInteger::ZERO,BigInteger::ZERO};
    while(1)
    {
        d.a=BigInteger::randomNumber(q.bitLength());
        if(d.a.compareTo(q)<0)
            break;
    }
    d.beta=d.alpha.modPow(d.a,p);
    return d;
}

BigInteger hash2BigInteger(const uint8_t* hash)
{
    uint8_t temp[32];
    uint32_t size=sizeof(uint8_t)<<2,t;
    for(int i=0;i<32;i+=4)
    {
        memcpy(&t,hash+28-i,size);
        t=swap_uint32(t);
        memcpy(temp+i,&t,size);
    }
    return BigInteger(temp,32);
}

unsigned swap_uint32(const unsigned &x)
{
    unsigned a = x;
    a = ((a&(0x0000FFFF)) << 16) | ((a&(0xFFFF0000)) >> 16);
    a = ((a&(0x00FF00FF)) << 8) | ((a&(0xFF00FF00)) >> 8);
    return a;
}

bool verify(const char*s,struct sig& si,const data&d)
{
    uint8_t hash[32];
    sha256(s,strlen(s),hash);
    BigInteger sha=hash2BigInteger(hash);
    BigInteger delta_inv=si.delta.modInverse(d.q);
    BigInteger e1=sha.multiply(delta_inv);
    BigInteger e2=si.gamma.multiply(delta_inv);
    BigInteger alpha_e1=d.alpha.modPow(e1,d.p);
    printf("\ne1:\t");
    e1.show();
    printf("\ne2:\t");
    e2.show();
    BigInteger beta_e2=d.beta.modPow(e2,d.p);
    BigInteger ans= alpha_e1.multiply(beta_e2).mod(d.p).mod(d.q);
    printf("\n(alpha^e1)(beta^e2):\t");
    ans.show();
    return ans.compareTo(si.gamma)==0;
}
