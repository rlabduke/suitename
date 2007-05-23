/*               suitenout.h                  */

#ifdef  SUITENOUT
#undef  SUITENOUT
#define SUITENOUT
#else
#define SUITENOUT extern
#endif

SUITENOUT int Ltriage,Loutlier; /*diagnostic flags*/

SUITENOUT int L33out,L32out,L23out,L22out,Ltriageout;           /*groups*/
SUITENOUT int binout[MAXBINS]; /*incl 1--12 named bins */       /*subgroups*/
SUITENOUT int clusterout[MAXBINS][MAXCLST];                     /*lists*/
SUITENOUT double suitenesssum[MAXBINS][MAXCLST],     suitenesssumall;
SUITENOUT int    suitenesscnt[MAXBINS][MAXCLST][12]; 
            /*3rd index for suiteness intervals: 10 at 10ths + 2 extras   */
            /*   extra at zero counts valid suites with suiteness == 0    */
            /*        e.g. 4D distance < 1 but 7D distance > 1            */
            /*   extra at 11   counts valid suites but triaged or outlier */
SUITENOUT int    binnedsuitecountall, reportcountall;
SUITENOUT int    triagecountall;  /*070328*/
SUITENOUT int    Lwannabeout;  /*wannabe in output flag 070429*/

SUITENOUT char temps[256];

static int nout=1; /*Lsringout residue/suite counter*/

/*prototypes*/
void writeoutput(void);
void binstuffout(int, int);
int  transferout(char*);
void kinemageheader(char*); /*070328 char* */
void kinemagestuffer(char* kinemagestuff[]);  /*070414, 070421*/
void  writesuite(int, int, char*, float, float, char*, char*);
void  clearbinout(void);
void  clearclusterout(void);
void  usageout(void);
void  suitenessaverage(int);

