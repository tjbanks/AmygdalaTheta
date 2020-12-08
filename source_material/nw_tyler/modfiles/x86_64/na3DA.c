/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__na3DA
#define _nrn_initial _nrn_initial__na3DA
#define nrn_cur _nrn_cur__na3DA
#define _nrn_current _nrn_current__na3DA
#define nrn_jacob _nrn_jacob__na3DA
#define nrn_state _nrn_state__na3DA
#define _net_receive _net_receive__na3DA 
#define states states__na3DA 
#define trates trates__na3DA 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define ar2 _p[1]
#define m _p[2]
#define h _p[3]
#define s _p[4]
#define ena _p[5]
#define ina _p[6]
#define thegna _p[7]
#define tha1 _p[8]
#define Dm _p[9]
#define Dh _p[10]
#define Ds _p[11]
#define _g _p[12]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_DA2(void);
 static void _hoc_DA1(void);
 static void _hoc_alps(void);
 static void _hoc_alpv(void);
 static void _hoc_bets(void);
 static void _hoc_trap0(void);
 static void _hoc_trates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_na3DA", _hoc_setdata,
 "DA2_na3DA", _hoc_DA2,
 "DA1_na3DA", _hoc_DA1,
 "alps_na3DA", _hoc_alps,
 "alpv_na3DA", _hoc_alpv,
 "bets_na3DA", _hoc_bets,
 "trap0_na3DA", _hoc_trap0,
 "trates_na3DA", _hoc_trates,
 0, 0
};
#define DA2 DA2_na3DA
#define DA1 DA1_na3DA
#define alps alps_na3DA
#define alpv alpv_na3DA
#define bets bets_na3DA
#define trap0 trap0_na3DA
 extern double DA2( double );
 extern double DA1( double );
 extern double alps( double );
 extern double alpv( double );
 extern double bets( double );
 extern double trap0( double , double , double , double );
 /* declare global and static user variables */
#define DA_t2 DA_t2_na3DA
 double DA_t2 = -1;
#define DA_start2 DA_start2_na3DA
 double DA_start2 = 36000;
#define DA_period2 DA_period2_na3DA
 double DA_period2 = 100;
#define DA_t1 DA_t1_na3DA
 double DA_t1 = -0.1;
#define DA_ext2 DA_ext2_na3DA
 double DA_ext2 = 212000;
#define DA_ext1 DA_ext1_na3DA
 double DA_ext1 = 196000;
#define DA_stop DA_stop_na3DA
 double DA_stop = 96000;
#define DA_start DA_start_na3DA
 double DA_start = 64000;
#define DA_period DA_period_na3DA
 double DA_period = 500;
#define Rd Rd_na3DA
 double Rd = 0.03;
#define Rg Rg_na3DA
 double Rg = 0.01;
#define Rb Rb_na3DA
 double Rb = 0.124;
#define Ra Ra_na3DA
 double Ra = 0.4;
#define a0s a0s_na3DA
 double a0s = 0.0003;
#define gms gms_na3DA
 double gms = 0.2;
#define hmin hmin_na3DA
 double hmin = 0.5;
#define htau htau_na3DA
 double htau = 0;
#define hinf hinf_na3DA
 double hinf = 0;
#define mmin mmin_na3DA
 double mmin = 0.02;
#define mtau mtau_na3DA
 double mtau = 0;
#define minf minf_na3DA
 double minf = 0;
#define qq qq_na3DA
 double qq = 10;
#define q10 q10_na3DA
 double q10 = 2;
#define qg qg_na3DA
 double qg = 1.5;
#define qd qd_na3DA
 double qd = 1.5;
#define qa qa_na3DA
 double qa = 7.2;
#define qinf qinf_na3DA
 double qinf = 4;
#define smax smax_na3DA
 double smax = 10;
#define sinf sinf_na3DA
 double sinf = 0;
#define tq tq_na3DA
 double tq = -55;
#define thi2 thi2_na3DA
 double thi2 = -45;
#define thi1 thi1_na3DA
 double thi1 = -45;
#define tha tha_na3DA
 double tha = -30;
#define tone_period tone_period_na3DA
 double tone_period = 4000;
#define thinf thinf_na3DA
 double thinf = -50;
#define taus taus_na3DA
 double taus = 0;
#define vvs vvs_na3DA
 double vvs = 2;
#define vvh vvh_na3DA
 double vvh = -58;
#define vhalfs vhalfs_na3DA
 double vhalfs = -60;
#define zetas zetas_na3DA
 double zetas = 12;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tha_na3DA", "mV",
 "qa_na3DA", "mV",
 "Ra_na3DA", "/ms",
 "Rb_na3DA", "/ms",
 "thi1_na3DA", "mV",
 "thi2_na3DA", "mV",
 "qd_na3DA", "mV",
 "qg_na3DA", "mV",
 "Rg_na3DA", "/ms",
 "Rd_na3DA", "/ms",
 "qq_na3DA", "mV",
 "tq_na3DA", "mV",
 "thinf_na3DA", "mV",
 "qinf_na3DA", "mV",
 "vhalfs_na3DA", "mV",
 "a0s_na3DA", "ms",
 "zetas_na3DA", "1",
 "gms_na3DA", "1",
 "smax_na3DA", "ms",
 "vvh_na3DA", "mV",
 "vvs_na3DA", "mV",
 "mtau_na3DA", "ms",
 "htau_na3DA", "ms",
 "sinf_na3DA", "ms",
 "taus_na3DA", "ms",
 "gbar_na3DA", "mho/cm2",
 "ar2_na3DA", "1",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tone_period_na3DA", &tone_period_na3DA,
 "DA_period_na3DA", &DA_period_na3DA,
 "DA_start_na3DA", &DA_start_na3DA,
 "DA_stop_na3DA", &DA_stop_na3DA,
 "DA_ext1_na3DA", &DA_ext1_na3DA,
 "DA_ext2_na3DA", &DA_ext2_na3DA,
 "DA_t1_na3DA", &DA_t1_na3DA,
 "DA_period2_na3DA", &DA_period2_na3DA,
 "DA_start2_na3DA", &DA_start2_na3DA,
 "DA_t2_na3DA", &DA_t2_na3DA,
 "tha_na3DA", &tha_na3DA,
 "qa_na3DA", &qa_na3DA,
 "Ra_na3DA", &Ra_na3DA,
 "Rb_na3DA", &Rb_na3DA,
 "thi1_na3DA", &thi1_na3DA,
 "thi2_na3DA", &thi2_na3DA,
 "qd_na3DA", &qd_na3DA,
 "qg_na3DA", &qg_na3DA,
 "mmin_na3DA", &mmin_na3DA,
 "hmin_na3DA", &hmin_na3DA,
 "q10_na3DA", &q10_na3DA,
 "Rg_na3DA", &Rg_na3DA,
 "Rd_na3DA", &Rd_na3DA,
 "qq_na3DA", &qq_na3DA,
 "tq_na3DA", &tq_na3DA,
 "thinf_na3DA", &thinf_na3DA,
 "qinf_na3DA", &qinf_na3DA,
 "vhalfs_na3DA", &vhalfs_na3DA,
 "a0s_na3DA", &a0s_na3DA,
 "zetas_na3DA", &zetas_na3DA,
 "gms_na3DA", &gms_na3DA,
 "smax_na3DA", &smax_na3DA,
 "vvh_na3DA", &vvh_na3DA,
 "vvs_na3DA", &vvs_na3DA,
 "minf_na3DA", &minf_na3DA,
 "hinf_na3DA", &hinf_na3DA,
 "mtau_na3DA", &mtau_na3DA,
 "htau_na3DA", &htau_na3DA,
 "sinf_na3DA", &sinf_na3DA,
 "taus_na3DA", &taus_na3DA,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"na3DA",
 "gbar_na3DA",
 "ar2_na3DA",
 0,
 0,
 "m_na3DA",
 "h_na3DA",
 "s_na3DA",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gbar = 0.01;
 	ar2 = 1;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _na3DA_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 na3DA /home/tbg28/git_stage/nw_tyler/modfiles/x86_64/na3DA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double _zmexp , _zhexp , _zsexp ;
static int _reset;
static char *modelname = "na3";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int trates(double, double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
double alpv (  double _lv ) {
   double _lalpv;
 _lalpv = 1.0 / ( 1.0 + exp ( ( _lv - vvh ) / vvs ) ) ;
   
return _lalpv;
 }
 
static void _hoc_alpv(void) {
  double _r;
   _r =  alpv (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alps (  double _lv ) {
   double _lalps;
 _lalps = exp ( 1.e-3 * zetas * ( _lv - vhalfs ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lalps;
 }
 
static void _hoc_alps(void) {
  double _r;
   _r =  alps (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double bets (  double _lv ) {
   double _lbets;
 _lbets = exp ( 1.e-3 * zetas * gms * ( _lv - vhalfs ) * 9.648e4 / ( 8.315 * ( 273.16 + celsius ) ) ) ;
   
return _lbets;
 }
 
static void _hoc_bets(void) {
  double _r;
   _r =  bets (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Ds = ( sinf - s ) / taus ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v , ar2 ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taus )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taus)))*(- ( ( ( sinf ) ) / taus ) / ( ( ( ( - 1.0 ) ) ) / taus ) - s) ;
   }
  return 0;
}
 
static int  trates (  double _lvm , double _la2 ) {
   double _la , _lb , _lc , _lqt ;
 _lqt = pow( q10 , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   tha1 = tha + DA1 ( _threadargscomma_ t ) + DA2 ( _threadargscomma_ t ) ;
   _la = trap0 ( _threadargscomma_ _lvm , tha1 , Ra , qa ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - tha1 , Rb , qa ) ;
   mtau = 1.0 / ( _la + _lb ) / _lqt ;
   if ( mtau < mmin ) {
     mtau = mmin ;
     }
   if ( v < - 57.5 ) {
     minf = 0.0 ;
     }
   else {
     minf = _la / ( _la + _lb ) ;
     }
   _la = trap0 ( _threadargscomma_ _lvm , thi1 , Rd , qd ) ;
   _lb = trap0 ( _threadargscomma_ - _lvm , - thi2 , Rg , qg ) ;
   htau = 1.0 / ( _la + _lb ) / _lqt ;
   if ( htau < hmin ) {
     htau = hmin ;
     }
   hinf = 1.0 / ( 1.0 + exp ( ( _lvm - thinf ) / qinf ) ) ;
   _lc = alpv ( _threadargscomma_ _lvm ) ;
   sinf = _lc + _la2 * ( 1.0 - _lc ) ;
   taus = bets ( _threadargscomma_ _lvm ) / ( a0s * ( 1.0 + alps ( _threadargscomma_ _lvm ) ) ) ;
   if ( taus < smax ) {
     taus = smax ;
     }
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   _r = 1.;
 trates (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double trap0 (  double _lv , double _lth , double _la , double _lq ) {
   double _ltrap0;
 if ( fabs ( _lv - _lth ) > 1e-6 ) {
     _ltrap0 = _la * ( _lv - _lth ) / ( 1.0 - exp ( - ( _lv - _lth ) / _lq ) ) ;
     }
   else {
     _ltrap0 = _la * _lq ;
     }
   
return _ltrap0;
 }
 
static void _hoc_trap0(void) {
  double _r;
   _r =  trap0 (  *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 hoc_retpushx(_r);
}
 
double DA1 (  double _lt ) {
   double _lDA1;
 if ( _lt >= DA_start  && _lt <= DA_stop ) {
     if ( ( _lt / tone_period - floor ( _lt / tone_period ) ) >= ( 1.0 - DA_period / tone_period ) ) {
       _lDA1 = DA_t1 ;
       }
     else if ( ( _lt / tone_period - floor ( _lt / tone_period ) )  == 0.0 ) {
       _lDA1 = DA_t1 ;
       }
     else {
       _lDA1 = 0.0 ;
       }
     }
   else if ( _lt >= DA_ext1  && _lt <= DA_ext2 ) {
     if ( ( _lt / tone_period - floor ( _lt / tone_period ) ) >= ( 1.0 - DA_period / tone_period ) ) {
       _lDA1 = DA_t1 ;
       }
     else if ( ( _lt / tone_period - floor ( _lt / tone_period ) )  == 0.0 ) {
       _lDA1 = DA_t1 ;
       }
     else {
       _lDA1 = 0.0 ;
       }
     }
   else {
     _lDA1 = 0.0 ;
     }
   
return _lDA1;
 }
 
static void _hoc_DA1(void) {
  double _r;
   _r =  DA1 (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double DA2 (  double _lt ) {
   double _lDA2;
 if ( _lt >= DA_start2  && _lt <= DA_stop ) {
     if ( ( _lt / tone_period - floor ( _lt / tone_period ) ) >= ( 1.0 - DA_period2 / tone_period ) ) {
       _lDA2 = DA_t2 ;
       }
     else if ( ( _lt / tone_period - floor ( _lt / tone_period ) )  == 0.0 ) {
       _lDA2 = DA_t2 ;
       }
     else {
       _lDA2 = 0.0 ;
       }
     }
   else {
     _lDA2 = 0.0 ;
     }
   
return _lDA2;
 }
 
static void _hoc_DA2(void) {
  double _r;
   _r =  DA2 (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  s = s0;
 {
   trates ( _threadargscomma_ v , ar2 ) ;
   m = minf ;
   h = hinf ;
   s = sinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   thegna = gbar * m * m * m * h * s ;
   ina = thegna * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
 _g = _nrn_current(_v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
 { error =  states();
 if(error){fprintf(stderr,"at line 82 in file na3DA.mod:\n        SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(s) - _p;  _dlist1[2] = &(Ds) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/tbg28/git_stage/nw_tyler/modfiles/na3DA.mod";
static const char* nmodl_file_text = 
  "TITLE na3\n"
  ": Na current \n"
  ": from Jeff M.\n"
  ":  ---------- modified -------M.Migliore may97\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX na3DA\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE  gbar, ar2\n"
  "	GLOBAL minf, hinf, mtau, htau, sinf, taus,qinf, thinf\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tone_period = 4000   \n"
  "	DA_period = 500\n"
  "	DA_start = 64000		    : D1R(Low Affinity) Dopamine Effect after 6 conditioning trials (14*4000) = 64000)\n"
  "	DA_stop = 96000\n"
  "	DA_ext1 = 196000\n"
  "	DA_ext2 = 212000\n"
  "	DA_t1 = -0.1 : -0.3 : -0.15            : Amount of DA effect- negative value decreases AP threshold / positive value increases threshold of AP\n"
  "\n"
  "	DA_period2 = 100\n"
  "	DA_start2 = 36000		   : shock Dopamine Effect during shock after 1 conditioning trial\n"
  "	DA_t2 = -1 : -1 : -0.1           : Amount(%) of DA effect- negative value decreases AP threshold / positive value increases threshold of AP	\n"
  "	\n"
  "	gbar = 0.010   	(mho/cm2)	\n"
  "								\n"
  "	tha  =  -30	(mV)		: v 1/2 for act	\n"
  "	qa   = 7.2	(mV)		: act slope (4.5)		\n"
  "	Ra   = 0.4	(/ms)		: open (v)		\n"
  "	Rb   = 0.124 	(/ms)		: close (v)		\n"
  "\n"
  "	thi1  = -45	(mV)		: v 1/2 for inact 	\n"
  "	thi2  = -45 	(mV)		: v 1/2 for inact 	\n"
  "	qd   = 1.5	(mV)	        : inact tau slope\n"
  "	qg   = 1.5      (mV)\n"
  "	mmin=0.02	\n"
  "	hmin=0.5			\n"
  "	q10=2\n"
  "	Rg   = 0.01 	(/ms)		: inact recov (v) 	\n"
  "	Rd   = .03 	(/ms)		: inact (v)	\n"
  "	qq   = 10        (mV)\n"
  "	tq   = -55      (mV)\n"
  "\n"
  "	thinf  = -50 	(mV)		: inact inf slope	\n"
  "	qinf  = 4 	(mV)		: inact inf slope \n"
  "\n"
  "        vhalfs=-60	(mV)		: slow inact.\n"
  "        a0s=0.0003	(ms)		: a0s=b0s\n"
  "        zetas=12	(1)\n"
  "        gms=0.2		(1)\n"
  "        smax=10		(ms)\n"
  "        vvh=-58		(mV) \n"
  "        vvs=2		(mV)\n"
  "        ar2=1		(1)		: 1=no inact., 0=max inact.\n"
  "	ena		(mV)            : must be explicitly def. in hoc\n"
  "	celsius\n"
  "	v 		(mV)\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "	ina 		(mA/cm2)\n"
  "	thegna		(mho/cm2)\n"
  "	minf 		hinf 		\n"
  "	mtau (ms)	htau (ms) 	\n"
  "	sinf (ms)	taus (ms)\n"
  "	tha1	\n"
  "}\n"
  " \n"
  "\n"
  "STATE { m h s}\n"
  "\n"
  "BREAKPOINT {\n"
  "        SOLVE states METHOD cnexp\n"
  "        thegna = gbar*m*m*m*h*s\n"
  "	ina = thegna * (v - ena)\n"
  "} \n"
  "\n"
  "INITIAL {\n"
  "	trates(v,ar2)\n"
  "	m=minf  \n"
  "	h=hinf\n"
  "	s=sinf\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION alpv(v(mV)) {\n"
  "         alpv = 1/(1+exp((v-vvh)/vvs))\n"
  "}\n"
  "        \n"
  "FUNCTION alps(v(mV)) {  \n"
  "  alps = exp(1.e-3*zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))\n"
  "}\n"
  "\n"
  "FUNCTION bets(v(mV)) {\n"
  "  bets = exp(1.e-3*zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))\n"
  "}\n"
  "\n"
  "LOCAL mexp, hexp, sexp\n"
  "\n"
  "DERIVATIVE states {   \n"
  "        trates(v,ar2)      \n"
  "        m' = (minf-m)/mtau\n"
  "        h' = (hinf-h)/htau\n"
  "        s' = (sinf - s)/taus\n"
  "}\n"
  "\n"
  "PROCEDURE trates(vm,a2) {  \n"
  "        LOCAL  a, b, c, qt\n"
  "        qt=q10^((celsius-24)/10)\n"
  "		tha1 = tha + DA1(t)	+ DA2(t)\n"
  "	a = trap0(vm,tha1,Ra,qa)\n"
  "	b = trap0(-vm,-tha1,Rb,qa)\n"
  "	mtau = 1/(a+b)/qt\n"
  "        if (mtau<mmin) {mtau=mmin}\n"
  "		\n"
  "	if (v < -57.5 ) {\n"
  "	minf = 0\n"
  "	} else{\n"
  "	minf = a/(a+b)\n"
  "	}\n"
  "	\n"
  "	a = trap0(vm,thi1,Rd,qd)\n"
  "	b = trap0(-vm,-thi2,Rg,qg)\n"
  "	htau =  1/(a+b)/qt\n"
  "        if (htau<hmin) {htau=hmin}\n"
  "	hinf = 1/(1+exp((vm-thinf)/qinf))\n"
  "	c=alpv(vm)\n"
  "        sinf = c+a2*(1-c)\n"
  "        taus = bets(vm)/(a0s*(1+alps(vm)))\n"
  "        if (taus<smax) {taus=smax}\n"
  "}\n"
  "\n"
  "FUNCTION trap0(v,th,a,q) {\n"
  "	if (fabs(v-th) > 1e-6) {\n"
  "	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))\n"
  "	} else {\n"
  "	        trap0 = a * q\n"
  " 	}\n"
  "}	\n"
  "FUNCTION DA1(t) {\n"
  "	    if (t >= DA_start && t <= DA_stop){ 									: During conditioning\n"
  "			if ((t/tone_period-floor(t/tone_period)) >= (1-DA_period/tone_period)) {DA1 = DA_t1}\n"
  "			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA1 = DA_t1}\n"
  "			else {DA1 = 0}}\n"
  "		else if (t >= DA_ext1 && t <= DA_ext2){								: During 4trials of Extinction\n"
  "			if ((t/tone_period-floor(t/tone_period)) >= (1-DA_period/tone_period)) {DA1 = DA_t1}\n"
  "			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA1 = DA_t1}\n"
  "			else {DA1 = 0}}		\n"
  "		else  {DA1 = 0}\n"
  "	}\n"
  "FUNCTION DA2(t) {\n"
  "	    if (t >= DA_start2 && t <= DA_stop){\n"
  "			if((t/tone_period-floor(t/tone_period)) >= (1-DA_period2/tone_period)) {DA2 = DA_t2}\n"
  "			else if ((t/tone_period-floor(t/tone_period)) == 0) {DA2 = DA_t2}\n"
  "			else  {DA2 = 0}}\n"
  "		else  {DA2 = 0}\n"
  "	}\n"
  ;
#endif
