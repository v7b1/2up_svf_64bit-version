//
// 2up_svf~.c
// state-variable filter module for Max/MSP. 
// based on original SVF algorithm by Hal Chamberlin. 
// oversampled and with other tweaks for musical sound.
// by Randy Jones rej@2uptech.com.
// 
// who	when 		what
// rej	11/22/00	created
// rej	7/29/03		Carbon OS X compatible, fast sin approx, improved math


#ifdef __cplusplus		
	extern "C" {
#endif


#ifndef __max_ext_h
#define __max_ext_h
	#include "ext.h"
#endif

#ifndef __max_ext_wind_h
#define __max_ext_wind_h
	#include "ext_wind.h"
#endif

#ifndef __max_z_dsp_h
#define __max_z_dsp_h
	#include "z_dsp.h"
#endif

#ifdef __cplusplus
	} /* end extern "C" */
#endif


#include <string.h>
#include <math.h>


void *svf_class;

typedef struct _svf
{
    t_pxobject x_obj;
    t_float freq;
    t_float q;
    int mode;
    long freq_connected;	// is freq. lead connected to signal?
    long q_connected;		// how about q?
    t_float b0;
    t_float b1;
    t_float b2;
	t_float fs;
	t_float one_over_fs;
	t_float qc;
} t_svf;

void *svf_new(Symbol *s, int ac, Atom *av);
t_int *offset_perform(t_int *w);
t_int *svf_perform(t_int *w);
void svf_float(t_svf *x, double f);
void svf_freq(t_svf *x, double f);
void svf_q(t_svf *x, double f);
void svf_mode(t_svf *x, long n);
void svf_clear(t_svf *x);
void svf_dsp(t_svf *x, t_signal **sp, short *count);
void svf_assist(t_svf *x, void *b, long m, long a, char *s);

int main(void)
{
    setup((struct messlist **)&svf_class, (method)svf_new, (method)dsp_free, (short)sizeof(t_svf), 0L, A_GIMME, 0);
    addmess((method)svf_dsp, "dsp", A_CANT, 0);
   
    addint((method)svf_mode);
    addfloat((method)svf_float);
    
    addmess((method)svf_assist,"assist",A_CANT,0);
    addmess((method)svf_clear,"clear",0);
    dsp_initclass();
    rescopy('STR#',2241);
	return(0);
}

void svf_assist(t_svf *x, void *b, long m, long a, char *s)
{
	assist_string(2241,m,a,1,4,s);
}

void *svf_new(Symbol *s, int ac, Atom *av)
{
    t_svf *x = (t_svf *)newobject(svf_class);
   	dsp_setup((t_pxobject *)x,3);    
    outlet_new((t_pxobject *)x, "signal");
    
    // defaults
    x->freq = 0.;
    x->q = 0.0;
    x->mode = 0;
    svf_clear(x);
    
    if (ac > 0)
    {
    	if(av[0].a_type == A_FLOAT)
    		x->freq = av[0].a_w.w_float;
    	else if(av[0].a_type == A_LONG)
    		x->freq = av[0].a_w.w_long;
    	else
    		post("2up_svf: expected number for param 1 (frequency)");
    }
    if (ac > 1)
    {
    	if(av[1].a_type == A_FLOAT)
    		x->q = av[1].a_w.w_float;
    	else
    		post("2up_svf: expected float for param 2 (resonance)");
    }
   	
    return (x);
}

void svf_float(t_svf *x, double f)
{
	int inlet = ((t_pxobject*)x)->z_in;
	
	if (inlet == 1)
	{
		x->freq = f;
	}
	else if (inlet == 2)
	{
		x->q = f; 
	}
}

void svf_freq(t_svf *x, double f)
{
	x->freq = f;
}

void svf_q(t_svf *x, double f)
{
	x->q = f;
}

void svf_mode(t_svf *x, long n)
{
	x->mode = n;
}

void svf_clear(t_svf *x) 
{
    x->b0 = 0.;
    x->b1 = 0.;
    x->b2 = 0.;
}


void svf_dsp(t_svf *x, t_signal **sp, short *count)
{
	x->fs = sp[0]->s_sr;
	x->one_over_fs = 1.0/x->fs;
	x->freq_connected = count[1];
	x->q_connected = count[2];
	
	// clear filter history
	x->b0 = x->b1 = x->b2 = 0.;
	
	dsp_add(svf_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

inline float sin_e(float x);
inline float sin_e(float x)
{
    return x - (x*x*x*0.1666666666666);
}

void calc_params(t_svf *x, t_int *w, int offset, float * pm_f, float * p_q);
void calc_params(t_svf *x, t_int *w, int offset, float * pm_f, float * p_q)
{
	float lq, fs, fc;
	float oversample = 1.0 / 2.0; // 1/2 = 2x
	fs      = x->fs;
	fc      = x->freq_connected ? ((float *)w[3])[offset] : x->freq;
	lq 		= (x->q_connected ? ((float *)w[4])[offset] : x->q);
	if (fc > fs*0.5) fc = fs*0.5;
	if (fc < 0) fc = 0;
	lq = 0.99 - lq;
	if (lq > 0.7) lq = 0.7;	// by experiment, prevent blowups at high fc and low q
	if (lq < -0.01) lq = -0.01; 
	*p_q = lq;
    *pm_f     = 2.0 * sin_e(PI*fc*x->one_over_fs * oversample);	
    
	// coeff for one-pole
	x->qc = (lq + 0.02)*40000.*x->one_over_fs;
}

t_int *svf_perform(t_int *w)
{
 	t_svf *x = (t_svf *)(w[1]);
   	t_float *p_in = (t_float *)(w[2]);
    t_float *p_out = (t_float *)(w[5]);
	long n = w[6];
	int i = 0;
	int do_signal_params = 0;

	// tiny bit of noise added to prevent denormals
	static const float denorm_noise[3] = {1e-10, 0, -1e-10};
	long nk = 0;

	register t_float in, m_in, wc; 
	register t_float m_h, m_b, m_l; 
	t_float m_f, q;
	
	if (x->x_obj.z_disabled)
		goto bail;
	
	do_signal_params = x->freq_connected || x->q_connected;
	if (!do_signal_params)
    	calc_params(x, w, 0, &m_f, &q);

    // restore history
 	m_b = x->b0;
	m_l = x->b1;
	m_in = x->b2;
	
	// intermediates
	wc = 4.*x->one_over_fs;

    switch (x->mode)
    {
		case 0: // lowpass
			while(n--)
			{	
				in = (p_in[i]);	 

				if (do_signal_params && (0 == (i & 0xF)))
			    	calc_params(x, w, i, &m_f, &q);
			    	
				// one pole LP filter input, based on q
				m_in = m_in*(1. - x->qc) + in*(x->qc);
				
				// regular state variable code here
				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				// repeat for stability (2x oversampling).
				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;
				
				m_l += denorm_noise[nk++];
				if (nk >= 3) nk=0;

				// waveshape
				m_b = m_b - m_b*m_b*m_b*wc;
				
				p_out[i] = m_l;

				i++;
			}
			break;
		case 1: // highpass
			while(n--)
			{	
				in = (p_in[i]);	

				if (do_signal_params && (0 == (i & 0xF)))
			    	calc_params(x, w, i, &m_f, &q);
			    	
				m_in = m_in*(1. - x->qc) + in*(x->qc);

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l += denorm_noise[nk++];
				if (nk >= 3) nk=0;

				// waveshape
				m_b = m_b - m_b*m_b*m_b*wc;
				
				p_out[i] = m_h;
				i++;
			}
			break;
		case 2: // bandpass
			while(n--)
			{	
				in = (p_in[i]);	

				if (do_signal_params && (0 == (i & 0xF)))
			    	calc_params(x, w, i, &m_f, &q);
			    	
				m_in = m_in*(1. - x->qc) + in*(x->qc);

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l += denorm_noise[nk++];
				if (nk >= 3) nk=0;

				// waveshape
				m_b = m_b - m_b*m_b*m_b*wc;
				
				p_out[i] = m_b;
				i++;
			}
			break;
		case 3: // notch
			while(n--)
			{	
				in = (p_in[i]);

				if (do_signal_params && (0 == (i & 0xF)))
			    	calc_params(x, w, i, &m_f, &q);
			    	
				m_in = m_in*(1. - x->qc) + in*(x->qc);

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l = m_l + (m_f * m_b);
				m_h = m_in - m_l - q * m_b;
				m_b = m_b + m_f * m_h;

				m_l += denorm_noise[nk++];
				if (nk >= 3) nk=0;

				// waveshape
				m_b = m_b - m_b*m_b*m_b*wc;

				p_out[i] = m_h + m_l;
				i++;
			}
			break;
	}
	
   
    // save history
    x->b0 = m_b;
    x->b1 = m_l;
    x->b2 = m_in;

bail:
    return (w+7);
}

