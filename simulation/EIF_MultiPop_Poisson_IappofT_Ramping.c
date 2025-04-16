/* To use function from matlab, first compile by entering this into the Matlab command window:
* mex EIF1DRF2tausyn.c
* Then call the function like this:
* [s,IsynX0,IsynE0,IsynI0,v]=EIF1DRFfastslowSyn(sx, Wrf,Wrr,param);
*/

/* sx is the feedforward presynaptic spike trains, which should be 3xNsx where Nsx is the number of spikes.
* sx(1,:) is spike times, sx(2,:) is the index of the neuron (from 1 to Nx)
* srr is the same as sx except it is poisson input from recurrent layer
* Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices,
* sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
* cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population.
* For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and
* Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop.
* Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
* The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.
*
* param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr,
* gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
* maxns, Irecord, Psyn
* Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
* Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii];
* taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
* taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
* Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
*
*
* Kab is the number of projections from each cell in pop b=e,i,x to all cells in pop a=e,i
* Ne, Ni are the number of excitatory, inhibitory cells, N=Ne+Ni.
* Nx is the number of neurons in feedforward layer.
* Jab is the synaptic strength of connections from b=e,i,x to a=e,i.
*
* C,gl,vl,DeltaT,VT,tref,Vth,Vre,Vlb are EIF neuron params
* They are each 2x1 vectors, for exc and inh neurons separately.
* For example, C(1) is the capacitance of exc neurons and C(2) of inh neurons
* tausynb is the time-constant of the synapses from neurons in population b=x,e,i.
* post-synaptic currents are of the form (1/tausynb)*exp(-t/tausynb) where t>0 is
* the time evolved since the presynaptic spike.
* V0 is vector of all membrane potential initial conditions.
* It should be Nx1 where N=Ne+Ni is the number of neurons in the recurrent network
* The first Ne elements are for exc neurons, the last Ni for inh neurons
* dt is bin size for time
* maxns is maximum number of spikes allowed for all neurons together
* Irecord is a 1x(Nrecord) matrix indicating for which neurons we should record
* the synaptic inputs and membrane potential. The first Ne elements are for exc neurons,
* the last Ni for inh neurons
*
*
* Outputs:
* s is a 2x(maxns) matrix of spikes
* s(1,:) contains spike times.
* s(2,:) contains indices of neurons that spike
* When there are fewer than maxns spikes, extra space in s will be filled
* with zeros.  This should be truncated in matlab by writing s=s(:,s(1,:)>0);
* IsynX0,IsynE0,IsynI0,v are the recorded
* feedforward synaptic input, exc synaptic input, inh synaptic input and voltage
* respectively.
*
*/



#include "mex.h"
#include "math.h"
#include "time.h"
#include "matrix.h"


/* A fast approximation of the exp function */
static union
{
    double d;
    struct {
        #ifdef LITTLE_ENDIAN
        int j,i;
        #else
        int i,j;
        #endif
    } n;
} _eco;
#define EXP_A (1048576/0.69314718055994530942)
#define EXP_C 60801
#define EXP(y) (_eco.n.i = EXP_A*(y) + (1072693248 - EXP_C), _eco.d)


//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    //int *Ntref,Npop,Nx,*Nsum,*Ksum,Kxsum,Kmax,post_pop,pre_pop,k,j,i,N,Nt,m1,m2,maxns,ns,Nrecord,jj,nskiprecord,tempflag,Nsx,Nw;
    int *Ntref,Npop,Nx,*Nsum,*Ksum,Kxsum,Kmax,post_pop,pre_pop,k,j,i,N,Nt,m1,m2,maxns,ns,Nrecord,jj,kk, nskiprecord,tempflag,Nsx,Nw,Nsrr;
    double dt,*s,*v,*v0,*Isynrecord, *Psyn, *Isynprimerecord;
    double *Isyn, *Isynprime,*Kr,*Kx,*Ncell;
    int *Wrr, *Wrf, *Pr;
    double *taudsyn,*taursyn,*vr,*sx,Jex,Jix,*iXspkInd,*Psyn2, *srr;
    double *Jr,*Jx,*Iapp,T,*Irecord,*C,*Vleak,*DeltaT,*VT,*tref,*gl,*Vth,*Vre,*Vlb,xloc,yloc;
    int *refstate,iXspike,jspike,*postcell, Nsyn, isyn, Nsyntype, *syntype;
    const mxArray  *mxTmp;
    double  *temp1,*temp2,  Isyntot, *Jext, *iRspkInd;
    int rspike, iRspike;



    double *IappScaleOutput,*Vout;
    int popwithinput, modnum, coeffstep;
    double *iappScale;
    int scaleCount = 0;

    double *Isyn_out, *Isynprime_out, *refstate_out;
    double *Isyn_init, *Isynprime_init, *refstate_init;
    double *Isyn_in, *Isynprime_in;


    //double *Isyn_last, *Isynprime_last;
    //double *Isyn_initial, *Isynprime_initial, *Isyn_UnChangeinitial, *Isynprime_UnChangeinitial;;
   // int modnum = 5000;
    //double coeffstep = +0.05; // written with + coeffiecnt in addition, if minus/down step then make coef negative value


    /******
    * Import variables from matlab
    * This is messy looking and is specific to mex.
    * Ignore if you're implementing this outside of mex.
    *******/
    sx = mxGetPr(prhs[0]); //feedforward presynatptic spike trains
    m1 = mxGetM(prhs[0]); // first dimension of sx, should equal 2, spike times and index of neurons
    Nsx = mxGetN(prhs[0]); // second dimension of sx, number of spikes
    if(m1!=2){
        mexErrMsgTxt("sx should be Nsx2"); // error message in MATLB, quits mex and returns to control to MATLAB
    }

    Wrf = (int *)mxGetData(prhs[1]);// vector of connections from feedforward layer to recurrent layer, sorted by index of presynatic cell. The block of postsynaptic cell indices for each presynatic cell is sorted as excitatiry followed by inhib cells
    Nw = mxGetM(prhs[1]); //first dimension of Wrf, number of connections
    m1 = mxGetN(prhs[1]); // second dimension of Wrf, should be 1 since Wrf is a vectors
    if(m1!=1){
        mexErrMsgTxt("Weight matrix Wrf must be Nwx1 where Nw is numer of connections.");
    }

    Wrr = (int *)mxGetData(prhs[2]); //same as Wrf, except for connections amoung recurrent layer, containing postdynaptic cell indicies, sorted by the index of the presynaptic cell, sorted as excitatory cells in block then inhib cells
    Nw = mxGetM(prhs[2]); // first dimension of Wrr, number of connects in layer
    m1 = mxGetN(prhs[2]); // second dimension of Wrr, should be 1
    if(m1!=1){
        mexErrMsgTxt("Weight matrix Wrr must be Nwx1 where Nw is numer of connections.");
    }

    /* check if prhs[3] is a struct */
    if(mxIsStruct(prhs[3])!=1){
        mexErrMsgTxt("The 4th input should be the parameter struct.");
    }

    srr = mxGetPr(prhs[4]); //poisson distributed spk trains in recurrent layer
    m1 = mxGetM(prhs[4]); // first dimension of srr, should equal 2, spike times and index of neurons
    Nsrr = mxGetN(prhs[4]); // second dimension of srr, number of spikes
    if(m1!=2){
        mexErrMsgTxt("srr should be Nsrr x 2"); // error message in MATLB, quits mex and returns to control to MATLAB
    }

    /* Number of neurons in each pop. */
    mxTmp = mxGetField(prhs[3],0,"Ncell"); // mxGetField(pm, index, fieldname): pm - pointer to a structure mxArry, index is index of desired element, fieldname is name of field whose value you want to extract
    Ncell=mxGetPr(mxTmp); // Number of cells per population in reccurent layer
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    Npop = m1*m2;

    /*    for(i=0;i<Npop; i++){
    mexPrintf("i=%d, Ncell=%d \n",i,(int)Ncell[i]);
    } // Nsum results as vector:
    */
    /* Total number of neurons  */
    mxTmp = mxGetField(prhs[3],0,"N");
    N=(int)mxGetPr(mxTmp)[0]; // total number of neurons in simulation

    /* Number of neurons in the ffwd layer  */
    mxTmp = mxGetField(prhs[3],0,"Nx");
    Nx=(int)mxGetPr(mxTmp)[0]; //number of neurons in ffwd layer

    mxTmp = mxGetField(prhs[3],0,"Jx");
    Jx=mxGetPr(mxTmp); //feedforward weight (to all neurons???)
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("Jx should be 1xNpop");

    mxTmp = mxGetField(prhs[3],0,"Jext");
    Jext=mxGetPr(mxTmp);
    //mexPrintf("%f %f %f %f",Jext[0],Jext[1],Jext[2],Jext[3]);
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("Jext should be 1xNpop");

    mxTmp = mxGetField(prhs[3],0,"Jr");
    Jr=mxGetPr(mxTmp); //recurrent weight from all neurons???
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1!=Npop||m2!=Npop)
        mexErrMsgTxt("Jr should be Npop x Npop");


    mxTmp = mxGetField(prhs[3],0,"Kx");
    Kx=mxGetPr(mxTmp); //feedfoward projection matrix
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("Kx should be 1xNpop");


    mxTmp = mxGetField(prhs[3],0,"Kr");
    Kr = mxGetPr(mxTmp); //recurrent projection matrix
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1!=Npop||m2!=Npop)
        mexErrMsgTxt("Kr should be Npop x Npop");


    mxTmp = mxGetField(prhs[3],0,"Iapp");
    Iapp=mxGetPr(mxTmp); // constant external input
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop) //check to be m1*m2 = dtstep *Npop
        mexErrMsgTxt("Iapp should be 1xNpop"); //change error message


    mxTmp = mxGetField(prhs[3],0,"Cm");
    C=mxGetPr(mxTmp); //membrane capacitance
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"gl");
    gl=mxGetPr(mxTmp); // leaky conductance
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"vl");
    Vleak=mxGetPr(mxTmp); // resting leaky potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");

    mxTmp = mxGetField(prhs[3],0,"DeltaT");
    DeltaT=mxGetPr(mxTmp); //neuron EIF slope score
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");

    mxTmp = mxGetField(prhs[3],0,"vT");
    VT=mxGetPr(mxTmp); // parameter for EIF exponential part
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"tref");
    tref=mxGetPr(mxTmp); //refractory period
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"vth");
    Vth=mxGetPr(mxTmp); // voltage threshold for spiking
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"vre");
    Vre=mxGetPr(mxTmp); // voltage reset potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");
    mxTmp = mxGetField(prhs[3],0,"vlb");
    Vlb=mxGetPr(mxTmp); // lower boundary of membrane voltage
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=Npop)
        mexErrMsgTxt("All neuron parameters should be 1xNpop");

    mxTmp = mxGetField(prhs[3],0,"taursyn");
    m1=mxGetN(mxTmp);
    m2=mxGetM(mxTmp);
    if(m2!=(Npop+1))
        mexErrMsgTxt("size(taursyn,1) should be Npop+1");
    Nsyntype=m1; // synatpic type
    taursyn=mxGetPr(mxTmp); //synaptic rise time constant

    mxTmp = mxGetField(prhs[3],0,"taudsyn");
    m1=mxGetN(mxTmp);
    m2=mxGetM(mxTmp);
    if(m2!=(Npop+1))
        mexErrMsgTxt("size(taudsyn,1) should be Npop+1");
    if(m1!=Nsyntype)
        mexErrMsgTxt("size(taursyn,1) should equal size(taudsyn,2)");
    taudsyn=mxGetPr(mxTmp); //synaptic decary time constant

    mxTmp = mxGetField(prhs[3],0,"Psyn");
    Psyn=mxGetPr(mxTmp); //percentage of synapse, based on slow/fast spiking profile
    m2=mxGetM(mxTmp);
    if(m2!=(Npop+1))
        mexErrMsgTxt("size(Psyn,1) should be Npop+1");
    m1=mxGetN(mxTmp);
    if(m1!=Nsyntype)
        mexErrMsgTxt("size(Psyn,2) should equal size(taursyn,2)");

    syntype=mxMalloc(m1*m2*sizeof(int));
    Psyn2=mxMalloc(m1*m2*sizeof(double));
    temp1=mxMalloc(m1*m2*sizeof(double)); /* temporary constant */
    temp2=mxMalloc(m1*m2*sizeof(double));
    Nsyn = 0;
    for(isyn=0;isyn<m1*m2; isyn++){
        if(Psyn[isyn]){
            syntype[Nsyn]=isyn%(Npop+1);
            Psyn2[Nsyn]=Psyn[isyn];
            temp1[Nsyn]=(1/taudsyn[isyn]+1/taursyn[isyn]); /* taursyn: syn rise time const, taudsyn: syn decay time const */
            temp2[Nsyn]=1/(taudsyn[isyn]*taursyn[isyn]);
            /*   mexPrintf("Nsyn=%d,syntype[Nsyn]=%d,Psyn2[Nsyn]=%.2f,temp1[Nsyn]=%.2f,temp2[Nsyn]=%.2f  \n",Nsyn,syntype[Nsyn],Psyn2[Nsyn],temp1[Nsyn],temp2[Nsyn]); */
            Nsyn++; }
    }
    mexPrintf("Nsyn:  %d\n", Nsyn);

    mxTmp = mxGetField(prhs[3],0,"V0");
    v0 = mxGetPr(mxTmp); // initial membrane potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=N)
        mexErrMsgTxt("size of V0 should be Nx1");

    mxTmp = mxGetField(prhs[3],0,"Isyn_init");
    Isyn_init = mxGetPr(mxTmp); // initial membrane potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=N*Nsyn)
        mexErrMsgTxt("size of Isynprime_init should be 1x (N* Nsyn)");
    
    mxTmp = mxGetField(prhs[3],0,"Isynprime_init");
    Isynprime_init = mxGetPr(mxTmp); // initial membrane potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    mexPrintf("Isynprime_init:  %d x %d", m1, m2);
    if(m1*m2!=N*Nsyn)
        mexErrMsgTxt("size of Isynprime_init should be 1x (N* Nsyn)");
    
    mxTmp = mxGetField(prhs[3],0,"refstate_init");
    refstate_init = mxGetPr(mxTmp); // initial membrane potential
    m1 = mxGetN(mxTmp);
    m2 = mxGetM(mxTmp);
    if(m1*m2!=N)
        mexErrMsgTxt("size of refstate_init should be Nx1");


    mxTmp = mxGetField(prhs[3],0,"T");
    T = mxGetPr(mxTmp)[0]; // total simulation time
    mxTmp = mxGetField(prhs[3],0,"dt");
    dt =mxGetPr(mxTmp)[0]; // time step

    mxTmp = mxGetField(prhs[3],0,"maxns");
    maxns = (int)mxGetPr(mxTmp)[0]; // max number of spikes allows for allowed for all time

    mxTmp = mxGetField(prhs[3],0,"Irecord");
    Irecord=mxGetPr(mxTmp); // 1x(Nrecord) matrix indicating for which neurons we should record the synaptic inputs and membrane potential. First Ne elements, then Ni elements
    Nrecord = mxGetN(mxTmp); // number of neurons being recorded
    m2 = mxGetM(mxTmp);
    if(m2!=1)
        mexErrMsgTxt("Irecord should be Nrecord x 1.");

    mxTmp = mxGetField(prhs[3],0,"popwithinput");
    popwithinput= mxGetPr(mxTmp)[0]; // total number of neurons in simulation
    mexPrintf("pop input = %d \n",popwithinput);

    /******
    * Finished importing variables.
    *******/

    /* Number of time bins */
    Nt=(int)(T/dt);

    mexPrintf("Npop=%d, N=%d, Nrecord=%d, Nsyn=%d, Nt=%d",Npop,N,Nrecord,Nsyn,Nt);

    mexPrintf("Kr(2,1)=%d, Kr(2,2)=%d, Kr(2,3)=%d, Kr(2,4)=%d  \n",(int)Kr[1],(int)Kr[1+1*Npop],(int)Kr[1+2*Npop],(int)Kr[1+3*Npop]);


    /******
    * Now allocate new variables.
    * This is also mex specific.  Use malloc in C, etc.
    *****/

    /* Allocate output vector */
    // s1
    plhs[0] = mxCreateDoubleMatrix(2, maxns, mxREAL);
    s=mxGetPr(plhs[0]); //2x(m x ns): matrix of spikes, spike times indices of neurons that spike

    // Isyn
    plhs[1] = mxCreateDoubleMatrix(N*Nsyn, 1, mxREAL);
    Isyn_out=mxGetPr(plhs[1]); // last timestep synaptic current
    //Isynprime
    plhs[2] = mxCreateDoubleMatrix(N*Nsyn, 1, mxREAL);
    Isynprime_out=mxGetPr(plhs[2]); // last timestep synaptic current


    plhs[3] = mxCreateDoubleMatrix(N*Nsyn, 1, mxREAL);
    Isyn_in=mxGetPr(plhs[3]); // last timestep synaptic current
    //Isynprime
    plhs[4] = mxCreateDoubleMatrix(N*Nsyn, 1, mxREAL);
    Isynprime_in=mxGetPr(plhs[4]); // last timestep synaptic current


    //Vm for all timesteps
    plhs[5] = mxCreateDoubleMatrix(Nrecord, Nt, mxREAL);
    vr=mxGetPr(plhs[5]); // record memebrane potential

    // last timestep Vm recorded
    plhs[6] = mxCreateDoubleMatrix(N,1,mxREAL);
    Vout = mxGetPr(plhs[6]);

    //  plhs[3] = mxCreateDoubleMatrix(4,T/modnum,mxREAL);
    //  iappScale=mxGetPr(plhs[3]);
    // ramped input values
    plhs[7] = mxCreateDoubleMatrix(4,Nt,mxREAL);
    IappScaleOutput = mxGetPr(plhs[7]);
    //iappScale=mxGetPr(plhs[3]);

    plhs[8] = mxCreateDoubleMatrix(N,1,mxREAL);
    refstate_out = mxGetPr(plhs[8]);

	plhs[9] = mxCreateDoubleMatrix(Nrecord*Nsyn, Nt, mxREAL);
    Isynrecord=mxGetPr(plhs[9]); // record synaptic current

	plhs[10] = mxCreateDoubleMatrix(Nrecord*Nsyn, Nt, mxREAL);
    Isynprimerecord=mxGetPr(plhs[10]); // record synaptic current

    



    /* Allocate membrane potential */
    v = mxMalloc(N*sizeof(double));; // membrane potential
    refstate=mxMalloc(N*sizeof(int)); // refractory state
    Ntref=mxMalloc(Npop*sizeof(int)); //(????)

   	Isyn = mxMalloc(Nsyn*N*sizeof(double));  /* synp. currents, NxNsyn, Nsyn: nnz of Psyn, col i corresponds to syntype[i]   */
   	Isynprime=mxMalloc(Nsyn*N*sizeof(double));


    Nsum=mxMalloc((Npop)*sizeof(int)); /* number of neurons before pop i */
    Nsum[0]=(int)Ncell[0];//set initial number of neurons to initial number of cells per population in reccurent layer
    mexPrintf("i=%d, Nsum=%d \n",0,(int)Nsum[0]);
    for(i=1;i<Npop; i++){
    //    mexPrintf("i=%d, Ncell=%d \n",i,(int)Ncell[i]);
        Nsum[i] = Nsum[i-1] + (int)Ncell[i];
     //   mexPrintf("i=%d, Nsum=%d \n",i,(int)Nsum[i]);
    } // Nsum results as vector: [Ncell[0], Ncell[0]+Ncell[1], Ncell[0]+Ncell[1]+Ncell[2],...]
    //  mexPrintf(N);
    /* Check for consistency with total number of neurons */
    if(N!=Nsum[Npop-1])
        mexErrMsgTxt("N should equal to the sum of Ncell");


    Ksum=mxMalloc((Npop)*sizeof(int)); /* total number of postsynaptic cells for pop i */
    Kmax = (int)Kr[0]; // set max out-deprees to initial outdegree rec.
    for(j=0;j<Npop; j++){
        Ksum[j]=0;
        for(i=0;i<Npop;i++){
            Ksum[j] = Ksum[j] + (int)Kr[i+j*Npop];
            if(Kmax<(int)Kr[i+j*Npop])
                Kmax = (int)Kr[i+j*Npop];
        }
        /*   mexPrintf("j=%d, Ksum=%d \n",j,Ksum[j]); */
    }
    Kxsum = 0; // total number of projects/out degrees from feedforward networks
    for(i=0;i<Npop; i++){
        Kxsum += (int)Kx[i];
        if(Kmax<(int)Kx[i])
            Kmax = (int)Kx[i];
    }

    /*  mexPrintf("Kxsum=%d, Kmax=%d \n",Kxsum,Kmax);  */

    postcell=mxMalloc(Kmax*sizeof(int)); /* index for postsynaptic cells */

    /*****
    * Finished allocating variables
    ****/

   
    /* Inititalize variables */
    for(j=0;j<N;j++){
        v[j]=v0[j];
        refstate[j]=refstate_init[j];
    }

    /* original loop
    for(jj=0;jj<Nsyn*N;jj++){
        Isyn[jj]=Isyn_init[jj];
        Isynprime[jj]=Isynprime_init[jj];
    }
    */

    for(jj=0;jj<N*Nsyn;jj++){
        Isyn[jj] = Isyn_init[jj];
        Isynprime[jj] = Isynprime_init[jj];



        Isyn_in[jj] = Isyn[jj];
        Isynprime_in[jj] = Isynprime[jj];
    }
    //for(jj=0;jj<N;jj++){
	//	for(isyn = 0; isyn<Nsyn;isyn++){
	//		Isyn[isyn*N+jj] = Isyn_init[jj*Nsyn+isyn];
	//		Isynprime[isyn*N+jj] = Isynprime_init[jj*Nsyn+isyn];


	//	}
	//}

    /* Record first time bin */
    for(jj=0;jj<Nrecord;jj++){
        if(Irecord[jj]<1 || Irecord[jj]>N+1)
            mexErrMsgTxt("Indices in Irecord must be between 1 and N");
        for(isyn=0;isyn<Nsyn;isyn++){
            Isynrecord[isyn*Nrecord+jj]=Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
            Isynprimerecord[isyn*Nrecord+jj]=Isynprime[(int)round(Irecord[jj]-1)*Nsyn+isyn]; // 

        }
        vr[jj]=v[(int)round(Irecord[jj]-1)]; 
    }

    /* Refractory states */
    for(i=0;i<Npop;i++){
        Ntref[i]=(int)round(tref[i]/dt);
        /*  mexPrintf("i=%d, Ntref=%d \n",i,Ntref[i]);*/
    }

    /* Initialize number of spikes */
    ns=0;

    /* Time loop */
    /* Exit loop and issue a warning if max number of spikes is exceeded */
    iXspike=0; // number of spikes recorded
    iXspkInd=&sx[0]; /* even index: spk times. odd index: neuron ID */

    iRspike=0; // number of spikes recorded
    iRspkInd=&srr[0]; /* even index: spk times. odd index: neuron ID */

    double IappScale[4] = {Iapp[0],Iapp[1],Iapp[2],Iapp[3]};
    IappScaleOutput[0] = 1.0*IappScale[0];
    IappScaleOutput[1] = 1.0*IappScale[1];
    IappScaleOutput[2] = 1.0*IappScale[2];
    IappScaleOutput[3] = 1.0*IappScale[3];
    for(i=1;i<=Nt && ns<maxns;i++){
        IappScaleOutput[i*4] = 1.0*IappScale[0];
        IappScaleOutput[i*4+1] = 1.0*IappScale[1];
        IappScaleOutput[i*4+2] = 1.0*IappScale[2];
        IappScaleOutput[i*4+3] = 1.0*IappScale[3];


        /* Update synaptic variables */
        for(jj=0;jj<N*Nsyn;jj++){
            isyn=jj%Nsyn; // type of synapse
            Isyn[jj]+=Isynprime[jj]*dt; //step forward in time
            Isynprime[jj]-=dt*(Isynprime[jj]*temp1[isyn]+Isyn[jj]*temp2[isyn]); // calculation so Isynprime in next timestep
        }

        /* Find all spikes in feedforward layer at this time bin */
        /* Add to corresponding elements of JnextX */
        while(*iXspkInd<=i*dt && iXspike<Nsx){
            iXspkInd++; /* point to neuron ID */
            // iXspkInd is a pointer. *iXspkIDn is the value this pointer is pointing to. iXspkInd is an address
            // Before this loop, iXspkInd points to the first element in s matrix, which is the first spike time.
            // When starting this loop, iXspkInd++ makes iXskpInd points to the first element in the second row,
            // which is the index of neuron of the first spike. The last line of this loop moves the points to the
            // next column, but first row. So when the next loop starts, this line makes iXspkInd points elements
            // in the second row. Basically, *iXspkInd is spike time.
            jspike=(int)round((*iXspkInd)-1); //neuron ID
            if(jspike<0 || jspike>=Nx){
                mexPrintf("\n %d %d %d %d %d\n",(int)round(*(iXspkInd-1)/dt),iXspike,i,jspike,(int)round((*iXspkInd)-1));
                mexErrMsgTxt("Out of bounds index in sx.");
            }
            Pr=&Wrf[jspike*Kxsum]; // All X neurons have Kxsum projections split over Npops based on prob and Kx
            for(post_pop=0;post_pop<Npop;post_pop++){
                if(Kx[post_pop]>0){ //if projection of feedfoward at post_pop >0, jump in loop
                    for(k=0;k<Kx[post_pop];k++){ // step through specific Kx population projection
                        postcell[k]=(int)((*Pr)-1); // kth postsynaptic neuron in this specific presynaptic neuron points to, assigns postsyn the presyn current
                        if(postcell[k]<0 || postcell[k]>=N){
                            mexPrintf("\n Wrf j=%d, postcell=%d\n",jspike, postcell[k]);
                            mexErrMsgTxt("postcell out of bounds");}
                        Pr++;}
                }
                for(isyn=0;isyn<Nsyn;isyn++){
                    if(syntype[isyn]==0){
                        for(k=0;k<Kx[post_pop];k++)
                            Isynprime[postcell[k]*Nsyn+isyn]+=Jx[post_pop]*temp2[isyn]; // update I current history
                    }
                }
            }
            iXspike++;
            iXspkInd++;
        }
        //mexPrintf("here lies the issue");
        //mexPrintf("%d",Nsrr);

        if (Nsrr > 0) {
            while(*iRspkInd<=i*dt && iRspike<Nsrr){
                iRspkInd++; /* point to neuron ID */
                rspike=(int)round((*iRspkInd)-1); //neuron ID
                if(rspike<0 || rspike>=N){ // chek within index of srr
                    mexPrintf("\n %d %d %d %d %d\n",(int)round(*(iRspkInd-1)/dt),iRspike,i,rspike,(int)round((*iRspkInd)-1));
                    mexErrMsgTxt("Out of bounds index in srr.");
                }
                for(isyn=0;isyn<Nsyn;isyn++){ // update I current history wtih temp 2, Jrr established in Wrr loop
                    if(syntype[isyn]==0){
                        Isynprime[rspike*Nsyn+isyn]+=Jext[0]*temp2[isyn]; 
                    }
                }
                iRspike++;
                iRspkInd++;
            }
        }
        /* loop over neurons  */
        Pr=&Wrr[0];
        pre_pop = 0;
        for(j=0;j<N;j++){
            
            /* Update membrane potential */
            /* Spikes will be propagated at the END of the time bin (see below)*/
            if(j>=Nsum[pre_pop]){
                pre_pop = pre_pop+1;
                if ((fmod(i*dt,modnum) == 0) && (pre_pop+1 == popwithinput)){
                    IappScale[pre_pop] = IappScale[pre_pop]+coeffstep;
                }

            }

          //  if ((fmod(i*dt,modnum) == 0) && (pre_pop+1 == popwithinput)){
          //      mexPrintf("%f\t,%f, IappScale %f %f %f %f\n",i*dt,fmod(i*dt,modnum),IappScale[0],IappScale[1],IappScale[2],IappScale[3]);
          //  }

            if(refstate[j]<=0){
                Isyntot=0;
                for(isyn=0;isyn<Nsyn;isyn++){
                    Isyntot+=Psyn2[isyn]*Isyn[j*Nsyn+isyn]; // update Isyntot
                } //Iapp[pre_pop * Nt + i]
                v[j]+=fmax((IappScale[pre_pop]+Isyntot-gl[pre_pop]*(v[j]-Vleak[pre_pop])+gl[pre_pop]*DeltaT[pre_pop]*EXP((v[j]-VT[pre_pop])/DeltaT[pre_pop]))*dt/C[pre_pop],Vlb[pre_pop]-v[j]);} //EIF voltage
            else{
                if(refstate[j]>1)
                    v[j]=Vth[pre_pop]; // v[j] is set to voltage threshold of pre-neuron
                else
                    v[j]=Vre[pre_pop]; //v[j] is set to voltage reset potential of pre-neuron
                refstate[j]--;
            }

            /* If a spike occurs */
            if(v[j]>=Vth[pre_pop] && refstate[j]<=0 && ns<maxns){
                // reset neuron such that it is ready to fire again
                refstate[j]=Ntref[pre_pop];

                v[j]=Vth[pre_pop];       /* reset membrane potential */
                s[0+2*ns]=i*dt; /* spike time */
                s[1+2*ns]=j+1;     /* neuron index 1 */
                ns++;           /* update total number of spikes */
                /* For each postsynaptic target, update synaptic inputs */
                for(post_pop=0;post_pop<Npop;post_pop++){
                    if(Kr[post_pop+pre_pop*Npop]>0){ //if recurrent out-degrees>0, jump in loop
                        for(k=0;k<Kr[post_pop+pre_pop*Npop];k++){ // given that k is less that the out-degrees - enter loop
                            postcell[k]=(int)(*Pr)-1; //post cell becomes presyn weight
                            if(postcell[k]<0 || postcell[k]>=N){ // ensure within bounds
                                mexPrintf("\n neuron j=%d, post_pop=%d, k=%d, postcell=%d\n",j,post_pop,k,postcell[k]);
                                mexErrMsgTxt("postcell out of bounds");
                            }
                            Pr++;
                        }
                        for(isyn=0;isyn<Nsyn;isyn++){
                            if(syntype[isyn]==(pre_pop+1)){ //if synaptic type is equal to firing neuron - enter loop
                                for(k=0;k<Kr[post_pop+pre_pop*Npop];k++)
                                    Isynprime[postcell[k]*Nsyn+isyn]+=Jr[post_pop+pre_pop*Npop]*temp2[isyn];}
                        }
                    }
                }
            }
            else{
                Pr=Pr+Ksum[pre_pop];}
        }
       
     /* Store recorded variables */
        for(jj=0;jj<Nrecord;jj++){
            if(Irecord[jj]<1 || Irecord[jj]>N+1)
                mexErrMsgTxt("Indices in Irecord must be between 1 and N");
            for(isyn=0;isyn<Nsyn;isyn++){
                Isynrecord[isyn*Nrecord+jj+i*Nrecord*Nsyn]=Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
                Isynprimerecord[isyn*Nrecord+jj+i*Nrecord*Nsyn]=Isynprime[(int)round(Irecord[jj]-1)*Nsyn+isyn];
                //Isynprime[postcell[k]*Nsyn+isyn]+=Jr[post_pop+pre_pop*Npop]*temp2[isyn];}
/*
                if (isyn == Nsyn-1){
                    ("\n neuron j=%d: refstate %d",j,(int)refstate[jj]);
                    ("\n neuron j=%d: Isynrecord %d",j,(int)Isynrecord[isyn*Nrecord+jj+i*Nrecord*Nsyn]);
                    ("\n neuron j=%d: Isynprimerecord %d",j,(int)Isynprimerecord[isyn*Nrecord+jj+i*Nrecord*Nsyn]);
                }
                */
              
            }
            vr[jj+Nrecord*i]=v[(int)round(Irecord[jj]-1)];
            //mexPrintf("Irecord[jj] %d , int(Irecord[jj]-1) %d \n", Irecord[jj], (int)round(Irecord[jj]-1));
        }
			//mexPrintf("i = %d\n",(int)i);

    //}  //// time loops ends


		if(i == (int)(Nt)){
			mexPrintf("saving loop \n",(int)i);
			mexPrintf("i = %d\n",(int)i);
			mexPrintf("Nt = %d\n",(int)Nt);

           /* for(jj=0;jj<Nrecord;jj++){
            //if(Irecord[jj]<1 || Irecord[jj]>N+1)
            //    mexErrMsgTxt("Indices in Irecord must be between 1 and N");
                for(isyn=0;isyn<Nsyn;isyn++){
                 //   mexPrintf("last time step record loop");
                    Isyn_last[isyn*Nrecord+jj]= Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
                    Isynprime_last[isyn*Nrecord+jj]=Isynprime[(int)round(Irecord[jj]-1)*Nsyn+isyn];
                }
            }*/

            //Isyn_initial[isyn*Nrecord+jj]=Isyn[(int)round(Irecord[jj]-1)*Nsyn+isyn];
            //Isynprime_initial[isyn*Nrecord+jj]=Isynprime[(int)round(Irecord[jj]-1)*Nsyn+isyn]; // come back


			//mexPrintf("i = %d\n",(int)i);
			for(jj=0;jj<N;jj++){
			//	mexPrintf("saving vout and refout \n");
				Vout[jj] = 1.0*v[jj];
				refstate_out[jj] = 1.0*refstate[jj];
			}
            for(jj=0;jj<N*Nsyn;jj++){
                    Isyn_out[jj] = Isyn[jj];
                    Isynprime_out[jj] = Isynprime[jj];

            }

			/*
			for(jj=0;jj<Nsyn*N;jj++){
        		Isyn_out[jj]=Isyn[jj];
        		Isynprime_out[jj]=Isynprime[jj];
    		}*/

			/*for(jj=0; jj<N; jj++){
				for(isyn=0; isyn<Nsyn; isyn++){
					Isyn_out[isyn*N+jj+N*Nsyn] = Isyn[jj*Nsyn+isyn];
					Isynprime_out[isyn*N+jj+N*Nsyn] = Isynprime[jj*Nsyn+isyn];
				}
			}*/

			//for(jj=0;jj<Nsyn*N;jj++){
        	//	Isyn_out[jj]=Isyn[jj];
        	//	Isynprime_out[jj]=Isynprime[jj];
    		//}
		/*    		        for(jj=0;jj<Nrecord;jj++){

    		for(isyn=0;isyn<Nsyn;isyn++){
                for(k=0;k<Kx[post_pop];k++)
                    Isynprime[postcell[k]*Nsyn+isyn]+=Jx[post_pop]*temp2[isyn]; // update I current history
                    }
               }
               */
		

				//mexPrintf("jj = %d\n",(int)jj);

		    // if(i==Nt-1){
                  //  if(jj == 1 || jj == N){
                  //      mexPrintf("final dt to record out i =%d, jj = %d",i,jj);       
           //Vout[jj] = 1.0*v[jj*Nt2];
           // mexPrintf("\n neuron jj=%d",jj);
       //     refstate_out[jj] = 1.0*refstate[jj];
           // mexPrintf("\n neuron jj=%d: refstate %d",jj,(int)round((refstate_out[jj]));
       //     Isynrecord_out[jj] = 1.0*Isynrecord[jj*Nt];
           // mexPrintf("\n neuron jj=%d: Isynrecord %d",jj,(int)round(Isynrecord[jj*Nt]);
      //      Isynprimerecord_out[jj] = 1.0*Isynprimerecord[jj*Nt];
          //  mexPrintf("\n neuron jj=%d: Isynprimerecord %d",jj,(int)round(Isynprimerecord[jj*Nt]));
       // }
	    
	//timeloop now ends here

			/*
			for(jj=0;jj<Nsyn*N;jj++){
			    Isyn_last[jj]=Isyn[jj];
			    Isynprime_last[jj]=Isynprime[jj];
			}
			*/
		}

	}


    /* Issue a warning if max number of spikes reached */
    if(ns>=maxns)
        mexWarnMsgTxt("Maximum number of spikes reached, simulation terminated.");

    /* Free allocated memory */
    mxFree(v);
    mxFree(refstate);
    mxFree(Isyn);
    mxFree(Isynprime);
    mxFree(temp1);
    mxFree(temp2);
    mxFree(Ksum);
    mxFree(Nsum);
    mxFree(postcell);
    mxFree(syntype);
    mxFree(Ntref);
}
