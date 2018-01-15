package machinelearning;

import java.io.BufferedReader;   
import java.io.FileReader;
import java.io.FileWriter;

// This Java code for RandomForest was written by Dr. Renxiang Yan in Fuzhou University.
// Please email bugs/problems to him (yanrenxiang@fzu.edu.cn).
// This code was heavily affected by the Fortran code of Breiman in
// https://www.stat.berkeley.edu/~breiman/RandomForests/cc_software.htm
// Compile the code by using: gfortran.exe -fno-range-check RRR1.f -o R
// 2017-7-15v0.6
// 119: /home/renxyan/RRR/new
// RRR1.f
// tRRR2.f
// final error test %      0.
// 888 Test set confusion matrix:
//          true class 
//          1     2
//    1     2     0
//    2     0     3

public class RandomForest7{
	static int mdim,ntrain,nclass,maxcat,ntest;
    static int labelts,labeltr,mtry0,ndsize,jbt,look,lookcls;
    static int jclasswt,mdim2nd,mselect,imp,interact,impn;
    static int nprox,nrnn,noutlier,nscale,nprot,missfill,iviz;
    static int isaverf,isavepar,isavefill,isaveprox;
    static int  irunrf,ireadpar,ireadfill,ireadprox;
    static float code;
    
    //OUTPUT CONTROLS
    static int isumout,idataout,impfastout,impout,impnout,interout;
    static int iprotout,iproxout,iscaleout,ioutlierout;
	
    //DERIVED PARAMETERS (DO NOT CHANGE)
    static int nsample,nrnodes,mimp,near;
    static int ifprot,ifscale,iftest,mdim0,ntest0,nprot0,nscale0;

    //  DIMENSIONING OF ARRAYS (End variables)
    static int nmfmax=-1;
    static int ncsplit,ncmax,nmissfill,ndimreps,nmf,nmd,iseed;
	static String text = "";
    
	// SCALAR DECLARATIONS
	static float errtr,errts,tavg,er,randomu;
	static int mtry,n,m,mdimt,m1,jb,nuse,ndbigtree,jj,mr;
    static int  n0,n1,n2,n3,n4,n5,n6,n7;
         
	static float[][] x;
	static float[][] xts;
	static float[] v5;
	static float[] v95;
	static float[] tgini;
	static float[] zt;
	static float[] avgini;	    
	static float[][] votes;
	static float[][] effect;
	static float[][] teffect;	    
	static float[][] hist;
	static float[] g;
	static float[] fill;
	static float[][] rinpop;	    
	static float[] dgini;
	static float[] xbestsplit;
	static float[] tnodewt;	    
	static float[] tw;
	static float[] tn;
	static float[] v;
	static float[] win;
	static float[] temp;	    
	static float[][] q;
	static float[] devout;
	static float[] classwt;
	static float[] wr;
	static float[] tmissts;
	static float[] tmiss;
	static float[] tclasspop;
	static float[] wl;
	static float[] rmedout;
	static float[][] tclasscat;
	static float[][] qts;	    
	static float[][] classpop;
	static float[] signif;
	static float[] zscore;
	static float[] sqsd;
	static float[] avimp;
	static float[] qimp;
	static float[][] qimpm;
	static float[] tout;
	static float[] outtr;
	static float[] xc;
	static float[] dn;
	static float[] cp;
	static float[] cm;	    
	static float[] votecat;
	static float[] freq;
	static float[] wc;
	static float[] outts;
	static float[][] popclass;
	static float[][][] protlow;
	static float[][][] prot;
	static float[][][] prothigh;
	static float[][][][] protfreq;
	static float[] rpop;
	static float[][][] protv;
	static float[] wtx;
	static float[][][] protvlow;
	static float[][][] protvhigh;	    
	    
	static int[] cat;		
	static int[] iv;
	static int[] msm;		
	static int[] muse;
	static int[][] irnk;
	static int[][] missing;
	static int[][] a;	    
	static int[][] asave;
	static int[][] b;	    
	static int[] cl;
	static int[] out;
	static int[] nodextr;
	static int[] nodexvr;
	static int[] jin;
	static int[] joob;
	static int[] pjoob;
	static int[][] ndbegin;	    
	static int[] jvr;
	static int[] jtr;
	static int[] jest;
	static int[] ibest;	    
	static int[] isort;
	static int[][] loz;	    
	static int[] ta;
	static int[] ncase;
	static int[] idmove;
	static int[] kpop;	    
	static int[] jests;
	static int[] jts;
	static int[] iwork;	    
	static int[] nodexts;
	static int[] clts;
	static int[] imax;
	static int[][] jinb;
	static int[] bestsplitnext;
	static int[] bestvar;
	static int[] bestsplit;
	static int[] nodestatus;
	static int[] nodepop;
	static int[] nodestart;
	static int[] nodeclass;
	static int[] parent;
	static int[][] treemap;
	static int[] ncts;
	static int[] nc;
	static int[][] mtab;
	static int[] ncn;
	static int[] its;
	static int[] jpur;
	static int[] npend;
	static int[] inear;
	static int[] nrcat;
	static int[] kcat;
	static int[] ncatsplit;
	static int[][] nbestcat;
	static int[] ncp;
	static int[][] nodexb;
	static int[][] npcase;
	static int[][] ncount;
	static int[] nod;	      
	
    // USED IN PROXIMITY AND SCALING	     
    static float[][] prox;
    static float[] y;
    static float[] u;
    static float[] dl;
    static float[][] xsc;
    static float[] red;
    static float[] ee;
    static float[][] ev;
    static float[] ppr;
    
    static int jstat=-999999;   
    static int ndendl=-999,ndstart=-999,ndend=-999;
    static float decsplit = 1;   
    static int msplit=0,nbest=0;
    
	public static void main(String[] args) throws Exception {	
		// my variable
		mdim=5;
	    ntrain=5;
	    nclass=2;
	    maxcat=1;
	    ntest=5;
	    labelts=1;
	    labeltr=1;
	    mtry0=3;
	    ndsize=1;
	    jbt=3;
	    look=10;
	    lookcls=1;
	    jclasswt=0;
	    mdim2nd=0;
	    mselect=0;
	    imp=0;
	    interact=0;
	    impn=0;
	    nprox=0;
	    nrnn=ntrain;
	    noutlier=0;
	    nscale=0;
	    nprot=0;
	    code=-999;
	    missfill=0;
	    iviz=0;
	    isaverf=1;
	    isavepar=1;
	    isavefill=1;
	    isaveprox=1;
	    irunrf=1;
	    ireadpar=0;
	    ireadfill=0;
	    ireadprox=0;	    
	    // output control
	    isumout = 1;  //!0/1          1=summary to screen
	    idataout= 2;  //!0/1/2        1=train,2=adds test(7)
	    impfastout=0; //!0/1          1=gini fastimp  (8)
	    impout    =1; //!0/1/2        1=imp,2=to screen(9)
	    impnout   =1; //!0/1          1=impn          (10)
	    interout  =0; //!0/1/2        1=interaction,2=screen(11)
	    iprotout  =1; //!0/1/2        1=prototypes,2=screen(12)
	    iproxout  =1; //!0/1/2        1=prox,2=adds test(13)
	    iscaleout =0; //!0/1          1=scaling coors (14)
	    ioutlierout=1;//) !0/1/2       1=train,2=adds test (15)
	    
	    // DERIVED PARAMETERS (DO NOT CHANGE)
	    nsample=(2-labeltr)*ntrain;
	    nrnodes=2*nsample+1;
	    mimp=imp*(mdim-1)+1;
	    ifprot=(int)(nprot/(nprot-.1));
	    ifscale=(int)(nscale/(nscale-.1));
	    iftest =(int)(ntest/(ntest-.1));
	    nprot0=(1-ifprot)+nprot;
	    nscale0=(1-ifscale)+nscale;
	    ntest0=(1-iftest)+ntest;
	    mdim0=interact*(mdim-1)+1;
	    near=nprox*(nsample-1)+1;	    
	    System.out.println("ntest0=" + ntest0);	    
	    // float good start	    
	    x = new float[mdim+1][nsample+1]; 
	    xts = new float[mdim+1][ntest0+1]; 
	    v5 = new float[mdim+1]; 
	    v95 = new float[mdim+1]; 
	    tgini = new float[mdim+1]; 
	    zt = new float[mdim+1]; 
	    avgini = new float[mdim+1]; 
	    votes = new float[mdim0+1][jbt+1]; 
	    effect = new float[mdim0+1][mdim0+1]; 
	    teffect=new float[mdim0+1][mdim0+1]; 
	    hist = new float[mdim0+1][mdim0+1]; 
	    g = new float[mdim0+1]; 
	    fill = new float[mdim+1]; 
	    rinpop = new float[near+1][jbt+1]; 
	    dgini = new float[nrnodes+1]; 
	    xbestsplit = new float[nrnodes+1]; 
	    tnodewt = new float[nrnodes+1]; 
	    tw = new float[nrnodes+1]; 
	    tn = new float[nrnodes+1]; 
	    v = new float[nsample+1]; 
	    win = new float[nsample+1]; 
	    temp = new float[nrnn+1]; 
	    q = new float[nclass+1][nsample+1]; 
	    devout = new float[nclass+1]; 
	    classwt = new float[nclass+1]; 
	    wr = new float[nclass+1]; 
	    tmissts = new float[nclass+1]; 
	    tmiss = new float[nclass+1]; 
	    tclasspop = new float[nclass+1]; 
	    wl = new float[nclass+1]; 
	    rmedout = new float[nclass+1]; 
	    tclasscat = new float[nclass+1][maxcat+1]; 
	    qts = new float[nclass+1][ntest0+1]; 
	    classpop = new float[nclass+1][nrnodes+1]; 
	    signif = new float[mimp+1]; 
	    zscore = new float[mimp+1]; 
	    sqsd = new float[mimp+1]; 
	    avimp = new float[mimp+1]; 
	    qimp = new float[nsample+1]; 
	    qimpm = new float[nsample+1][mimp+1]; 
	    tout = new float[near+1]; 
	    outtr = new float[near+1]; 
	    xc = new float[maxcat+1]; 
	    dn = new float[maxcat+1]; 
	    cp = new float[maxcat+1]; 
	    cm = new float[maxcat+1]; 
	    votecat = new float[maxcat+1]; 
	    freq = new float[maxcat+1]; 
	    wc = new float[nsample+1]; 
	    outts = new float[ntest0+1]; 
	    popclass = new float[nprot0+1][nclass+1]; 
	    protlow = new float[mdim+1][nprot0+1][nclass+1]; 
	    prot = new float[mdim+1][nprot0+1][nclass+1]; 
	    prothigh = new float[mdim+1][nprot0+1][nclass+1]; 
	    protfreq = new float[mdim+1][nprot0+1][nclass+1][maxcat+1]; 
	    rpop = new float[nrnodes+1]; 
	    protv = new float[mdim+1][nprot0+1][nclass+1]; 
	    wtx = new float[nsample+1]; 
	    protvlow = new float[mdim+1][nprot0+1][nclass+1]; 
	    protvhigh = new float[mdim+1][nprot0+1][nclass+1];   
	    
	    cat = new int[mdim+1]; 
	    iv = new int[mdim+1]; 
	    msm = new int[mdim+1]; 
	    muse = new int[mdim+1]; 
	    irnk = new int[mdim+1][jbt+1]; 
	    missing = new int[mdim+1][near+1]; 
	    a = new int[mdim+1][nsample+1]; 
	    asave = new int[mdim+1][nsample+1]; 
	    b = new int[mdim+1][nsample+1]; 
	    cl = new int[nsample+1]; 
	    out = new int[nsample+1]; 
	    nodextr = new int[nsample+1]; 
	    nodexvr = new int[nsample+1]; 
	    jin = new int[nsample+1]; 
	    joob = new int[nsample+1]; 
	    pjoob = new int[nsample+1]; 
	    ndbegin = new int[near+1][jbt+1]; 
	    jvr = new int[nsample+1]; 
	    jtr = new int[nsample+1]; 
	    jest= new int[nsample+1]; 
	    ibest = new int[nrnn+1]; 
	    isort = new int[nsample+1]; 
	    loz = new int[near+1][nrnn+1]; 
	    ta = new int[nsample+1]; 
	    ncase = new int[nsample+1]; 
	    idmove = new int[nsample+1]; 
	    kpop = new int[nrnodes+1]; 
	    jests = new int[ntest0+1]; 
	    jts = new int[ntest0+1]; 
	    iwork = new int[near+1]; 
	    nodexts = new int[ntest0+1]; 
	    clts = new int[ntest0+1]; 
	    imax = new int[ntest0+1]; 
	    jinb = new int[near+1][jbt+1]; 
	    bestsplitnext = new int[nrnodes+1]; 
	    bestvar = new int[nrnodes+1]; 
	    bestsplit = new int[nrnodes+1]; 
	    nodestatus = new int[nrnodes+1]; 
	    nodepop = new int[nrnodes+1]; 
	    nodestart = new int[nrnodes+1]; 
	    nodeclass = new int[nrnodes+1]; 
	    parent = new int[nrnodes+1]; 
	    treemap = new int[2+1][nrnodes+1]; 
	    ncts = new int[nclass+1]; 
	    nc = new int[nclass+1]; 
	    mtab = new int[nclass+1][nclass+1]; 
	    ncn = new int[near+1]; 
	    its = new int[nsample+1]; 
	    jpur = new int[nrnn+1]; 
	    npend = new int[nclass+1]; 
	    inear = new int[nrnn+1]; 
	    nrcat = new int[maxcat+1]; 
	    kcat = new int[maxcat+1]; 
	    ncatsplit = new int[maxcat+1]; 
	    nbestcat = new int[maxcat+1][nrnodes+1]; 
	    ncp = new int[near+1]; 
	    nodexb = new int[near+1][jbt+1]; 
	    npcase = new int[near+1][jbt+1]; 
	    ncount = new int[near+1][jbt+1]; 
	    nod = new int[nrnodes+1];     
		    
	    // USED IN PROXIMITY AND SCALING	     
	    prox = new float[near+1][nrnn+1]; 
	    y = new float[near+1]; 
	    u = new float[near+1]; 
	    dl= new float[nscale0+1]; 
	    xsc = new float[near+1][nscale0+1]; 
	    red = new float[near+1]; 
	    ee = new float[near+1]; 
	    ev = new float[near+1][nscale0+1]; 
	    ppr = new float[near+1];
	    // SCALAR DECLARATIONS    
	    String text = "";		

	    int i,j,k;
		String line = "";		
		boolean test = true;
		boolean train = false;
		if(test==true){		
			ntest=5;
			ntest0=5;
			mdim=5;		
			if(ntest>0){
				BufferedReader br_f17 = new BufferedReader(new FileReader("C:\\bioinformatics_package\\javaRF\\satimage.tes"));	       
		        for(n=1;n<=ntest;n++){
		        	    line = br_f17.readLine();
		        	    line = line.replaceFirst("^\\s+", "");
		        	    String[] wds = line.split("\\s+");
		        	    System.out.println("Size=" + wds.length);
		        	    for(m=1;m<=mdim;m++){
		        	    	xts[m][n] = Float.valueOf(wds[m-1]).floatValue();;
		        	    	// System.out.println("Yes: " +xts[m][n]);
		        	    }
		        	    clts[n] = Integer.valueOf(wds[mdim]).intValue();
		        	    // System.out.println("OK:" + clts[n]);
		                // read(17,*) (xts(m,n),m=1,mdim),clts(n)
		        }
		        br_f17.close();
		        br_f17=null;		       		        
			 }			
		     if(mselect==0){
		         mdimt=mdim;
		         for(k=1;k<=mdim;k++){
		                 msm[k]=k;
		         }
		     }		     
		     // SET CATEGORICAL VALUES
		     for(m=1;m<=mdim;m++){
		         cat[m]=1;
		     }
		     if(jclasswt==0){
		          for(j=1;j<=nclass;j++){
		                  classwt[j]=1;
		          }
		     }
		     // MISC PARAMETERS (SHOULD NOT USUALLY NEED ADJUSTMENT)
		     iseed=4351;
		     nmfmax = 5;     //  !number of iterations (iterative fillin)
		     ncsplit = 25;   //  !number of random splits (big categoricals)
		     ncmax = 25;       //  !big categorical has more than ncmax levels
		     if(irunrf==1&&ntest>0){
		    	 System.out.println("Running random forest\n");
		    	 runforest(mdim,ntest,nclass,maxcat,nrnodes,labelts,jbt,clts,xts,xbestsplit,qts,treemap,nbestcat,
		    		               nodestatus,cat,jts,jests,bestvar,tmissts,ncts,fill,missfill,code,
		    		               tnodewt,outts,idataout,imax,look,lookcls,nodexts,isumout,mtab,
		    		               "C:\\bioinformatics_package\\2017RRR\\savedforest");	 
		    	 for(int m=1;m<=5;m++){
		    		 System.out.println("m=" + m + " pro1=" + qts[1][m] + " pro2=" + qts[2][m]);
		    	 }
		     }			
		}else if(train==false){
			//System.exit(1);
		}			
		System.exit(1);
		// Below is the training procedure
		// open(16, file='satimage.tra', status='old')
		// do n=1,ntrain
		// 	read(16,*) (x(m,n),m=1,mdim),cl(n)
		// enddo
		// close(16)		
		String[] wds;
		int n = 1;
		BufferedReader br = new BufferedReader(new FileReader("C:\\bioinformatics_package\\javaRF\\satimage.tra"));
			while((line=br.readLine())!=null){
				//line = line.replaceFirst("^\\s+", "");
				wds = line.split("\\s+");
				for(m=1;m<=mdim;m++){
					x[m][n] = Float.valueOf(wds[m-1]).floatValue();
				}
				cl[n] = Integer.valueOf(wds[mdim]).intValue();
				n++;
			}
		br.close();
		br = null;		
		n = 1;
		ntest=5;
		if(ntest>0){
			// open(17, file='satimage.tes', status='old')
			// do n=1,ntest0
			//	read(17,*) (xts(m,n),m=1,mdim),clts(n)
			// enddo
			// close(17)
			br = new BufferedReader(new FileReader("C:\\bioinformatics_package\\javaRF\\satimage.tes"));
			while((line=br.readLine())!=null){
				//line = line.replaceFirst("^\\s+", "");
				wds = line.split("\\s+");
				for(m=1;m<=mdim;m++){
					xts[m][n] = Float.valueOf(wds[m-1]).floatValue();
				}
				clts[n] = Integer.valueOf(wds[mdim]).intValue();
				n++;
			}
			br.close();
			br = null;
		}		
		
		// INITIALIZATION
		// mtry will change if mdim2nd>0, so make it a variable
		mtry=mtry0;
		// nmissfill is the number of missing-value-fillin loops
		nmissfill=1;
		nmfmax=-1000;
		if(missfill==2) nmissfill=nmfmax;
		// ndimreps is 2 if we want to do variable selection,
		// otherwise it's 1
		ndimreps=1;
		if(mdim2nd>0) ndimreps=2;
		// assign weights for equal weighting case 
		// (use getweights for unequal weighting case)
		if(jclasswt==0){
			for(k=1;k<=nrnodes;k++){
				tnodewt[k]=1;
			}
		}
		if(labeltr==1){
			for(n=1;n<=nsample;n++){
				wtx[n]=classwt[cl[n]];
			}
		}else{
			for(n=1;n<=nsample;n++){
				wtx[n]=1.0f;
			}
		}
		// IF DATA ARE UNLABELED, ADD A SYNTHETIC CLASS
		if(labeltr==0){
			createclass(x,cl,ntrain,nsample,mdim);
		}
		// COUNT CLASS POPULATIONS		
	    zerv(nc,nclass);
		for(n=1;n<=nsample;n++){
			nc[cl[n]]=nc[cl[n]]+1;
		}
		if(ntest>0&&labelts==1){
			zerv(ncts,nclass);
			for(n=1;n<=ntest0;n++){
				ncts[clts[n]]=ncts[clts[n]]+1;
			}
		}
		// DO PRELIMINARY MISSING DATA		
		if(missfill==1){
		    roughfix(x,v,ncase,mdim,nsample,cat,code,nrcat,maxcat,fill);
			if(ntest>0){
				xfill(xts,ntest,mdim,fill,code);
			}
		}
		
		if(missfill==2){
			for(n=1;n<=near;n++){
				for(m=1;m<=mdim;m++){
					if(Math.abs(code-x[m][n])<8.232D-11){
						missing[m][n]=1;
					}else{ 
						missing[m][n]=0;
					}
				}
			}
			roughfix(x,v,ncase,mdim,nsample,cat,code,nrcat,maxcat,fill);
		}		
		
		iseed=1009990;
		ncsplit=1000;
		ncmax=-1000;		
		FileWriter rf_fw = new FileWriter("C:\\MHMM\\rf.txt");
		// TOP OF LOOP FOR ITERATIVE MISSING VALUE FILL OR SUBSET SELECTION
		for(nmd=1;nmd<=ndimreps;nmd++){
				if(imp>0){
					zervr(sqsd,mimp);
					zervr(avimp,mimp);
				}
				zervr(avgini,mdim);
				if(impn==1){
					zervr(qimp,nsample);
					zermr(qimpm,nsample,mdim);
				}
				if(interact==1) zermr(votes,mdim,jbt);
				for(nmf=1;nmf<=nmissfill;nmf++){			
					//INITIALIZE FOR RUN			
					makea(x,mdim,nsample,cat,isort,v,asave,b,mdimt,msm,v5,v95,maxcat);			
					zermr(q,nclass,nsample);
					zerv(out,nsample);
					if(nprox>0) zermd(prox,near,nrnn);
					if(ntest>0){
				 		zermr(qts,nclass,ntest);
					}			
					//   BEGIN MAIN LOOP  			
					for(jb=1;jb<=jbt;jb++){					 
						//INITIALIZE				
						zervr(win,nsample);
						zerv(jin,nsample);
						zervr(tclasspop,nclass);		
						//TAKE A BOOTSTRAP SAMPLE 			
						for(n=1;n<=nsample;n++){
							k=(int)(randomu()*nsample);
							if(k<nsample) k=k+1;
							win[k]=win[k]+classwt[cl[k]];
							jin[k]=jin[k]+1;
							tclasspop[cl[k]]=tclasspop[cl[k]]+classwt[cl[k]];
						}
						for(n=1;n<=nsample;n++){
							if(jin[n]==0) out[n]=out[n]+1;
							//out(n)=number of trees for which the nth 
							//obs has been out-of-bag
						}					
						//PREPARE TO BUILD TREE		
						moda(asave,a,nuse,nsample,mdim,cat,maxcat,ncase,jin,mdimt,msm);				
						//MAIN TREE-BUILDING 			
						
						// delete ndbigtree in the Function and set it as a global variable
						buildtree(a,b,cl,cat,mdim,nsample,nclass,treemap,bestvar,bestsplit,bestsplitnext,dgini,
					     	nodestart,classpop,tclasspop,tclasscat,ta,idmove,ndsize,ncase,parent,mtry,
					     	win,wr,wl,nuse,kcat,ncatsplit,xc,dn,cp,cm,nbestcat,msm,mdimt);
						
						System.out.println("building tree..  ndbigtree=" + ndbigtree);
					    // SPLIT X: delete ndbigtree in the Function and set it as a global variable		
						// delte bestvar and set it as a global variable
						xtranslate(x,mdim,nrnodes,nsample,bestsplit,bestsplitnext,xbestsplit,nodestatus,cat);		
						
						// ASSIGN CLASSWEIGHTS TO NODES: delete ndbigtree in the Function and set it as a global variable			
						if(jclasswt==1){
							getweights(x,nsample,mdim,treemap,nodestatus,xbestsplit,bestvar,nrnodes,cat,maxcat,nbestcat,jin,win,tw,tn,tnodewt);
						}									
						
						// GET OUT-OF-BAG ESTIMATES			
						testreebag(x,nsample,mdim,treemap,nodestatus,xbestsplit,nrnodes,kpop,
					     	cat,jtr,nodextr,maxcat,nbestcat,rpop,dgini,tgini,jin,wtx);
						
						for(n=1;n<=nsample;n++){
							if(jin[n]==0){
								//this case is out-of-bag
								q[jtr[n]][n]=q[jtr[n]][n]+tnodewt[nodextr[n]];
							}
						}
						for(k=1;k<=mdimt;k++){
							m=msm[k];
							avgini[m]=avgini[m]+tgini[m];
							if(interact==1) votes[m][jb]=tgini[m];
						}		
						
						// DO PRE-PROX COMPUTATION						
						if(nprox==1){
							preprox(near,nrnodes,jbt,nodestatus,ncount,jb,
					     		nod,nodextr,nodexb,jin,jinb,ncn,ndbegin,kpop,rinpop,npcase,rpop,nsample);
						}			
						
						// GET TEST SET ERROR ESTIMATES
						if(ntest>0){
							System.out.println("testreelite start....");
							testreelite(xts,mdim,treemap,nodestatus,xbestsplit,nrnodes,cat,jts,maxcat,nbestcat,nodexts);
							for(n=1;n<=ntest0;n++){
								qts[jts[n]][n]=qts[jts[n]][n]+tnodewt[nodexts[n]];
							}
							System.out.println("testreelite end....");
						}	
						
						// GIVE RUNNING OUTPUT			
						if(lookcls==1&&jb==1){
							System.out.println("class counts-training data");
							//write(*,*) (nc(j),j=1,nclass)
							for(j=1;j<=nclass;j++){
								System.out.print(nc[j] + " ");
							}
							System.out.println("");
							if(ntest>0&&labelts==1){
								System.out.println("class counts-test data");
								//write(*,*) (ncts(j),j=1,nclass)
								for(j=1;j<=nclass;j++){
									System.out.print(ncts[j] + " ");
								}
								System.out.println("");
							}
						}
						
						// if(jb%look==0||jb==jbt){
						if(jb==jbt){
							comperrtr(q,cl,nsample,nclass,tmiss,nc,jest,out);
							if(look>0){
								if(lookcls==1){
									//write(*,'(i8,100f10.2)')  jb,100*errtr,(100*tmiss(j),j=1,nclass)
									System.out.println(jb +" " +100*errtr);
									for(j=1;j<=nclass;j++){
										System.out.print(100*tmiss[j] + " ");
									}
								}else{
									//write(*,'(i8,2f10.2)') jb,100*errtr
									System.out.println(jb +" " +100*errtr);
								}
							}
							if(ntest>0){
								comperrts(qts,clts,ntest,nclass,tmissts,ncts,jests,labelts);
								if(look>0){
									if(labelts==1){
									   if(lookcls==1){
										   //write(*,'(i8,20f10.2)') jb,100*errts,(100*tmissts(j),j=1,nclass)
										    System.out.println(jb +" " +100*errtr);
											for(j=1;j<=nclass;j++){
												System.out.print(100*tmiss[j] + " ");
											}
									   }else{
										   //write(*,'(i8,2f10.2)') jb,100*errts
										   System.out.println(jb +" " +100*errtr);
									   }
									   //print *
									}
								}
							}
						}
						
						// VARIABLE IMPORTANCE 			
						if(imp==1&&nmf==nmissfill){
						    varimp(x,nsample,mdim,cl,nclass,jin,jtr,impn,interact,msm,mdimt,qimp,qimpm,avimp,sqsd,
						     	treemap,nodestatus,xbestsplit,bestvar,nodeclass,nrnodes,cat,jvr,nodexvr,maxcat,nbestcat,tnodewt,
						     	nodextr,joob,pjoob,iv);
						}		
						
						System.out.println("ABC");
						// SEND SAVETREE DATA TO FILE (IF THIS IS THE FINAL FOREST)			
						if(isaverf==1&&nmd==ndimreps&&nmf==nmissfill){
							System.out.println("Initial Check point here");							
							if(jb==1){
								//write(1,*) (cat(m),m=1,mdim)
								for(m=1;m<=mdim;m++){
									rf_fw.append(cat[m] + " ");
								}
								rf_fw.append("\n");
								if(missfill==1){ // write(1,*) (fill(m),m=1,mdim)
									for(m=1;m<=mdim;m++){
										rf_fw.append(fill[m] + " ");
									}
									rf_fw.append("\n");
								}
							}
							System.out.println("ndbigtree: " + ndbigtree);
							rf_fw.append(ndbigtree + "\n");
					   		for(n=1;n<=ndbigtree;n++){
								//write(1,*) n,nodestatus(n),bestvar(n),treemap(1,n),
					     		//	treemap(2,n),nodeclass(n),xbestsplit(n),tnodewt(n),
					     		//	(nbestcat(k,n),k=1,maxcat)
					   			System.out.println("Check point here");
					   			// Model file variables 2 and 5 may have problems 
					   			//           0          1                    2                  3
					   			rf_fw.append(n + " " + nodestatus[n] + " " + bestvar[n] + " " + treemap[1][n] + " "
					   		    //           4               5                    6                    7
					     			+ treemap[2][n] + " " + nodeclass[n] + " " + xbestsplit[n] + " " + tnodewt[n] + " ");
					   			for(k=1;k<=maxcat;k++){
					   				rf_fw.append(nbestcat[k][n] + " ");
					   			}		
					   			rf_fw.append("\n");
					   		}
						}								
						//  END MAIN LOOP	
					} // jb
					
					//FIND PROXIMITIES, FILL IN MISSING VALUES AND ITERATE			
					if(nmf<nmissfill){
						System.out.println("nrep" + nmf);					
						//compute proximities between each obs in the training 
						//set and its nrnn closest neighbors				
						comprox(prox,nodexb,jinb,ndbegin,npcase,ppr,rinpop,near,jbt,noutlier,outtr,cl,loz,nrnn,wtx,nsample,iwork,ibest);					
						//use proximities to impute missing values
						impute(x,prox,near,mdim,maxcat,votecat,cat,nrnn,loz,missing);
					}
			} //!nmf
			rf_fw.close();
			rf_fw=null;	
			
			//	COMPUTE IMPORTANCES, SELECT SUBSET AND LOOP BACK			
			if(imp==1){
					finishimp(mdim,sqsd,avimp,signif,zscore,jbt,mdimt,msm);
					for(m=1;m<=mdimt;m++){
						//zt(m) is the negative z-score for variable msm(m)
						zt[m]=-zscore[msm[m]];
						muse[m]=msm[m];
					}
					//sort zt from smallest to largest,and sort muse accordingly
					Y_quicksort(zt,muse,1,mdimt,mdim);
					//muse(m) refers to the variable that has the mth-smallest zt
					//select the mdim2nd most important variables and iterate
					if(mdim2nd>0){
						mdimt=mdim2nd;
						for(m=1;m<=mdimt;m++){
							msm[m]=muse[m];
						}
						mtry=nint((float)(Math.sqrt((mdimt))));
					}
				}
		} //!nmd		
		// END OF ITERATIONS - NOW ENDGAME		
		// NORMALIZE VOTES
		for(j=1;j<=nclass;j++){
			for(n=1;n<=nsample;n++){
				if(q[j][n]>0.0&&out[n]>0) q[j][n]=q[j][n]/(float)(out[n]);
			}
			if(ntest>0){
				for(n=1;n<=ntest0;n++){
					qts[j][n]=qts[j][n]/jbt;
				}
			}
		}
		
		// COMPUTE PROXIMITIES AND SEND TO FILE
		if(nprox>=1){
			// compute proximities between each obs in the training set and its nrnn closest neighbors				
			comprox(prox,nodexb,jinb,ndbegin,npcase,ppr,rinpop,near,jbt,noutlier,outtr,cl,loz,nrnn,wtx,nsample,iwork,ibest);			   
			if(iproxout>=1){ 
				FileWriter fw = new FileWriter("save-run-proximities");
				for(n=1;n<=near;n++){
					//write(13,'(i5,500(i5,f10.3))') n,(loz(n,k),prox(n,k),k=1,nrnn)
					fw.append(n + " ");
					for(k=1;k<=nrnn;k++){
						fw.append(loz[n][k] + " " + prox[n][k]);
					}
					fw.append("\n");
				}
				fw.close();
				fw=null;
			}
		}
		
		// COMPUTE SCALING COORDINATES AND SEND TO FILE
		if(nprox==1&&nscale>0){
			myscale(loz,prox,xsc,y,u,near,nscale,red,nrnn,ee,ev,dl);
			if(iscaleout==1){
				FileWriter fw = new FileWriter("save-scale");
				for(n=1;n<=nscale0;n++){
					//write(14,'(i5,f10.3)') n,dl(n)
					fw.append(n + " " + dl[n] + " ");
				}
				for(n=1;n<=near;n++){
					//write(14,'(3i5,15f10.3)') n,cl(n),jest(n),(xsc(n,k),k=1,nscale0)
					fw.append(n+ " " + cl[n] + " " + jest[n] + " ");
					for(k=1;k<=nscale0;k++){
						fw.append(xsc[n][k] + " ");
					}
					fw.append("\n");
				}							
				fw.close();
				fw=null;
			}
		}
		
		// COMPUTE CASEWISE VARIABLE IMPORTANCE (FOR GRAPHICS) AND SEND TO FILE
		if(impn==1){
			FileWriter fw = new FileWriter("save-caseimp-data");
			for(n=1;n<=nsample;n++){
				for(m1=1;m1<=mdimt;m1++){
					mr=msm[m1];
					qimpm[n][mr]=100*(qimp[n]-qimpm[n][mr])/jbt;
				}
				if(impnout==1){   //write(10,'(100f10.3)')	(qimpm(n,msm(m1)),m1=1,mdimt)
					for(m1=1;m1<=mdimt;m1++){
						fw.append(qimpm[n][msm[m1]] + " ");
					}
				}
			}
			fw.close();
			fw=null;
		}
		
		// SEND IMPORTANCES TO FILE OR SCREEN
		//if(imp==1){
		//	for(k=1;k<=min0(mdimt,25);k++){
		//		m=muse[k];
		//		//if(impout.eq.1) write(9,'(i5,10f10.3)')m,100*avimp(m),zscore(m),signif(m)
		//		if(impout==2){  // write(*,'(i5,10f10.3)')m,100*avimp(m),zscore(m),signif(m)
		//			
		//		}
		//	}
		//	//close(9)
		//}
		
		// COMPUTE INTERACTIONS AND SEND TO FILE OR SCREEN 
		if(interact==1){
			FileWriter fw = new FileWriter("save-pairwise-effects");
			compinteract(votes,effect,msm,mdim,mdimt,jbt,g,iv,irnk,hist,teffect);
			if(interout==1){
				//write(11,*)'CODE'
				fw.append("CODE" + "\n");
				for(i=1;i<=mdimt;i++){
					//write(11,*) i,msm(i)
					fw.append(i + " " + msm[i] + "\n");
				}
				for(i=1;i<=mdimt;i++){
					m=msm[i];
					effect[m][m]=0;
					// write(11,'(40i5)') i,(nint(effect(m,msm(j))),j=1,mdimt)
					fw.append(i + " ");
					for(j=1;j<=mdimt;j++){
						//fw.append((String)(nint(effect[m][msm[j]])));
						int aa = nint(effect[m][msm[j]]);
						fw.append(String.valueOf(aa));
					}
					fw.append("\n");
				}
				//close(11)
				fw.close();
				fw=null;
			}			
			if(interout==2){
				System.out.println("CODE");
				for(i=1;i<=mdimt;i++){
					System.out.println(i + " " + msm[i] + " ");
				}
				for(i=1;i<=mdimt;i++){
					m=msm[i];
					effect[m][m]=0;
					//write(*,'(20i5)') i,(nint(effect(m,msm(j))),j=1,mdimt)
					System.out.print(i + " ");
					for(j=1;j<=mdimt;j++){
						System.out.print(nint(effect[m][msm[j]]));
					}
					
				}
			}
		}
		
		// COMPUTE FASTIMP AND SEND TO FILE
		if(impfastout==1){
			tavg=0;
			for(k=1;k<=mdimt;k++){
				m=msm[k];
				tavg=tavg+avgini[m];
			}
			tavg=tavg/mdimt;
			FileWriter fw = new FileWriter("save-impfast");
			for(k=1;k<=mdimt;k++){
				m=msm[k];
				fw.append(m + " " + avgini[m]/tavg + " \n");
			}
			fw.close();
			fw=null;
		}
		
		// COMPUTE PROTOTYPES AND SEND TO FILE
		if(nprox>=1&&nprot>0){
			compprot(loz,nrnn,nsample,mdim,its,cl,wc,nclass,x,mdimt,msm,temp,cat,maxcat,
      		         jpur,inear,nprot,protlow,prothigh,prot,protfreq,protvlow,protvhigh,protv,
     		         popclass,npend,freq,v5,v95);
//		if (iprotout==1){
//			FileWriter fw = new FileWriter("save-protos");
//			write(12,'(a5,50i10)') '     ',(
//	     &	    	(nint(popclass(i,j)),i=1,npend(j)),j=1,nclass)
//			write(12,'(a5,50i10)') '     ',(
//	     &	    	(i,i=1,npend(j)),j=1,nclass)
//			write(12,'(a5,50i10)') '     ',(
//	     &	    	(j,i=1,npend(j)),j=1,nclass)
//			do k=1,mdimt
//				m=msm(k)
//				if(cat(m).eq.1) then
//		 			write(12,'(i5,50f10.3)') k,
//	     &		     	 	((prot(m,i,j),protlow(m,i,j),
//	     & 		 		prothigh(m,i,j),i=1,npend(j)),j=1,nclass)
//				else 
//		 			write(12,'(i5,50f10.3)') k,
//	     &		     	 	(((protfreq(m,i,j,jj),jj=1,cat(m)),
//	     &				i=1,npend(j)),j=1,nclass)
//				endif
//		 	enddo
//		 	fw.close();
//			fw=null;
//		}		
//	 	if(iprotout.eq.2) then
//	 	write(*,'(a5,50i10)') '     ',((nint(popclass(i,j)),
//     &	    	i=1,npend(j)),j=1,nclass)
//	 	do k=1,mdimt
//	 		m=msm(k)
//			if(cat(m).eq.1) then
//				write(*,'(i5,50f10.3)') k,
//     &				((prot(m,i,j),protlow(m,i,j),
//     &				prothigh(m,i,j),i=1,npend(j)),j=1,nclass)
//			else
//	 			write(*,'(i5,50f10.3)') k,
//     &		     	 	(((protfreq(m,i,j,jj),jj=1,cat(m)),
//     &				i=1,npend(j)),j=1,nclass)
//			endif
//		enddo
//	 	endif
	  }
	
	  // COMPUTE OUTLIER MEASURE AND SEND TO FILE
	  if(noutlier>=1){
		locateout(cl,tout,outtr,ncp,isort,devout,near,nsample,nclass,rmedout);
		if(ioutlierout>=1){
			for(n=1;n<=near;n++){
				//write(15,'(2i5,f10.3)') n,cl(n),amax1(outtr(n),0.0)
				
			}
		}
	  }
	
	  // SUMMARY OUTPUT
	  if(isumout==1){
			System.out.println("final error rate %    " + 100*errtr);
			if(ntest>0&&labelts!=0){
				System.out.println("final error test %    " + 100*errts); 
			}			
			System.out.println("Training set confusion matrix (OOB):");
			zerm(mtab,nclass,nclass);
			for(n=1;n<=nsample;n++){
				if(jest[n]>0) mtab[cl[n]][jest[n]]=mtab[cl[n]][jest[n]]+1;
			}
			System.out.print("	    true class \n ");			
			//write(*,'(20i6)')  (i,i=1,nclass)
			for(i=1;i<=nclass;i++){
				System.out.print(" " + i);
			}
			System.out.println("");
			for(j=1;j<=nclass;j++){
				//write(*,'(20i6)')  j,(mtab(i,j),i=1,nclass)
				System.out.print(j + " ");
				for(i=1;i<=nclass;i++){
					System.out.print(mtab[i][j] + " ");
				}
				System.out.println("");
			}
			System.out.println("");
	 		if(ntest>0&&labelts!=0){
	 		   zerm(mtab,nclass,nclass);
	 		   for(n=1;n<=ntest0;n++){
	 			   mtab[clts[n]][jests[n]]=mtab[clts[n]][jests[n]]+1;
	 		   }
	 		   System.out.println("Test set confusion matrix:");
	 		   System.out.print("	    true class \n ");
	 		   for(i=1;i<=nclass;i++){
	 			   System.out.print(" " + i);
	 		   }
	 		   System.out.println("");
	 		   for(j=1;j<=nclass;j++){
	 			   //write(*,'(20i6)')  j,(mtab(i,j),i=1,nclass)
	 			   System.out.print(j + " ");
	 			   for(i=1;i<=nclass;i++){
	 				   System.out.print(mtab[i][j] + " ");
	 			   }
	 			  System.out.println("");
	 		   }
	 		   //print *
	 		   System.out.println("");
	 		}
		}	
		
	    // SEND INFO ON TRAINING AND/OR TEST SET DATA TO FILE
		if(idataout>=1){
			FileWriter fw = new FileWriter("save-data-from-run.txt");
			for(n=1;n<=nsample;n++){
				 // write(7,'(3i5,5000f10.3)') n,cl(n),jest(n),(q(j,n),j=1,nclass),(x(m,n),m=1,mdim)
				 fw.append(n + " " + cl[n] + " " + jest[n] +" ");
				 for(j=1;j<=nclass;j++){
					 fw.append(q[j][n] + " ");
				 }
				 fw.append("\n");
			}
			
			if(idataout==2&&ntest>0){
				if(labelts==1){
					for(n=1;n<=ntest0;n++){
						//write(7,'(3i5,1000f10.3)') n,clts(n),jests(n),(qts(j,n),j=1,nclass),(xts(m,n),m=1,mdim)
						fw.append(n + " " + clts[n] + " " + jests[n] + " ");
						for(j=1;j<=nclass;j++){
							fw.append(qts[j][n] + " ");
						}
						fw.append("\n");
					}
				}else{
					for(n=1;n<=ntest0;n++){
						//write(7,'(2i5,1000f10.3)') n,jests(n),(qts(j,n),j=1,nclass),(xts(m,n),m=1,mdim)
						fw.append(n + " " + jests[n] + " ");
						for(j=1;j<=nclass;j++){
							fw.append(qts[j][n] + " ");
						}
						fw.append("\n");
					}
				}
		   	}		
			fw.close();
			fw=null;
		} // end idataout>=1
	
		// SEND RUN PARAMETERS TO FILE 
		if(isavepar==1){
			FileWriter fw = new FileWriter("savedparams.txt");
			//write(2,*) ntrain,mdim,maxcat,nclass,jbt,jclasswt,missfill,code,nrnodes,100*errtr
			// write(2,*) 'this is a test run to verify that my descriptive output works.'
			fw.append(ntrain+ " " + mdim + " " + maxcat + " " + nclass + " " + jbt + " " + jclasswt + " "+ missfill + 
					" " + code + " " + nrnodes + " " + "error_rate=" + 100*errtr);
			fw.close();
			fw=null;
		}		
	} 
	
	// fun1: runforest, delete errts variable from function, set it as global variable
    public static void runforest(int mdim,int ntest,int nclass,int maxcat, int nrnodes,
     int labelts,int jbt,int[] clts,float[][] xts,float[] xbestsplit,float[][] qts,int[][] treemap,int[][] nbestcat,
     int[] nodestatus,int[] cat,int[] jts,int[] jests,int[] bestvar,float[] tmissts,int[] ncts,
     float[] fill,int missfill,float code,float[] tnodewt,float[] outts,int idataout,int[] imax,
     int look,int lookcls,int[] nodexts,int isumout,int[][] mtab,String rf) throws Exception{   	
	     int m,j,n,jb,idummy;
	     zermr(qts,nclass,ntest);	     
	     for(m=1;m<=mdim;m++){
	    	 cat[m]=1;
	     }	     
	     if(labelts==1){
            for(n=1;n<=ntest;n++){
                    ncts[clts[n]]=ncts[clts[n]]+1;
            }
	     }
	     String line = "";
	     BufferedReader br_f1 = new BufferedReader(new FileReader(rf)); //  open(1,file='savedforest',status='old')
	     //BufferedReader br_f1 = new BufferedReader(new FileReader("C:\\MHMM\\rf.txt")); 
	     // START DOWN FOREST
	     for(jb=1;jb<=jbt;jb++){	    	 
	    	 	line = br_f1.readLine();
	    	 	line = line.replaceAll("\\s+","");
	    	 	ndbigtree = Integer.valueOf(line).intValue();
	            for(n=1;n<=ndbigtree;n++){	                	
	                line = br_f1.readLine(); 
	                line = line.replaceFirst("^\\s+", "");
	                String[] wds = line.split("\\s+");	                		
	                idummy = Integer.valueOf(wds[0]).intValue();
	                nodestatus[n] = Integer.valueOf(wds[1]).intValue();
	                System.out.println("nodestatus[n] "  + n + " " + Integer.valueOf(wds[1]).intValue());
	                bestvar[n] = Integer.valueOf(wds[2]).intValue();
	                treemap[1][n] = Integer.valueOf(wds[3]).intValue();
	                treemap[2][n] = Integer.valueOf(wds[4]).intValue();
	                nodeclass[n]  = Integer.valueOf(wds[5]).intValue();
	                xbestsplit[n] = Float.valueOf(wds[6]).floatValue();
	                tnodewt[n]    = Float.valueOf(wds[7]).floatValue();
	                for(j=1;j<=maxcat;j++){
	                	nbestcat[j][n] = Integer.valueOf(wds[7+j]).intValue();
	                }	                		
	                //read(1,*) idummy,nodestatus(n),bestvar(n),
	                //treemap(1,n),treemap(2,n),nodeclass(n),
	                //xbestsplit(n),tnodewt(n),(nbestcat(j,n),j=1,maxcat)
	             }	                
	             System.out.println("OOO ndbtree:" + ndbigtree);
	             testreelite(xts,mdim,treemap,nodestatus,xbestsplit,nrnodes,cat,jts,maxcat,nbestcat,nodexts);
	             for(n=1;n<=ntest;n++){
	                 qts[jts[n]][n]=qts[jts[n]][n]+tnodewt[nodexts[n]];
	             }
	             if(labelts==1){
	                  //if(jb%look==0&&jb<=jbt){
	                  if(jb<=jbt){
	                	    //delete errts variable from function, set it as global variable
	                      comperrts(qts,clts,ntest,nclass,tmissts,ncts,jests,labelts);
	                      if(lookcls==1){
	                         //write(*,'(i8,100f10.2)')
	                         //jb,100*errts,(100*tmissts(j),j=1,nclass)
                             System.out.println("jb=" + jb + " 100*errts=" + 100*errts + " " + 100*tmissts[1] + " " + 100*tmissts[2]);                               	        
	                      }else{
	                         System.out.println("jb=" + jb + " 100*errts=" + 100*errts);
	                      }
	                   }
	             }         	
	    } // end of jb

		///////////////////////////////////////////		         
		if(idataout==2){
			if(labelts==1){
				for(n=1;n<=ntest;n++){
				//  write(7,'(3i5,1000f10.3)') n,clts(n),jests(n),
				//  (qts(j,n),j=1,nclass),(xts(m,n),m=1,mdim)
				    System.out.println(n + " " + clts[n] + " " + jests[n]);
					for(j=1;j<=nclass;j++){
						System.out.print(qts[j][n] + " ");
					}
					//System.out.print(" ");
					for(m=1;m<=mdim;m++){
						System.out.print(xts[m][n] + " ");
					}
				}
			}else{
				for(n=1;n<=ntest;n++){
					System.out.println(n  + " " + jests[n]);
				    //write(7,'(3i5,1000f10.3)') n,jests(n),
					for(j=1;j<=nclass;j++){
						System.out.print(qts[j][n] + " ");
					}
					System.out.println("");
					for(m=1;m<=mdim;m++){
						System.out.print(xts[m][n] + " ");
					}
					System.out.println("");
				    //(qts(j,n),j=1,nclass),(xts(m,n),m=1,mdim)
				}
			}
		}               	                
		///////////////////////////////////////////		
		if (isumout==1){
	        if(labelts==1){
	           System.out.println("final error test %    " + 100*errts);
	           zerm(mtab,nclass,nclass);
	           for(n=1;n<=ntest;n++){
	                mtab[clts[n]][jests[n]]=mtab[clts[n]][jests[n]]+1;
	           }
	           System.out.println("Test set confusion matrix:");
	           System.out.println("     true class ");
	           System.out.println("");
	           //write(*,'(20i6)')  (j,j=1,nclass)
	           System.out.print("   ");
	           for(j=1;j<=nclass;j++){
	        	   System.out.print(j + "  ");
	           }
	           System.out.println("");
	           for(n=1;n<=nclass;n++){
	               System.out.print(n + " ");
	               for(j=1;j<=nclass;j++){
	            	   System.out.print(" " + mtab[j][n] + " ");
	               }
	               System.out.println("");
	        	   //write(*,'(20i6)')  n,(mtab(j,n),j=1,nclass)
	           }
	           System.out.println("");
	        }
		}	     
	    br_f1.close();
	    br_f1=null;
    }
    
	// fun2: makea, start makea	
	public static void makea(float[][] x,int mdim,int nsample,int[] cat,int[] isort,float[] v,int[][] asave,int[][] b,int mdimt,
		     	int[] msm,float[] v5,float[] v95,int maxcat){		
			// real x(mdim,nsample),v(nsample),v5(mdim),v95(mdim)
			// integer cat(mdim),isort(nsample),asave(mdim,nsample),
		    // &	b(mdim,nsample),msm(mdim)
			// integer mdim,nsample,mdimt,maxcat
			int k,mvar,n,n1,n2,ncat,jj;		
//		c	submakea constructs the mdim x nsample integer array a.
//		c	If there are less than 32,000 cases, this can be declared 
//		c	integer*2,otherwise integer*4. For each numerical variable 
//		c	with values x(m,n),n=1,...,nsample, the x-values are sorted 
//		c	from lowest to highest. Denote these by xs(m,n). 
//		c	Then asave(m,n) is the case  number in which 
//		c	xs(m,n) occurs. The b matrix is also constructed here. If the mth 
//		c	variable is categorical, then asave(m,n) is the category of the nth 
//		c	case number.  
//		c
//		c	input: x,cat
//		c	output: a,b
//		c	work: v,isort		
			for(k=1;k<=mdimt;k++){
			  mvar=msm[k];
			  if(cat[mvar]==1){
				for(n=1;n<=nsample;n++){
					v[n]=x[mvar][n];
					isort[n]=n;
				}				
				Y_quicksort(v,isort,1,nsample,nsample);				
				// this sorts the v(n) in ascending order. isort(n) is the 
				// case number of that v(n) nth from the lowest (assume 
				// the original case numbers are 1,2,...).  		
				n1=nint(0.05f*nsample);
				if(n1<1) n1=1;
				v5[mvar]=v[n1];
				n2=nint(0.95f*nsample);
				if(n2>nsample) n2=nsample;
				v95[mvar]=v[n2];
				for(n=1;n<=nsample-1;n++){
					n1=isort[n];
					n2=isort[n+1];
					asave[mvar][n]=n1;
					if(n==1) b[mvar][n1]=1;
					if (v[n]<v[n+1]){
						b[mvar][n2]=b[mvar][n1]+1;
					}else{
						b[mvar][n2]=b[mvar][n1];
					}
				}
				asave[mvar][nsample]=isort[nsample];
			  }else{
				for(ncat=1;ncat<=nsample;ncat++){
					jj=nint(x[mvar][ncat]);
					asave[mvar][ncat]=jj;
				}
			  }
			}		
	}	
	
	// fun3.1 new
	public static void moda(int[][] asave,int[][] a,int nuse,int nsample,int mdim,int[] cat,int maxcat,
		     	int[] ncase,int[] jin,int mdimt,int[] msm){
		
			// integer asave(mdim,nsample),a(mdim,nsample),cat(mdim),
		    // &	jin(nsample),ncase(nsample),msm(mdim)
			// integer nuse,nsample,mdim,mdimt,maxcat
			int n,jj,m,k,nt,j;		
			// copy rows msm(1),...,msm(mdimt) of asave into the same 
			// rows of a
		
			nuse=0;
			for(n=1;n<=nsample;n++){
				for(k=1;k<=mdimt;k++){
					m=msm[k];
					a[m][n]=asave[m][n];
				}
				if(jin[n]>=1) nuse=nuse+1;
			}
			
			for(jj=1;jj<=mdimt;jj++){
				m=msm[jj];
				k=1;
				nt=1;
				if(cat[m]==1){
	   AA37:		for(n=1;n<=nsample;n++){
						if(k>nsample) break AA37;
						    if(jin[a[m][k]]>=1){
								a[m][nt]=a[m][k];
								k=k+1;
						    }else{
				AA28:			for(j=1;j<=nsample-k;j++){
									if(jin[a[m][k+j]]>=1){
										a[m][nt]=a[m][k+j];
										k=k+j+1;
										break AA28;
										//goto 28
									}
								}
							}
			//8				continue
							nt=nt+1;
							if(nt>nuse) break AA37;   // goto 37
					}
		//37			continue
				}
			}
			
			if(maxcat>1){
				k=1;
				nt=1;
		AA85:	for(n=1;n<=nsample;n++){
					if(jin[k]>=1){
						ncase[nt]=k;
						k=k+1;
					}
					else{
				AA58:	for(jj=1;jj<=nsample-k;jj++){
							if(jin[k+jj]>=1){
								ncase[nt]=k+jj;
								k=k+jj+1;
								break AA58;
								//goto 58
							}
						}
					}
		//58			continue
					nt=nt+1;
					if(nt>nuse) break AA85; //goto 85
				}
		//85		continue
			}
	}
	
//	// fun3 moda£º this function should be re-written
//	public static void old_moda(int[][] asave,int[][] a,int nuse,int nsample,int mdim,int[] cat,int maxcat,
//		     	int[] ncase,int[] jin,int mdimt,int[] msm){		
//			// integer asave(mdim,nsample),a(mdim,nsample),cat(mdim),
//		    //	jin(nsample),ncase(nsample),msm(mdim)
//			//integer nuse,nsample,mdim,mdimt,maxcat
//			int n,jj,m,k,nt,j;		
//			// copy rows msm(1),...,msm(mdimt) of asave into the same 
//			// rows of a		
//			nuse=0;
//			for(n=1;n<=nsample;n++){
//				for(k=1;k<=mdimt;k++){
//					m=msm[k];
//					a[m][n]=asave[m][n];
//				}
//				if(jin[n]>=1) nuse=nuse+1;
//			}
//    AA37:	for(jj=1;jj<=mdimt;jj++){
//				m=msm[jj];
//				k=1;
//				nt=1;
//				if(cat[m]==1){
//			AA28:for(n=1;n<=nsample;n++){
//					if(k>nsample) continue AA37;
//					    if(jin[a[m][k]]>=1){
//							a[m][nt]=a[m][k];
//							k=k+1;
//					    }
//						else{
//							for(j=1;j<=nsample-k;j++){
//								if(jin[a[m][k+j]]>=1){
//									a[m][nt]=a[m][k+j];
//									k=k+j+1;
//									continue AA28;
//								}
//							}
//						}
//		//28			continue
//						nt=nt+1;
//						if(nt>nuse) continue AA37; //goto 37
//			 	    }
//		//37			//continue
//				}
//			}
//			if(maxcat>1){
//				k=1;
//				nt=1;
//		AA85:   for(n=1;n<=nsample;n++){
//					if(jin[k]>=1){
//						ncase[nt]=k;
//						k=k+1;
//					}else{
//				AA58:	for(jj=1;jj<=nsample-k;jj++){
//							if(jin[k+jj]>=1){
//								ncase[nt]=k+jj;
//								k=k+jj+1;
//								continue AA58;
//							}
//						}
//					}
//		//58			continue
//					nt=nt+1;
//					if(nt>nuse)  continue AA85; //goto 85
//				}
//		//85		continue
//			}
//	}
	
	// fun4: buildtree
	// delete ndbigtree in the Function and set it as global variable
	public static void buildtree(int[][] a,int[][] b,int[] cl,int[] cat,int mdim,int nsample,int nclass,int[][] treemap,
			      int[] bestvar,int[] bestsplit,int[] bestsplitnext,float[] dgini,
			      int[] nodestart,float[][] classpop,float[] tclasspop,float[][] tclasscat,int[] ta,
			      int[] idmove,int ndsize,int[] ncase,int[] parent,int mtry,
			      float[] win,float[] wr,float[] wl,int nuse,int[] kcat,int[] ncatsplit,float[] xc,float[] dn,
			      float[] cp,float[] cm,int[][] nbestcat,int[] msm,int mdimt){			
				// Buildtree consists of repeated calls to findbestsplit and movedata.
				// Findbestsplit does just that--it finds the best split of the current 
				// node. Movedata moves the data in the split node right and left so 
				// that the data corresponding to each child node is contiguous.  		
				// The buildtree bookkeeping is different from that in Friedman's 
				// original CART program: 
				//	ncur is the total number of nodes to date
				//	nodestatus(k)=1 if the kth node has been split.
				//	nodestatus(k)=2 if the node exists but has not yet been split
				//		  and=-1 of the node is terminal.  
				// A node is terminal if its size is below a threshold value, or if it 
				// is all one class,or if all the x-values are equal. If the current 
				// node k is split,then its children are numbered ncur+1 (left), and 
				// ncur+2(right),ncur increases to ncur+2 and the next node to be split 
				// is numbered k+1. When no more nodes can be split,buildtree
				// returns to the main program.
			
				// integer cl(nsample),cat(mdim),ncatsplit(maxcat),
			    //   treemap(2,nrnodes),bestvar(nrnodes),nodeclass(nrnodes),
			    //  	bestsplit(nrnodes),nodestatus(nrnodes),ta(nsample),
			    //      nodepop(nrnodes),nodestart(nrnodes),idmove(nsample),
			    //   	bestsplitnext(nrnodes),ncase(nsample),parent(nrnodes),
			    //   	kcat(maxcat),msm(mdim),iseed,ncsplit,ncmax
			    	
				// integer a(mdim,nsample),b(mdim,nsample)
				// real tclasspop(nclass),classpop(nclass,nrnodes),
			    //  tclasscat(nclass,maxcat),win(nsample),wr(nclass),
			    //  wl(nclass),dgini(nrnodes),xc(maxcat),dn(maxcat),
			    // 	cp(maxcat),cm(maxcat)
				// integer mdim,nsample,nclass,nrnodes,ndsize,
			    //  mtry,ndbigtree,nuse,maxcat,mdimt,j,ncur,
			    //  kbuild,ndstart,ndend,jstat
			
				int ncur,kbuild;
				int i,j,k;
				int lcat,kn,n,nc;
				float popt1,popt2,pp;	
				zerv(nodestatus,nrnodes);
				zerv(nodestart,nrnodes);
				zerv(nodepop,nrnodes);
				zermr(classpop,nclass,nrnodes);
				zerm(treemap,2,nrnodes);
				zerm(nbestcat,maxcat,nrnodes);				
				for(j=1;j<=nclass;j++){
					classpop[j][1]=tclasspop[j];
				}			
				ncur=1;
				nodestart[1]=1;
				nodepop[1]=nuse;
				nodestatus[1]=2;			
				//start main loop		
	  		    //do 30 kbuild=1,nrnodes
       AA30:	for(kbuild=1;kbuild<=nrnodes;kbuild++){
					if(kbuild>ncur)  break AA30; //goto 50
					if(nodestatus[kbuild]!=2) continue AA30; // goto 30					
					// initialize for next call to findbestsplit		
					ndstart=nodestart[kbuild];
					ndend=ndstart+nodepop[kbuild]-1;
					for(j=1;j<=nclass;j++){
						tclasspop[j]=classpop[j][kbuild];
					}
					jstat=0;						
					findbestsplit(a,b,cl,mdim,nsample,nclass,cat,
			     		tclasspop,tclasscat,
			     		ncase,mtry,win,wr,wl,kcat,ncatsplit,xc,dn,
			     		cp,cm,msm,mdimt);			
					if(jstat==1){
						nodestatus[kbuild]=-1;
						continue AA30;
						//goto 30
					}else{
						bestvar[kbuild]=msplit;
						dgini[kbuild]=decsplit;					
						if (cat[msplit]==1){
							//continuous
							bestsplit[kbuild]=a[msplit][nbest];
							bestsplitnext[kbuild]=a[msplit][nbest+1];
						}else{
							//categorical
							lcat=cat[msplit];
							for(i=1;i<=lcat;i++){
								nbestcat[i][kbuild]=ncatsplit[i];
							}
						}
					}
					// delete ndendl  ndstart ndend   , set them as global variables
					movedata(a,ta,mdim,nsample,idmove,ncase,msplit,cat,nbest,ncatsplit,maxcat,mdimt,msm);		
					//leftnode no.=ncur+1,rightnode no.=ncur+2.		
					nodepop[ncur+1]=ndendl-ndstart+1;
					nodepop[ncur+2]=ndend-ndendl;
					nodestart[ncur+1]=ndstart;
					nodestart[ncur+2]=ndendl+1;
					//find class populations in both nodes
					for(n=ndstart;n<=ndendl;n++){
						nc=ncase[n];
						j=cl[nc];
						classpop[j][ncur+1]=classpop[j][ncur+1]+win[nc];
					}
					for(n=ndendl+1;n<=ndend;n++){
						nc=ncase[n];
						j=cl[nc];
						classpop[j][ncur+2]=classpop[j][ncur+2]+win[nc];
					}			
					//check on nodestatus		
					nodestatus[ncur+1]=2;
					nodestatus[ncur+2]=2;
					if(nodepop[ncur+1]<=ndsize) nodestatus[ncur+1]=-1;
					if(nodepop[ncur+2]<=ndsize) nodestatus[ncur+2]=-1;
					popt1=0;
					popt2=0;
					for(j=1;j<=nclass;j++){
						popt1=popt1+classpop[j][ncur+1];
						popt2=popt2+classpop[j][ncur+2];
					}
					for(j=1;j<=nclass;j++){
						if(Math.abs(classpop[j][ncur+1]-popt1)<0.00001f) nodestatus[ncur+1]=-1;
						if(Math.abs(classpop[j][ncur+2]-popt2)<0.00001f) nodestatus[ncur+2]=-1;
					}
					treemap[1][kbuild]=ncur+1;
					treemap[2][kbuild]=ncur+2;
					parent[ncur+1]=kbuild;
					parent[ncur+2]=kbuild;
					nodestatus[kbuild]=1;
					ncur=ncur+2;
					if(ncur>=nrnodes) break AA30;//goto 50		
				} //30	continue
			//50	continue				
			
				System.out.println("Oracle: ndbigtree: "+ ndbigtree);
				ndbigtree=nrnodes;
				System.out.println("Oracle: after ndbigtree: "+ ndbigtree);
				System.out.println("CCC" + " " + ndbigtree + " nrnodes" + nrnodes);
				//do k=nrnodes,1,-1
				System.out.println("00000000000000000000000 nrnodes: " + nrnodes);
				for(k=nrnodes;k>=1;k--){   // do k=nrnodes,1,-1
					if(nodestatus[k]==0) ndbigtree=ndbigtree-1;
					if(nodestatus[k]==2) nodestatus[k]=-1;
				}
				for(kn=1;kn<=ndbigtree;kn++){
					if(nodestatus[kn]==-1){
						pp=0;
						for(j=1;j<=nclass;j++){
							if(classpop[j][kn]>pp){
								nodeclass[kn]=j;
								pp=classpop[j][kn];
							}
						}
					}
				}				
	}
	
	// fun5: findbestsplit: jstat should be removed and set as a global variable	
	// delete decsplit and as a global variable
	public static void findbestsplit(int[][] a,int[][] b,int[] cl,int mdim,int nsample,int nclass,int[] cat,
		        float[] tclasspop,float[][] tclasscat,
		        int[] ncase,int mtry,float[] win,float[] wr,float[] wl,int[] kcat,int[] ncatsplit,float[] xc,float[] dn,
		        float[] cp,float[] cm,int[] msm,int mdimt){	
			//For the best split,msplit is the variable split on. decsplit is the dec. in impurity.
			//If msplit is numerical,nsplit is the case number of value of msplit split on,
			//and nsplitnext is the case number of the next larger value of msplit.  If msplit is
			//categorical,then nsplit is the coding into an integer of the categories going left.		
			//integer a(mdim,nsample),b(mdim,nsample),iseed,ncsplit,ncmax		
			//integer cl(nsample),cat(mdim),ncase(nsample),msm(mdim),
		    // &	kcat(maxcat),ncatsplit(maxcat),icat(32)		
			//integer mdim,nsample,nclass,ndstart,ndend,mdimt,mtry,
			int mv,mvar,nc;
			int i,j,k;
		    // &	msplit,nbest,jstat,maxcat,i,j,k,mv,mvar,nc,
			int nsp;
		    // &	lcat,nnz,nhit,ncatsp,nsp		
			//real tclasspop(nclass),tclasscat(nclass,maxcat),
		    // &	win(nsample),wr(nclass),wl(nclass),xc(maxcat),
		    // &	dn(maxcat),cp(maxcat),cm(maxcat)		
			float pno,pdo,rrn,rrd,rln,rld,u;
		    float crit0,critmax,crit;		
//			float randomu;		
			//external unpack		
			//compute initial values of numerator and denominator of Gini			
			pno=0;
			pdo=0;
			for(j=1;j<=nclass;j++){
				pno=pno+tclasspop[j]*tclasspop[j];
				pdo=pdo+tclasspop[j];
			}
			crit0=pno/pdo;
			jstat=0;			
			//start main loop through variables to find best split			
			critmax=-999999;		
			for(mv=1;mv<=mtry;mv++){
				k=(int)(mdimt*randomu())+1;
				//k=3;
				mvar=msm[k];
				if(cat[mvar]==1){
					//it's not a categorical variable:		
					rrn=pno;
					rrd=pdo;
					rln=0;
					rld=0;
					zervr(wl,nclass);
					for(j=1;j<=nclass;j++){
						wr[j]=tclasspop[j];
					}
					for(nsp=ndstart;nsp<=ndend-1;nsp++){
						nc=a[mvar][nsp];
						u=win[nc];
						k=cl[nc];
						rln=rln+u*(2*wl[k]+u);
						rrn=rrn+u*(-2*wr[k]+u);
						rld=rld+u;
						rrd=rrd-u;
						wl[k]=wl[k]+u;
						wr[k]=wr[k]-u;
		  				if(b[mvar][nc]<b[mvar][a[mvar][nsp+1]]){
							if(amin1(rrd,rld)>1.0e-5){
								crit=(rln/rld)+(rrn/rrd);
								if (crit>critmax){
									nbest=nsp;
									critmax=crit;
									msplit=mvar;									
								}
							}
		  				}
					}	
				}
			} // end if		
			decsplit=critmax-crit0;
			if(critmax<-1.0e6) jstat=1;
	}	
	
	// fun6: catmaxr
	public static void catmaxr(int ncsplit,float[][] tclasscat,float[] tclasspop,int[] icat,
		     int nclass,int lcat,int maxcat,int[] ncatsplit,float critmax,float pdo,int nhit,int iseed){
			// this routine takes the best of ncsplit random splits		
			//  	real tclasscat(nclass,maxcat),tclasspop(nclass),tmpclass(100)
			//		real critmax,pdo
			//		integer icat(32),iseed
			//		integer ncsplit,nclass,lcat,maxcat,ncatsplit(maxcat),nhit		
			float[] tmpclass = new float[100+1];
			int irbit;
			float pln,pld,prn,tdec;
			int j,n,i,k;
			nhit=0;
			for(n=1;n<=ncsplit;n++){
				// generate random split
				for(k=1;k<=lcat;k++){
					icat[k]=irbit(iseed);	// !icat(k) is bernouilli
				}
				for(j=1;j<=nclass;j++){
					tmpclass[j]=0;
					for(k=1;k<=lcat;k++){
						if(icat[k]==1){
							tmpclass[j]=tmpclass[j]+tclasscat[j][k];
						}
					}
				}
				pln=0;
				pld=0;
				for(j=1;j<=nclass;j++){
					pln=pln+tmpclass[j]*tmpclass[j];
					pld=pld+tmpclass[j];
				}
				prn=0;
				for(j=1;j<=nclass;j++){
					tmpclass[j]=tclasspop[j]-tmpclass[j];
					prn=prn+tmpclass[j]*tmpclass[j];
				}
				tdec=(pln/pld)+(prn/(pdo-pld));
				if (tdec>critmax){
					critmax=tdec;
					nhit=1;
					for(k=1;k<=lcat;k++){
						ncatsplit[k]=icat[k];
					}
				}
			}
	}
	
	// fun7: catmax
	public static void catmax(float pdo,float[][] tclasscat,float[] tclasspop,int nclass,int lcat,
		       int ncatsp,float critmax,int nhit,int maxcat){		
			// this finds the best split of a categorical variable
			// with lcat categories and nclass classes, where 
			// tclasscat(j,k) is the number of cases in
		    // class j with category value k. The method uses an 
			// exhaustive search over all partitions of the category 
			// values. For the two class problem,there is a faster 
			// exact algorithm. If lcat.ge.10,the exhaustive search
			// gets slow and there is a faster iterative algorithm.		
			// real tclasscat(nclass,maxcat),tclasspop(nclass),tmpclass(100)
			// integer icat(32),n,lcat,nhit,l,j,ncatsp,nclass,maxcat
			// real critmax,pdo,pln,pld,prn,tdec
			float[] tmpclass = new float[100+1];
			int[] icat = new int[32+1];
			int n,l,j;
			float pln,pld,prn,tdec;			
			//external unpack			
			nhit=0;
			//do n=1,(2**(lcat-1))-1
			int tmp = (int)Math.pow(2, lcat-1)-1;
			for(n=1;n<=tmp;n++){
				unpack(lcat,n,icat);
				for(j=1;j<=nclass;j++){
					tmpclass[j]=0;
					for(l=1;l<=lcat;l++){
						if(icat[l]==1){
							tmpclass[j]=tmpclass[j]+tclasscat[j][l];
						}
					}
				}
				pln=0;
				pld=0;
				for(j=1;j<=nclass;j++){
					pln=pln+tmpclass[j]*tmpclass[j];
					pld=pld+tmpclass[j];
				}
				prn=0;
				for(j=1;j<=nclass;j++){
					tmpclass[j]=tclasspop[j]-tmpclass[j];
					prn=prn+tmpclass[j]*tmpclass[j];
				}
				tdec=(pln/pld)+(prn/(pdo-pld));
				if (tdec>critmax){
					critmax=tdec;
					ncatsp=n;
					nhit=1;
				}
			}
	}
	
	// fun8	
	public static void movedata(int[][] a,int[] ta,int mdim,int nsample,int[] idmove,
			       int[] ncase,int msplit,int[] cat,int nbest,int[] ncatsplit,int maxcat,int mdimt,int[] msm){		
				// movedata is the heart of the buildtree construction. 
				// Based on the best split the data corresponding to the 
				// current node is moved to the left if it belongs to the 
				// left child and right if it belongs to the right child.		
				// integer a(mdim,nsample),ta(nsample),idmove(nsample),
			    //	  ncase(ndend),cat(mdim),ncatsplit(maxcat),msm(mdimt)
				//integer mdim,ndstart,ndend,nsample,msplit,nbest,ndendl
				int nsp,nc,ms,msh,k,n,ih;			
				// compute idmove=indicator of case nos. going left		
				if(cat[msplit]==1){
					for(nsp=ndstart;nsp<=nbest;nsp++){
						nc=a[msplit][nsp];
						idmove[nc]=1;
					}
					for(nsp=nbest+1;nsp<=ndend;nsp++){
						nc=a[msplit][nsp];
						idmove[nc]=0;
					}
					ndendl=nbest;
				}else{
					ndendl=ndstart-1;
					for(nsp=ndstart;nsp<=ndend;nsp++){
						nc=ncase[nsp];
						if(ncatsplit[a[msplit][nc]]==1){
							idmove[nc]=1;
							ndendl=ndendl+1;
						}else{
							idmove[nc]=0;
						}
					}
				 }			
				// shift case. nos. right and left for numerical variables.				
				for(ms=1;ms<=mdimt;ms++){
					msh=msm[ms];
					if(cat[msh]==1){
						k=ndstart-1;
						for(n=ndstart;n<=ndend;n++){
							ih=a[msh][n];
							if(idmove[ih]==1){
								k=k+1;
								ta[k]=a[msh][n];
							}
						}
						for(n=ndstart;n<=ndend;n++){
							ih=a[msh][n];
							if (idmove[ih]==0){ 
								k=k+1;
								ta[k]=a[msh][n];
							}
						}
						for(k=ndstart;k<=ndend;k++){
							a[msh][k]=ta[k];
						}
					}
				}
				// compute case nos. for right and left nodes.			
				if(cat[msplit]==1){
					for(n=ndstart;n<=ndend;n++){
						ncase[n]=a[msplit][n];
					}
				}else{
					k=ndstart-1;
					for(n=ndstart;n<=ndend;n++){
						if(idmove[ncase[n]]==1){
							k=k+1;
							ta[k]=ncase[n];
						}
					}
					for(n=ndstart;n<=ndend;n++){
						if (idmove[ncase[n]]==0){
							k=k+1;
							ta[k]=ncase[n];
						}
					}
					for(k=ndstart;k<=ndend;k++){
						ncase[k]=ta[k];
					}
				}
	}
	
	// fun9: xtranslate
	public static void xtranslate(float[][] x,int mdim,int nrnodes,int nsample,
		     int[]  bestsplit,int[] bestsplitnext,float[] xbestsplit,int[] nodestatus,int[] cat){		
			// xtranslate takes the splits on numerical variables and translates them
			// back into x-values. It also unpacks each categorical split into a 32-
			// dimensional vector with components of zero or one--a one indicates 
			// that the corresponding category goes left in the split.		
			//integer cat(mdim),bestvar(nrnodes),bestsplitnext(nrnodes),
		    //	nodestatus(nrnodes),bestsplit(nrnodes)
			// real x(mdim,nsample),xbestsplit(nrnodes)
			int k,m;
			for(k=1;k<=ndbigtree;k++){
				if(nodestatus[k]==1){
					m=bestvar[k];
					if(cat[m]==1){
						xbestsplit[k]=(x[m][bestsplit[k]] + x[m][bestsplitnext[k]])/2;
					}else{
						xbestsplit[k]=1.0f*(bestsplit[k]);
					}
				}
			}
	}	
	
	// fun10: getweights
	public static void getweights(float[][] x,int nsample,int mdim,int[][] treemap,int[] nodestatus,
		        float[] xbestsplit,int[] bestvar,int nrnodes,
		     	int[] cat,int maxcat,int[][] nbestcat,int[] jin,float[] win,float[] tw,float[] tn,float[] tnodewt){		
			// real x(mdim,nsample),xbestsplit(nrnodes),
		    //  win(nsample),tw(ndbigtree),tn(ndbigtree),tnodewt(ndbigtree)
			// integer nsample,mdim,nrnodes,ndbigtree,maxcat
			// integer treemap(2,nrnodes),bestvar(nrnodes),
		    // cat(mdim),nodestatus(nrnodes),jin(nsample)		
			//integer nbestcat(maxcat,nrnodes)			
			int jcat,n,kt,k,m;			
			zervr(tw,ndbigtree);			
			zervr(tn,ndbigtree);			
AA100:		for(n=1;n<=nsample;n++){
				if(jin[n]>=1){
					kt=1;
					for(k=1;k<=ndbigtree;k++){
						if (nodestatus[kt]==-1){
							tw[kt]=tw[kt]+win[n];
							tn[kt]=tn[kt]+jin[n];
							continue AA100;
							//goto 100
						}
						m=bestvar[kt];
						if (cat[m]==1){
							if(x[m][n]<=xbestsplit[kt]){ 
								kt=treemap[1][kt];
							}else{
								kt=treemap[2][kt];
							}
						}else{
							jcat=nint(x[m][n]);
							if (nbestcat[jcat][kt]==1){
								kt=treemap[1][kt];
							}else{
								kt=treemap[2][kt];
							}
						}
					}
		//100 continue
				}
			}			
			for(n=1;n<=ndbigtree;n++){
				if(nodestatus[n]==-1) tnodewt[n]=tw[n]/tn[n];
			}
	}	
	
	// fun11: testreebag, delete ndbigtree and set it as global variable
	public static void testreebag(float[][] x,int nsample,int mdim,int[][] treemap,int[] nodestatus,
		     float[]  xbestsplit,int nrnodes,int[] kpop,
		     int[] cat,int[] jtr,int[] nodextr,int maxcat,int[][] nbestcat,float[] rpop,float[] dgini,float[] tgini,
		     int[] jin,float[] wtx){		
			//	predicts the class of all objects in x		
			//	input:
			//	real x(mdim,nsample),xbestsplit(nrnodes),dgini(nrnodes),
			//     &	rpop(nrnodes),wtx(nsample)
			//	integer treemap(2,nrnodes),bestvar(nrnodes),
			//     &  nodeclass(nrnodes),cat(mdim),nodestatus(nrnodes),
			//     &	nbestcat(maxcat,nrnodes),jin(nsample),kpop(nrnodes)		
			//	integer nsample,mdim,nrnodes,ndbigtree,maxcat		
			//	output:
			//real tgini(mdim)
			//integer jtr(nsample),nodextr(nsample)		
			//	local
			int n,kt,k,m,jcat;
			zerv(jtr,nsample);
			zerv(nodextr,nsample);
			zervr(tgini,mdim);
			zervr(rpop,nrnodes);
			zerv(kpop,nrnodes);		
			n=1;			
		    //903	kt=1
AAdo:   	do{							
				kt=1;			
				for(k=1;k<=ndbigtree;k++){
					if(nodestatus[kt]==-1){
						jtr[n]=nodeclass[kt];
						nodextr[n]=kt;
						if(jin[n]>0){
							rpop[kt]=rpop[kt]+jin[n];
							kpop[kt]=kpop[kt]+1;
						}
						{
							n = n + 1;
							if(n<=nsample){
								continue AAdo;
							}else{
								break AAdo;
							}
						}
						//goto 100
					}
					m=bestvar[kt];
					tgini[m]=tgini[m]+dgini[kt];
					if(cat[m]==1){
						if(x[m][n]<=xbestsplit[kt]){ 
							kt=treemap[1][kt];
						}else{
							kt=treemap[2][kt];
						}
					}
					if(cat[m]>1){
						jcat=nint(x[m][n]);
						if (nbestcat[jcat][kt]==1){
							kt=treemap[1][kt];
						}else{
							kt=treemap[2][kt];
						}
					}
				} // !k
				n=n+1;
			//100	n=n+1		
			}while(n<=nsample); //goto 903
			
			for(m=1;m<=mdim;m++){
				tgini[m]=tgini[m]/nsample;
			}			
	}
	
	// fun12: testreelite
	// can not delete the variable ndbigtree and set it as global variable
	public static void testreelite(float[][] xts,int mdim,int[][] treemap,int[] nodestatus,
		     float[] xbestsplit,int nrnodes,int[] cat,int[] jts,int maxcat,int[][] nbestcat,int[] nodexts){			
//	input:
//	        real xts(mdim,ntest),xbestsplit(nrnodes)
//	        integer treemap(2,nrnodes),bestvar(nrnodes),nodeclass(nrnodes),
//	        cat(mdim),nodestatus(nrnodes),nbestcat(maxcat,nrnodes)
//	        integer ntest,mdim,nrnodes,ndbigtree,maxcat
//  output: integer jts(ntest),nodexts(ntest)		
		int n,kt,k,m,jcat;
		n=1;
		do{
		//903     
		    kt=1;
		    KK:for(k=1;k<=ndbigtree;k++){		    	 
				 if(nodestatus[kt]==-1){
				       jts[n]=nodeclass[kt];
				       nodexts[n]=kt;			
				       //goto 100
				       break KK;				                         
				 }
				 m=bestvar[kt];
				 if(cat[m]==1){
				     if (xts[m][n]<=xbestsplit[kt]){
				           kt=treemap[1][kt];
				      }else{
				           kt=treemap[2][kt];
				      }
				  }else{
					  
				  }
				 
				  if(cat[m]>1){
				       jcat=nint((xts[m][n])); // jcat=nint(xts(m,n))  Please check this function
				       if(nbestcat[jcat][kt]==1){
				            kt=treemap[1][kt];
				       }else{
				            kt=treemap[2][kt];
				       }
				  }
		   }//enddo !k
		   //100     n=n+1
		   n=n+1; 
		}while(n<=ntest);
		
		// output: integer jts(ntest),nodexts(ntest)
		
		for(int i=1;i<=ntest;i++){
			System.out.println("jts[i]: " + jts[i]);
		}
		
		for(int i=1;i<=ntest;i++){
			// System.out.println("nodexts[i]: " + nodexts[i]);
		}
		
	}
	
	// fun13: testreeimp
	public static void testreeimp(float[][] x,int nsample,int mdim,int[] joob,int[] pjoob,int nout,int mr,
		     	int[][] treemap,int[] nodestatus,float[] xbestsplit,int[] bestvar,int[] nodeclass,int nrnodes,
		     	int[] cat,int[] jvr,int[] nodexvr,int maxcat,int[][] nbestcat){		
			// predicts the class of out-of-bag-cases for variable importance
			// also computes nodexvr		
			// input:
			//	real x(mdim,nsample),xbestsplit(nrnodes)
			//	integer treemap(2,nrnodes),bestvar(nrnodes),nodeclass(nrnodes),
			//     &	cat(mdim),nodestatus(nrnodes),nbestcat(maxcat,nrnodes),
			//     &	joob(nout),pjoob(nout)
			//	integer nsample,mdim,nout,mr,nrnodes,ndbigtree,maxcat		
			// output: integer jvr(nout),nodexvr(nout)		
			int n,kt,k,m,jcat;
			float xmn;		
			permobmr(joob,pjoob,nout);
			n=1;
AAdo:		do{
		//904	kt=1;
			kt=1;
			for(k=1;k<=ndbigtree;k++){
				if (nodestatus[kt]==-1){
					jvr[n]=nodeclass[kt];
					nodexvr[n]=kt;
					{						  
						n=n+1;
						continue AAdo;
						// equal to goto 100
					}
					//goto 100					
				}
				m=bestvar[kt];
				if(m==mr){
					xmn=x[m][pjoob[n]]; // ! permuted value
				}else{
					xmn=x[m][joob[n]];
				}
				if (cat[m]==1){
					if (xmn<=xbestsplit[kt]){ 
						kt=treemap[1][kt];
					}else{
						kt=treemap[2][kt];
					}
				}
				if(cat[m]>1){
					jcat=nint(xmn);
					if (nbestcat[jcat][kt]==1){
						kt=treemap[1][kt];
					}else{
						kt=treemap[2][kt];
					}
				}
			} // !k
			n=n+1;
		//100	n=n+1
			} while(n<=nout);//if(n.le.nout) goto 904			
	}
	
	// fun14: permobmr
	public static void permobmr(int[] joob,int[] pjoob,int nout){	
			// randomly permute the elements of joob and put them in pjoob	
			// input:
			// integer joob(nout),nout
			// output:
			// integer pjoob(nout)
			int j,k,jt;
			float rnd,randomu;	
			for(j=1;j<=nout;j++){
				pjoob[j]=joob[j];
			}
			j=nout;		
			do{
				rnd=randomu();
				k=(int)(j*rnd);
				if(k<j) k=k+1;
				//switch j and k
				jt=pjoob[j];
				pjoob[j]=pjoob[k];
				pjoob[k]=jt;
				j=j-1;
			}while(j>1); // go to 11	
	}
	
	// fun15: comperrtr
	public static void comperrtr(float[][] q,int[] cl,int nsample,int nclass,float[] tmiss,int[] nc,int[] jest,int[] out){			
			// integer cl(nsample),nc(nclass),jest(nsample),out(nsample)
			// real q(nclass,nsample),tmiss(nclass),errtr
			// integer nsample,nclass
			float cmax,ctemp;
			int n,j,jmax=-1;
			zervr(tmiss,nclass);
			errtr=0;
			for(n=1;n<=nsample;n++){
				cmax=0;
				if(out[n]>0){
					for(j=1;j<=nclass;j++){
						ctemp=q[j][n]/out[n];
						if(ctemp>cmax){
							jmax=j;
							cmax=ctemp;
						}
					}
				}else{
					jmax=1;
				}				
				jest[n]=jmax;
				if(jmax!=cl[n]){
					tmiss[cl[n]]=tmiss[cl[n]]+1;
					errtr=errtr+1;
				}
			}
			errtr=errtr/nsample;
			for(j=1;j<=nclass;j++){
				tmiss[j]=tmiss[j]/nc[j];
			}
	}
	
	// fun16: comperrts:  errts: modified to global variable    
	// q(j,n): the proportion of votes for the jth class out of the total of all the votes for the nth case 
    public static void comperrts(float[][] qts,int[] clts,int ntest,int nclass,float[] tmissts,int[] ncts,int[] jests,int labelts){     
        //integer clts(ntest),ncts(nclass),jests(ntest)
        //real qts(nclass,ntest),tmissts(nclass)
        //integer ntest,nclass,labelts         
    	//real errts,cmax
   	 	float cmax;
        int n,j,jmax=0;
        zervr(tmissts,nclass);
        errts=0;
        for(n=1;n<=ntest;n++){
              cmax=0;
              for(j=1;j<=nclass;j++){
                   if(qts[j][n]>cmax){
                       jmax=j;
                       cmax=qts[j][n];
                   }
              }
              jests[n]=jmax;
              if(labelts==1){
                   if(jmax!=clts[n]){
                       tmissts[clts[n]]=tmissts[clts[n]]+1;
                       errts=errts+1;
                   }
              }
        }         
        if(labelts==1){
              errts=errts/ntest;
              for(j=1;j<=nclass;j++){
                  tmissts[j]=tmissts[j]/ncts[j];
              }
        }            	 
    }
    
	// fun17: createclass
	public static void createclass(float[][] x,int[] cl,int ns,int nsample,int mdim){
		//real x(mdim,nsample)
		//integer cl(nsample)
		//integer ns,nsample,mdim
		float randomu;
		int n,m,k;	
		for(n=1;n<=ns;n++){
			cl[n]=1;
		}		
		for(n=ns+1;n<=nsample;n++){
			cl[n]=2;
		}
		for(n=ns+1;n<=nsample;n++){
			for(m=1;m<=mdim;m++){
				k=(int)(randomu()*ns);
				if(k<ns) k=k+1;
				x[m][n]=x[m][k];
			}
		}
	}
	
	// fun18: varimp
	public static void varimp(float[][] x,int nsample,int mdim,int[] cl,int nclass,int[] jin,int[] jtr,int impn,
			     	int interact,int[] msm,int mdimt,float[] qimp,float[][] qimpm,float[] avimp,float[] sqsd,
			     	int[][] treemap,int[] nodestatus,float[] xbestsplit,int[] bestvar,int[] nodeclass,int nrnodes,
			     	int[] cat,int[] jvr,int[] nodexvr,int maxcat,int[][] nbestcat,float[] tnodewt,
			     	int[] nodextr,int[] joob,int[] pjoob,int[] iv){		
					// real x(mdim,nsample),xbestsplit(nrnodes),
			     	// avimp(mdim),sqsd(mdim),
			     	// qimp(nsample),qimpm(nsample,mdim),tnodewt(nrnodes)
					// integer cl(nsample),jin(nsample),jtr(nsample),
			     	// msm(mdimt),treemap(2,nrnodes),nodestatus(nrnodes),
			     	// bestvar(nrnodes),nodeclass(nrnodes),cat(mdim),jvr(nsample),
			     	// nodexvr(nsample),nbestcat(maxcat,nrnodes),
			     	// nodextr(nsample),joob(nsample),pjoob(nsample),iv(mdim)
					//integer nsample,mdim,nrnodes,ndbigtree,maxcat,
			     	// nclass,impn,interact,mdimt
					int nout,mr,nn,n,jj,k;
					float right,rightimp;		
					nout=0;
					right=0;				
					for(n=1;n<=nsample;n++){
						if(jin[n]==0){
							//	this case is out-of-bag
							//	update count of correct oob classifications
							//	(jtr(n)=cl(n) if case n is correctly classified)
							if(jtr[n]==cl[n]) right=right+tnodewt[nodextr[n]];
							//	nout=number of obs out-of-bag for THIS tree
							nout=nout+1;
							joob[nout]=n;
						}
					}						
					if(impn==1){
						for(n=1;n<=nout;n++){
							nn=joob[n];
							if(jtr[nn]==cl[nn]){
								qimp[nn]=qimp[nn]+tnodewt[nodextr[nn]]/nout;
							}
						}
					}
					zerv(iv,mdim);
					for(jj=1;jj<=ndbigtree;jj++){
						//iv(j)=1 if variable j was used to split on
						if(nodestatus[jj]!=-1) iv[bestvar[jj]]=1;
					}
					for(k=1;k<=mdimt;k++){
						mr=msm[k];
						//choose only those that used a split on variable mr 
						if(iv[mr]==1){
							testreeimp(x,nsample,mdim,joob,pjoob,nout,mr,
				     			treemap,nodestatus,xbestsplit,bestvar,nodeclass,nrnodes,
				     			cat,jvr,nodexvr,maxcat,nbestcat);
							rightimp=0;
							for(n=1;n<=nout;n++){
								//the nth out-of-bag case is the nnth original case
								nn=joob[n];
								if(impn==1){
									if(jvr[n]==cl[nn]){
										qimpm[nn][mr]=qimpm[nn][mr]+
				     						tnodewt[nodexvr[n]]/nout;
									}
								}
								if(jvr[n]==cl[nn]) rightimp=rightimp+
				     						tnodewt[nodexvr[n]];
							}
							avimp[mr]=avimp[mr]+(right-rightimp)/nout;
							//sqsd(mr)=sqsd(mr)+((right-rightimp)**2)/(nout**2)
							sqsd[mr]=sqsd[mr]+((right-rightimp)*(right-rightimp))/(nout*nout);
						}else{
							for(n=1;n<=nout;n++){
								//the nth out-of-bag case is 
								//the nnth original case
								nn=joob[n];
								if(impn==1){
									if(jtr[nn]==cl[nn]){
										qimpm[nn][mr]=qimpm[nn][mr]+tnodewt[nodextr[nn]]/nout;
									}
								}
							}
						}
					} // !k
	}
	
	// fun19: finishimp
	public static void finishimp(int mdim,float[] sqsd,float[] avimp,float[] signif,float[] zscore,
		     	int jbt,int mdimt,int[] msm){		
			//real sqsd(mdim),avimp(mdim),zscore(mdim),signif(mdim)
		   	//integer msm(mdim)
			//integer mdim,mdimt,k,m1,jbt
			float v,av,se;
			float erfcc;
			int m1 = 0;
			for(int k=1;k<=mdimt;k++){
				m1=msm[k];
				avimp[m1]=avimp[m1]/jbt;
				av=avimp[m1];
				se=(sqsd[m1]/jbt)-av*av;
				se=1.0f*(float)Math.sqrt(se/jbt);
				if(se>0.0){
					zscore[m1]=avimp[m1]/se;
					v=zscore[m1];
					signif[m1]=erfcc(v);
				}else{
					zscore[m1]=-5;
					signif[m1]=1;
				}
			}
	}
	
	// fun20: compinteract
	public static void compinteract(float[][] votes,float[][] effect,int[] msm,int mdim,int mdimt,
		     	int jbt,float[] g,int[] iv,int[][] irnk,float[][] hist,float[][] teffect){		
			// real votes(mdim,jbt),effect(mdim,mdim),g(mdim),
		    //	hist(0:mdim,mdim),teffect(mdim,mdim)
			// integer msm(mdimt),iv(mdim),irnk(mdim,jbt)
			int  jb,i,nt,irk,j,ii;
		    int  jj,ij,mmin,m,k;
		    mmin = -1;
			float gmin,rcor;			
AAjbt:		for(jb=1;jb<=jbt;jb++){
				nt=0;
				zerv(iv,mdim);
				for(i=1;i<=mdimt;i++){
					m=msm[i];
					g[m]=votes[m][jb];
					if(Math.abs(g[m])<8.232D-11){
						irnk[m][jb]=0;
						iv[m]=1;
						nt=nt+1;
					}
				}
				irk=0;
				for(j=1;j<=8000;j++){
					gmin=10000;
					for(i=1;i<=mdimt;i++){
						m=msm[i];
						if(iv[m]==0&&g[m]<gmin){
							gmin=g[m];
							mmin=m;
						}
					}
					iv[mmin]=1;
					irk=irk+1;
					irnk[mmin][jb]=irk;
					nt=nt+1;
					if(nt>=mdimt) continue AAjbt; //goto 79
				} // !j
		//79		continue
			} // !jb
			
			for(j=0;j<=mdimt;j++){
				for(i=1;i<=mdimt;i++){
					m=msm[i];
					hist[j][m]=0;
				}
			}
			
			for(i=1;i<=mdimt;i++){
				m=msm[i];
				for(jb=1;jb<=jbt;jb++){
					hist[irnk[m][jb]][m]=hist[irnk[m][jb]][m]+1;
				}
				for(j=0;j<=mdimt;j++){
					hist[j][m]=hist[j][m]/jbt;
				}
			} // !m
			
			zermr(effect,mdim,mdim);			
			for(i=1;i<=mdimt;i++){
				for(j=1;j<=mdimt;j++){
					m=msm[i];
					k=msm[j];
					for(jb=1;jb<=jbt;jb++){
						effect[m][k]=effect[m][k]+Math.abs(irnk[m][jb]-irnk[k][jb]);
					}
					effect[m][k]=effect[m][k]/jbt;
				}
			}
			
			zermr(teffect,mdim,mdim);
			for(i=1;i<=mdimt;i++){
				for(j=1;j<=mdimt;j++){
					m=msm[i];
					k=msm[j];
					for(ii=0;ii<=mdimt;ii++){
						for(jj=0;jj<=mdimt;jj++){
							teffect[m][k]=teffect[m][k]+Math.abs(ii-jj)*hist[jj][m]*hist[ii][k];
						}
					}
					rcor=0;
					for(ij=1;ij<=mdimt;ij++){
						rcor=rcor+hist[ij][m]*hist[ij][k];
					}
					teffect[m][k]=teffect[m][k]/(1-rcor);
				}
			}

			for(i=1;i<=mdimt;i++){
				for(j=1;j<=mdimt;j++){
					m=msm[i];
					k=msm[j];
					effect[m][k]=100*(effect[m][k]-teffect[m][k]);
				}
			}
	}
	
	// fun21: compprot
	public static void compprot(int[][] loz,int nrnn,int ns,int mdim,int[] its,
		int[] jest,float[] wc,int nclass,float[][] x,int mdimt,int[] msm,float[] temp,int[] cat,int maxcat,
		int[] jpur,int[] inear,int nprot,float[][][] protlow,float[][][] prothigh,float[][][] prot,float[][][][] protfreq,
		float[][][] protvlow,float[][][] protvhigh,float[][][] protv,float[][] popclass,int[] npend,float[] freq,float[] v5,float[] v95){		
				// integer nrnn,ns,mdim,nclass,mdimt,nprot,maxcat
				// integer loz(ns,nrnn),jest(ns),msm(mdim),its(ns),
			    // 	jpur(nrnn),inear(nrnn),npend(nclass),cat(mdim)
				// real wc(ns),prot(mdim,nprot,nclass),
			    // 	protlow(mdim,nprot,nclass),prothigh(mdim,nprot,nclass),
			    //	protfreq(mdim,nprot,nclass,maxcat),
			    // 	protvlow(mdim,nprot,nclass),protvhigh(mdim,nprot,nclass),
			    //	x(mdim,ns),temp(nrnn),protv(mdim,nprot,nclass),
			    // 	popclass(nprot,nclass),freq(maxcat),v5(mdim),v95(mdim)
				int ii,i,k,n,jp,npu,mm,m,nclose,jj,ll,jmax,nn;
				float fmax,dt;
	AAjp:		for(jp=1;jp<=nclass;jp++){
					zerv(its,ns);
					//we try to find nprot prototypes for this class:
					npend[jp]=nprot;
					for(i=1;i<=nprot;i++){
						zervr(wc,ns);
						for(n=1;n<=ns;n++){
							if(its[n]==0){
								//	wc(n) is the number of unseen neighbors 
								//	of case n that are predicted to be 
								//	in class jp
								//	loz(n,1),...,loz(n,nrnn) point to 
								//	the nrnn nearest neighbors of 
								//	case n 
								for(k=1;k<=nrnn;k++){
									nn=loz[n][k];
									if(its[nn]==0){
										ii=jest[nn];
										if(ii==jp) wc[n]=wc[n]+1;
									}
								}
							}
						}
						//find the unseen case with the largest number 
						//of unseen predicted-class-jp neighbors
						nclose=0;
						npu=0;
						for(n=1;n<=ns;n++){
							if(wc[n]>=nclose&&its[n]==0){
								npu=n;
								nclose=(int)wc[n];
							}
						}
						//if nclose=0,no case has any unseen predicted-class-jp neighbors
						//can't find another prototype for this class - reduce npend by 1 and
						//start finding prototypes for the next class
						if(nclose==0){
							npend[jp]=i-1;
							//goto 93
							continue AAjp;
						}
						//case npu has the largest number 
						//of unseen predicted-class-jp neighbors
						//put these neighbors in a list of length nclose
						ii=0;	
						for(k=1;k<=nrnn;k++){
							nn=loz[npu][k];
							if(its[nn]==0&&jest[nn]==jp){
								ii=ii+1;
								inear[ii]=nn;	
							}
						}
						//popclass is a measure of the size of the cluster around 
						//this prototype
						popclass[i][jp]=nclose;
						for(mm=1;mm<=mdimt;mm++){
							//m is the index of the mmth variable
							m=msm[mm];	
							if(cat[m]==1){
								dt=v95[m]-v5[m];
								for(ii=1;ii<=nclose;ii++){
									//put the value of the mmth variable into the list
									temp[ii]=x[m][inear[ii]];
								} // !ii
								//sort the list
								Y_quicksort(temp,jpur,1,nclose,nclose);
								ii=nclose/4;
								if(ii==0) ii=1;
								//find the 25th percentile
								protvlow[m][i][jp]=temp[ii];
								protlow[m][i][jp]=(temp[ii]-v5[m])/dt;
								ii=nclose/4;
								ii=(3*nclose)/4;
								if(ii==0) ii=1;
								//find the 75th percentile
								protvhigh[m][i][jp]=temp[ii];
								prothigh[m][i][jp]=(temp[ii]-v5[m])/dt;
								ii=nclose/2;
								if(ii==0) ii=1;
								//find the median
								protv[m][i][jp]=temp[ii];
								prot[m][i][jp]=(temp[ii]-v5[m])/dt;
							}	
							if(cat[m]>=2){
								//for categorical variables,choose the most frequent class
								zervr(freq,maxcat);	
								for(k=1;k<=nclose;k++){
									jj=nint(x[m][loz[npu][k]]);	
									freq[jj]=freq[jj]+1;
								}
								jmax=1;
								fmax=freq[1];
								for(ll=2;ll<=cat[m];ll++){
									if(freq[ll]>fmax){
										jmax=ll;
										fmax=freq[ll];
									}
								}
								protv[m][i][jp]=jmax;
								protvlow[m][i][jp]=jmax;
								protvhigh[m][i][jp]=jmax;
								for(ll=1;ll<=cat[m];ll++){
									protfreq[m][i][jp][ll]=freq[ll];
								}
							}
						} // !m
						//record that npu and it's neighbors have been 'seen'
						its[npu]=1;
						for(k=1;k<=nclose;k++){
							nn=loz[npu][k];
							its[nn]=1;
						}
					} // !nprot
			//93		continue
				} // !jp	
	}
	
	// fun22: preprox
	public static void preprox(int near,int nrnodes,int jbt,int[] nodestatus,int[][] ncount,int jb,
		     	int[] nod,int[] nodextr,int[][] nodexb,int[] jin,int[][] jinb,int[] ncn,int[][] ndbegin,int[] kpop,
		     	float[][] rinpop,int[][] npcase,float[] rpop,int nsample){
		     	// integer near,nrnodes,jbt,jb,nsample,
		     	// nodestatus(nrnodes),ncount(near,jbt),nod(nrnodes),
		     	// nodextr(near),nodexb(near,jbt),jin(nsample),jinb(near,jbt),
		     	// ncn(near),kpop(near),npcase(near,jbt),ndbegin(near,jbt)
		     	// rpop(near),rinpop(near,jbt)
				int n, ntt, k, nterm, kn;			
				for(n=1;n<=near;n++){
					ncount[n][jb]=0;
					ndbegin[n][jb]=0;
				}				
				ntt=0;				
				for(k=1;k<=nrnodes;k++){
					if(nodestatus[k]==-1){
						ntt=ntt+1;
						nod[k]=ntt;
					}
				}			
				nterm=ntt;			
				for(n=1;n<=near;n++){
					rinpop[n][jb]=rpop[nodextr[n]];
					nodexb[n][jb]=nod[nodextr[n]];
					jinb[n][jb]=jin[n];
					k=nodexb[n][jb];
					ncount[k][jb]=ncount[k][jb]+1;
					ncn[n]=ncount[k][jb];
				}
				ndbegin[1][jb]=1;
				for(k=2;k<=(nterm+1);k++){
					ndbegin[k][jb]=ndbegin[k-1][jb]+ncount[k-1][jb];
				}
				for(n=1;n<=near;n++){
					kn=ndbegin[nodexb[n][jb]][jb]+ncn[n]-1;
					npcase[kn][jb]=n;
				}
	}
	
	// fun23: comprox
	public static void comprox(float[][] prox,int[][] nodexb,int[][] jinb,int[][] ndbegin,
		     	int[][] npcase,float[] ppr,float[][] rinpop,int near,int jbt,int noutlier,float[] outtr,int[] cl,
		     	int[][] loz,int nrnn,float[] wtx,int nsample,int[] iwork,int[] ibest){		
			//double precision prox(near,nrnn),ppr(near)
			//real outtr(near),rinpop(near,jbt),wtx(nsample)
			//integer nodexb(near,jbt),jinb(near,jbt),
		    //	ndbegin(near,jbt),npcase(near,jbt),cl(nsample),
		    //	loz(near,nrnn),iwork(near),ibest(nrnn),
		    //	near,jbt,noutlier,nrnn,nsample
			int n,jb,k,j,kk;
			float rsq;		
			for(n=1;n<=near;n++){
				zervd(ppr,near);
				for(jb=1;jb<=jbt;jb++){
					k=nodexb[n][jb];
					if(jinb[n][jb]>0){
						for(j=ndbegin[k][jb];j<=(ndbegin[k+1][jb]-1);j++){
							kk=npcase[j][jb];
							if(jinb[kk][jb]==0){
								ppr[kk]=ppr[kk]+(wtx[n]/rinpop[n][jb]);
							}
						}
					}
					if(jinb[n][jb]==0){
						for(j=ndbegin[k][jb];j<=(ndbegin[k+1][jb]-1);j++){
							kk=npcase[j][jb];
							if(jinb[kk][jb]>0){
								ppr[kk]=ppr[kk]+(wtx[kk]/rinpop[kk][jb]);
							}
						}
					}
				} // !jbt
				if(noutlier==1){
					rsq=0;
					for(k=1;k<=near;k++){
						if(ppr[k]>0&&cl[k]==cl[n]) rsq=rsq+ppr[k]*ppr[k];
					}
					if(rsq==0) rsq=1;
					outtr[n]=near/rsq;
				}
			
				if(nrnn==near){
					for(k=1;k<=near;k++){
						prox[n][k]=ppr[k];
						loz[n][k]=k;
					}
				}else{
					biggest(ppr,near,nrnn,ibest,iwork);
					for(k=1;k<=nrnn;k++){
						prox[n][k]=ppr[ibest[k]];
						loz[n][k]=ibest[k];
					}
				}
			} // !n
	}	
	
	// fun24: biggest
	public static void biggest(float[] x,int n,int nrnn,int[] ibest,int[] iwork){
		//double precision x(n)
		//integer n,nrnn,ibest(nrnn),iwork(n)
		int i,j,ihalfn,jsave;	
		// finds the nrnn largest values in the vector x and 
		// returns their positions in the vector ibest(1),...,ibest(nrnn):
		//	x(ibest(1)) is the largest
		//	...
		//	x(ibest(nrnn)) is the nrnn-th-largest
		// the vector x is not disturbed 
		// the vector iwork is used as workspace	
		ihalfn=(int)(n/2);
		for(i=1;i<=n;i++){
			iwork[i]=i;
		}
		for(j=1;j<=ihalfn;j++){
			i=ihalfn + 1 - j;
			sift(x,iwork,n,n,i);
		}
		for(j=1;j<=(nrnn-1);j++){
			i=n-j+1;
			ibest[j]=iwork[1];
			jsave=iwork[i];
			iwork[i]=iwork[1];
			iwork[1]=jsave;
			sift(x,iwork,n,n-j,1);
		}
		ibest[nrnn]=iwork[1];
	}
	
	// fun25: sift
	public static void sift(float[] x,int[] iwork,int n,int m,int i){
		//double precision x(n)
		//integer iwork(m),n,m,i
		float xsave;
		int j,k,jsave,ksave;	
		//used by subroutine biggest,to bring the largest element to the 
		//top of the heap	
		xsave=x[iwork[i]];
		ksave=iwork[i];
		jsave=i;
		j=i + i;
AA10:	for(k=1;k<=m;k++){
			if(j>m) break AA10;
			if(j<m){
				if(x[iwork[j]]<x[iwork[j+1]]) j =j+1;
			}
			if(xsave>=x[iwork[j]]) break AA10; //goto 10
			iwork[jsave]=iwork[j];
			jsave=j;
			j=j + j;
		}
	//10	continue
		iwork[jsave]=ksave;
		//return
	}
	
	// fun26: locateout	
	public static void locateout(int[] cl,float[] tout,float[] outtr,int[] ncp,int[] isort,float[] devout,
		     	float near,int nsample,int nclass,float[] rmedout){		
			//real outtr(near),tout(near),devout(nclass),rmedout(nclass)
		    //integer cl(nsample),isort(nsample),ncp(near),near,nsample,nclass
			float rmed, dev;
			int jp,nt,n,i;		     
			for(jp=1;jp<=nclass;jp++){
				nt=0;
				for(n=1;n<=near;n++){
					if(cl[n]==jp){
						nt=nt+1;
						tout[nt]=outtr[n];
						ncp[nt]=n;
					}
				}
				Y_quicksort(tout,isort,1,nt,nsample);
				rmed=tout[(1+nt)/2];
				dev=0;
				for(i=1;i<=nt;i++){
					dev=dev+amin1((float)Math.abs(tout[i]-rmed),5*rmed);
				}
				dev=dev/nt;
				devout[jp]=dev;
				rmedout[jp]=rmed;
				for(i=1;i<=nt;i++){
					outtr[ncp[i]]=amin1((outtr[ncp[i]]-rmed)/dev,20.0f);
				}
			} // !jp
	}
	
	// fun27: myscale
	public static void myscale(int[][] loz,float[][] prox,float[][] xsc,float[] y,float[] u,int near,int nscale,float[] red,int nrnn,
		     	float[] ee,float[][] ev,float[] dl){		
			// double precision prox(near,nrnn),y(near),u(near),dl(nscale),
		    //  xsc(near,nscale),red(near),ee(near),ev(near,nscale)
			float[] bl = new float[10+1];
			//integer loz(near,nrnn)		
			int j,i,it,n,jit,k;
			float dotd,y2,sred,eu,ru,ra,ynorm,sa;
			for(j=1;j<=near;j++){
				ee[j]=1f;
			}		
			for(j=1;j<=near;j++){
				red[j]=0;
				for(i=1;i<=nrnn;i++){
					red[j]=red[j]+prox[j][i];
				}
				red[j]=red[j]/near;
			}
			sred=dotd(ee,red,near);
			sred=sred/near;
			for(it=1;it<=nscale;it++){
				for(n=1;n<=near;n++){
					//if(mod(n,2).eq.0) then
					if(n%2==0){
						y[n]=1;
					}else{
						y[n]=-1;
					}
				}
   AA101:		for(jit=1;jit<=1000;jit++){
					y2=dotd(y,y,near);
					y2=(float)Math.sqrt(y2);
					for(n=1;n<=near;n++){
						u[n]=y[n]/y2;
					}
					for(n=1;n<=near;n++){
						y[n]=0;
						for(k=1;k<=nrnn;k++){
							y[n]=y[n]+prox[n][k]*u[loz[n][k]];
						}
					}
					eu=dotd(ee,u,near);
					ru=dotd(red,u,near);
					for(n=1;n<=near;n++){
						y[n]=y[n]-(red[n]-sred)*eu-ru;
						y[n]=0.5f*y[n];
					}
					if(it>1){
						for(j=1;j<=(it-1);j++){
							bl[j]=0;
							for(n=1;n<=near;n++){
								bl[j]=bl[j]+ev[n][j]*u[n];
							}
							for(n=1;n<=near;n++){
								y[n]=y[n]-bl[j]*dl[j]*ev[n][j];
							}
						}
					}
					ra=dotd(y,u,near);
					ynorm=0;
					for(n=1;n<=near;n++){
						//ynorm=ynorm+(y[n]-ra*u[n])**2
						ynorm=ynorm+(y[n]-ra*u[n])*(y[n]-ra*u[n]);
					}
					//sa=dabs(ra)
					sa=Math.abs(ra);
					if(ynorm<(sa*1.0e-7)){
						for(n=1;n<=near;n++){
							xsc[n][it]=(float)(Math.sqrt(sa))*u[n];
							ev[n][it]=u[n];
						}
						dl[it]=ra;
						break AA101;
					}
				}
		//101		continue
			} //!nn
	}
	
	// fun28: xfill
	public static void xfill(float[][] x,int nsample,int mdim,float[] fill,float code){
		//	input:
		//		real code,fill(mdim)
		//		integer mdim,nsample
		//	output: real x(mdim,nsample)
		int n,m;
		for(n=1;n<=nsample;n++){
			for(m=1;m<=mdim;m++){
				if(Math.abs(x[m][n]-code)<8.232D-11){ 
	     				x[m][n]=fill[m];
				}
			}
		}
	}
	
	// fun29: roughfix
	public static void roughfix(float[][] x,float[] v,int[] ncase,int mdim,int nsample,int[] cat,float code,
		     	int[] nrcat,int maxcat,float[] fill){		
			//real x(mdim,nsample),v(nsample),fill(mdim),code
			//integer ncase(nsample),cat(mdim),nrcat(maxcat)
			//integer mdim,nsample,maxcat
			int m,n,nt,j,jmax,lcat,nmax;
			float rmed;
		
			for(m=1;m<=mdim;m++){
				if(cat[m]==1){
					nt=0;
					for(n=1;n<=nsample;n++){
						if(Math.abs(x[m][n]-code)>=8.232D-11){
							nt=nt+1;
							v[nt]=x[m][n];
						}
					}
			   		Y_quicksort (v,ncase,1,nt,nsample);
					if(nt>0){
						rmed=v[(nt+1)/2];
					}else{
						rmed=0;
					}
					fill[m]=rmed;
				}else{
					lcat=cat[m];
					zerv(nrcat,maxcat);
					for(n=1;n<=nsample;n++){
						if(Math.abs(x[m][n]-code)>=8.232D-11){
							j=nint(x[m][n]);
							nrcat[j]=nrcat[j]+1;
						}
					}
					nmax=0;
					jmax=1;
					for(j=1;j<=lcat;j++){
						if(nrcat[j]>nmax){
							nmax=nrcat[j];
							jmax=j;
						}
					}
					fill[m]=(float)(jmax);
				}
			} // !m
			for(n=1;n<=nsample;n++){
				for(m=1;m<=mdim;m++){
					if(Math.abs(x[m][n]-code)<8.232D-11) x[m][n]=fill[m];
				}
			}
	}	
	
	// fun30: impute
	public static void impute(float[][] x,float[][] prox,int near,int mdim,
			     	int maxcat,float[] votecat,int[] cat,int nrnn,int[][] loz,int[][] missing){			
				//	real x(mdim,near),votecat(maxcat)
				//	double precision  prox(near,nrnn)
				//	integer near,mdim,maxcat,nrnn
				//	integer cat(mdim),loz(near,nrnn),
				//  missing(mdim,near)		     	
				int i,j,jmax,m,n,k;
				jmax=0;
				float sx,dt,rmax;
				for(m=1;m<=mdim;m++){
					if(cat[m]==1){
						for(n=1;n<=near;n++){
							if(missing[m][n]==1){
								sx=0;
								dt=0;
								for(k=1;k<=nrnn;k++){
									if(missing[m][loz[n][k]]!=1){
										sx=sx+1f*prox[n][k]*x[m][loz[n][k]];
										dt=dt+1f*prox[n][k];
									}
								}
								if(dt>0) x[m][n]=sx/dt;
							}
						}  //n
					}
				} //m
				for(m=1;m<=mdim;m++){
					if(cat[m]>1){
						for(n=1;n<=near;n++){
							if(missing[m][n]==1){
								zervr(votecat,maxcat);
								for(k=1;k<=nrnn;k++){
									if (missing[m][loz[n][k]]!=1){
										j=nint(x[m][loz[n][k]]);
										votecat[j]=votecat[j]+1f*prox[n][k];
									}
								}  //k
								rmax=-1;
								for(i=1;i<=cat[m];i++){
									if(votecat[i]>rmax){
										rmax=votecat[i];
										jmax=i;
									}
								}
								x[m][n]=1.0f*jmax;
							}
						} //n
					}
				} //m
	}	
	
	// fun31: checkin	
	public static void checkin(int labelts,int labeltr,int nclass,int lookcls,int jclasswt,
				     	int mselect,int mdim2nd,int mdim,int imp,int impn,int interact,int nprox,int nrnn,
				     	int nsample,int noutlier,int nscale,int nprot,int missfill,int iviz,int isaverf,
				     	int irunrf,int isavepar,int ireadpar,int isavefill,int ireadfill,int isaveprox,
				     	int ireadprox,int isumout,int idataout,int impfastout,int impout,int interout,
				     	int iprotout,int iproxout,int iscaleout,int ioutlierout,int[] cat,int maxcat,int[] cl){			
						//	integer labelts,labeltr,nclass,lookcls,jclasswt,mselect,
						//	&	mdim2nd,mdim,imp,impn,interact,nprox,nrnn,nsample,noutlier,
						//	&	nscale,nprot,missfill,iviz,isaverf,irunrf,isavepar,ireadpar,
						//	&	isavefill,ireadfill,isaveprox,ireadprox,isumout,idataout,
						//	&	impfastout,impout,interout,iprotout,iproxout,iscaleout,
						//	&	ioutlierout,cat(mdim),maxcat,cl(nsample)
					int n,m;
					if(labelts!=0&&labelts!=1){
						System.out.println("error labelts: "+labelts);
						System.exit(1);
					}
					if(labeltr!=0&&labeltr!=1){
						System.out.println("error labeltr: "+labeltr);
						System.exit(1);
					}
					if(labeltr==0&&nclass!=2){
						System.out.println("error,nclass should be 2 if labeltr=0");
						System.exit(1);
					}
					if(lookcls!=0&&lookcls!=1){
						System.out.println("error lookcls: " + lookcls);
						System.exit(1);
					}
					if(jclasswt!=0&&jclasswt!=1){
						System.out.println("error jclasswt: " + jclasswt);
						System.exit(1);
					}
					if(mselect!=0&&mselect!=1){
						System.out.println("error mselect: " + mselect);
						System.exit(1);
					}
					if(mdim2nd<0||mdim2nd>mdim){
						System.out.println("error mdim2nd: " + mdim2nd);
						System.exit(1);
					}
					if(imp!=0&&imp!=1){
						System.out.println("error imp: " + imp);
						System.exit(1);
					}
					if(impn!=0&&impn!=1){
						System.out.println("error impn: " + impn);
						System.exit(1);
					}
					if(interact!=0&&interact!=1){
						System.out.println("error interact: " + interact);
						System.exit(1);
					}
					if(nprox<0||nprox>1){
						System.out.println("error nprox: " + nprox);
						System.exit(1);
					}
					if(nrnn<0||nrnn>nsample){
						System.out.println("error - nrnn: " + nrnn);
						System.exit(1);
					}
					if(noutlier<0||noutlier>2){
						System.out.println("error - noutlier: " + noutlier);
						System.exit(1);
					}
					if(nscale<0||nscale>mdim){
						System.out.println("error - nscale:  " + nscale);
						System.exit(1);
					}
					if(nprot<0||nprot>nsample){
						System.out.println("error - nprot:" + nprot);
						System.exit(1);
					}
					if(missfill<0||missfill>2) {
						System.out.println("error missfill:  " +missfill);
						System.exit(1);
					}
					if(iviz!=0&&iviz!=1) {
						System.out.println("error iviz:  " +iviz);
						System.exit(1);
					}
					if(iviz==1&&impn!=1) {
						System.out.println("error iviz=1 and impn=:  " +impn);
						System.exit(1);
					}
					if(iviz==1&&imp!=1) {
						System.out.println("error iviz=1 and imp=:  " +imp);
						System.exit(1);
					}
					if(iviz==1&&nprox!=1) {
						System.out.println("error iviz=1 and nprox=:  " +nprox);
						System.exit(1);
					}
					if(iviz==1&&nscale!=3) {
						System.out.println("error iviz=1 and nscale=:  " +nscale);
						System.exit(1);
					}
					if(isaverf!=0&&isaverf!=1) {
						System.out.println("error isaverf:  " +isaverf);
						System.exit(1);
					}
					if(isaverf==1&&missfill==2) {
						System.out.println("error - only rough fix can be saved");
						System.exit(1);
					}
					if(irunrf!=0&&irunrf!=1) {
						System.out.println("error irunrf:  " +irunrf);
						System.exit(1);
					}
					if(isavepar!=0&&isavepar!=1) {
						System.out.println("error isavepar:  " +isavepar);
						System.exit(1);
					}
					if(ireadpar!=0&&ireadpar!=1) {
						System.out.println("error ireadpar:  " +ireadpar);
						System.exit(1);
					}
					if(isavefill!=0&&isavefill!=1) {
						System.out.println("error isavefill:  " +isavefill);
						System.exit(1);
					}
					if(ireadfill!=0&&ireadfill!=1) {
						System.out.println("error ireadfill:  " +ireadfill);
						System.exit(1);
					}
					if(isaveprox!=0&&isaveprox!=1) {
						System.out.println("error isaveprox:  " +isaveprox);
						System.exit(1);
					}
					if(ireadprox!=0&&ireadprox!=1) {
						System.out.println("error ireadprox:  " +ireadprox);
						System.exit(1);
					}
					if(isumout<0||isumout>1) {
						System.out.println("error isumout:  " +isumout);
						System.exit(1);
					}
					if(idataout<0||idataout>2) {
						System.out.println("error idataout:  " +idataout);
						System.exit(1);
					}
					if(impfastout<0||impfastout>1) {
						System.out.println("error impfastout:  " +impfastout);
						System.exit(1);
					}
					if(impout<0||impout>2) {
						System.out.println("error impout:  " +impout);
						System.exit(1);
					}
					if(interout<0||interout>2) {
						System.out.println("error interout:  " +interout);
						System.exit(1);
					}
					if(iprotout<0||iprotout>2) {
						System.out.println("error iprotout:  " +iprotout);
						System.exit(1);
					}
					if(iproxout<0||iproxout>2) {
						System.out.println("error iproxout:  " +iproxout);
						System.exit(1);
					}
					if(iscaleout<0||iscaleout>1) {
						System.out.println("error iscaleout:  " + iscaleout);
						System.exit(1);
					}
					if(ioutlierout<0||ioutlierout>2) {
						System.out.println("error ioutlierout:  " + ioutlierout);
						System.exit(1);
					}
					if(noutlier>0&&nprox==0) {
						System.out.println("error - noutlier>0 and nprox=0");
						System.exit(1);
					}
					if(nscale>0&&nprox==0) {
						System.out.println("error - nscale>0 and nprox=0");
						System.exit(1);
					}
					if(nprot>0&&nprox==0) {
						System.out.println("error - nprot>0 and nprox=0");
						System.exit(1);
					}
					if(mdim2nd>0&&imp==0) {
						System.out.println("error - mdim2nd>0 and imp=0");
						System.exit(1);
					}
					if(impn>0&&imp==0) {
						System.out.println("error - impn>0 and imp=0");
						System.exit(1);
					}
					for(m=1;m<=mdim;m++){
						if(cat[m]>maxcat) {
							System.out.println("error in cat:  " +m + " " + cat[m]);
							System.exit(1);
						}
					}
					if(labeltr==1) {
						for(n=1;n<=nsample;n++){
							if(cl[n]<1||cl[n]>nclass) { 
								System.out.println("error in class label:  " +n + " " + cl[n]);
								System.exit(1);
							}
						}
					}
	}		
	
	// fun32: quick fun
//	public static void MY_quicksort(float[] v,int[] iperm,int ii,int jj,int kk){
//		//	puts into iperm the permutation vector which sorts v into
//		//	increasing order. only elementest from ii to jj are considered.
//		//	array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
//		//	this is a modification of acm algorithm #347 by r. c. singleton,
//		//	which is a modified hoare quicksort.		
//		// start from position 1		
//		int len = v.length;
//		int lenR = iperm.length;		
//		len = ii-jj+1;
//		lenR= ii-jj+1;			
//		float[] new_v = new float[len];
//		int[] new_p = new int[lenR];
//		int i,j;
//		for(i=ii;i<=jj;i++){
//			new_v[i]=v[i];
//		}
//		for(i=ii;i<=jj;i++){
//			new_p[i]=iperm[i];
//		}
//		new_v[0]=-1000;
//		Arrays.sort(new_v); // sort the vector new_v
//		for(i=0;i<len;i++){          // old position
//			for(j=0;j<len;j++){      // new position
//				if(Math.abs(v[i]-new_v[j])<0.001){
//					new_p[j] = iperm[i];
//				}
//			}
//		}		
//		v[0]=-1000;
//		Arrays.sort(v);
//		for(i=ii;i<=jj;i++){
//			iperm[i]=new_p[i];
//		}
//		new_p=null;
//		new_v=null;		
//	}
	// end of qucick sort
	
	// fun32: quick fun
	public static void Y_quicksort(float[] v,int[] iperm,int ii,int jj,int kk){
			//	puts into iperm the permutation vector which sorts v into
			//	increasing order. only elementest from ii to jj are considered.
			//	array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
			//	this is a modification of acm algorithm #347 by r. c. singleton,
			//	which is a modified hoare quicksort.		
			// start from position 1				
			 int i,j;
			 float temp;
			 int len = iperm.length;
			 int lenv = v.length;
			 int[] old_iperm = new int[len+1];
			 float[] old_v = new float[lenv+1];
			 //initilize the values
			 for(i=0;i<len;i++){
				 old_v[i] = v[i];
			 }
			 for(i=0;i<lenv;i++){
				 old_iperm[i] = iperm[i];
			 }			 
			 //sorts v into increasing order
			 for(i=ii;i<=jj;i++){
				 for(j=i+1;j<=jj;j++){				 
					 if(v[i]>v[j]){
						 temp = v[i];
						 v[i] = v[j];
						 v[j] = temp;
					 }				 
				 }
			 }
			 for(i=1;i<lenv;i++){          // old position
					for(j=1;j<lenv;j++){      // new position
						if(Math.abs(v[i]-old_v[j])<0.001){
							iperm[i] = old_iperm[j];
						}
					}
			 }					
	}
	// end of qucick sort
	
	// fun32: quick fun
	/* R_qsort_I(v, index, 1, nsample)  This sorts the v(n) in ascending order. index(n) is the case
    number of that v(n) nth from the lowest (assume the original
    case numbers are 1,2,...).  */	
//	public static void quicksort(float[] v,int[] iperm,int ii,int jj,int kk){
//		//	puts into iperm the permutation vector which sorts v into
//		//	increasing order. only elementest from ii to jj are considered.
//		//	array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
//		//	this is a modification of acm algorithm #347 by r. c. singleton,
//		//	which is a modified hoare quicksort.
//		float vt,vtt;
//		//integer t,tt,iperm(kk),iu(32),il(32)
//		int t, tt;
//		iperm = new int[kk+1];
//		int m,i,j,k,ij,l;
//		m=1;
//		i=ii;
//		j=jj;	
//	}
	
	// fun33: unpack
	public static void unpack(int l,int npack,int[] icat){
		int j,n,k;
		if(l>32){
			System.out.println("error in unpack,l=" + l);
			System.exit(1);
		}
		for(j=1;j<=32;j++){
			icat[j]=0;
		}
		n=npack;
		icat[1]=n%2;
		for(k=2;k<=l;k++){
			n=(n-icat[k-1])/2;
			icat[k]=n%2;
		}
	}
	
	// fun34: zerv
	public static void zerv(int[] ix,int m1){
			int n;
			for(n=0;n<=m1;n++){
				ix[n]=0;
			}
	}
	
	 
    // fun35: zervr
    public static void zervr(float[] rx,int m1){
	   	 for(int n=1;n<=m1;n++){
	            rx[n]=0;
	   	 }
    }
	
	// fun36: zervd
	public static void zervd(float[] rx,int m1){
		int n;
		for(n=1;n<=m1;n++){
			rx[n]=0;
		}
	}
	
	// fun37: zerm
    public static void zerm(int[][] mx,int m1,int m2){
	     for(int j=1;j<=m2;j++){
	             for(int i=1;i<=m1;i++){
	                     mx[i][j]=0;
	             }
	     }
    }  
    
	// fun38: zermr
    public static void zermr(float[][] rx,int m1,int m2){
	    for(int j=1;j<=m2;j++){
	            for(int i=1;i<=m1;i++){
	                    rx[i][j]=0;
	            }
	    }
    }

	// fun39: zermd
	public static void zermd(float[][] rx,int m1,int m2){
		int i,j;
		for(j=1;j<=m2;j++){
			for(i=1;i<=m1;i++){
				rx[i][j]=0;
			}
		}
	}
   
 	// fun40: erfcc
	public static float erfcc(float x){
		float t,z, erfcc=-1;
		z=Math.abs(x)/1.41421356f;
		t=1f/(1f+0.5f*z);			
		erfcc=(float)(t*Math.exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+(t*
			     (0.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+(t*
			     (1.48851587+t*(-.82215223+t*.17087277)))				     
			    		 )		     
			    		 ))))				    		 
		))));					
		erfcc=erfcc/2;
		if (x<0) erfcc=2-erfcc;
		return erfcc;
	}
	
	// fun41: irbit	: this function only useful for categorical variable
	public static int irbit(int iseed){
		int irbit=0;
		int ib1, ib2, ib5, ib18, mask;
		ib1=1;
		ib2=2;
		ib5=16;
		ib18=131072;
		mask=ib1+ib2+ib5;
		if((iseed&ib18)!=0){
			//iseed=ior(ishft(ieor(iseed,mask),1),ib1)
			iseed=(((iseed^mask)<<1)|ib1);
			irbit=1;
		}else{
			//iseed=iand(ishft(iseed,1),not(ib1))
			iseed=((iseed<<1)&(~(ib1)));
			irbit=0;
		}
		return irbit;
	}	
	
	// fun42 randomu
	public static float randomu(){
		return (float)Math.random();
		//return (float)0.6f;
	}
       
	// fun43: rnorm
	public static float rnorm(int i){
		return (float)(Math.sqrt(-2*Math.log(randomu()))*Math.cos(6.283185*randomu()));
	}
     
	// fun44: dotd
 	public static float dotd(float[] u,float[] v,int ns){
 		//computes the inner product
 		float dotd=0;
 		for(int n=1;n<=ns;n++){
 				dotd=dotd+u[n]*v[n];
 		} 	
 		return dotd;
    }
	
 	// fun45: nint
    public static int nint(float a){
    	 if(a>0){
    		 return (int)(a+0.5);
    	 }else{
    		 return (int)(a-0.5);
    	 }
    }
 		
	// fun46: amin1	
	public static float amin1(float a, float b){
		if(a<b)
			return a;
		return b;
	}
}