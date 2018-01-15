package machinelearning;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
/*
 *      protein: /home/rxyan/bin/umdhmm-v1.02 ( delte *.o files before make).
 *      
 * 		Hidden Markov Model code was written by Renxiang Yan in the University of Michigan, Ann Anbor
 *      Date: 2011-10-14
 *      Date: 2011-10-27:change the variable names of methods
 *      
 * 		Initial framework finished Date: 2011-10-17
 * 
 * 		// Problem 1 Given the observation sequence O and model , calculate the P(O|model)
 *      // Forwardscale method
 *      
 * 		// Problem 2 Given the observation sequence O and model, identify the most possible states 
 * 		// Viterbi algorithm: Identify the most possible by recursion (dynamic programming)
 * 		
 * 		// Problem 3 How do we adjust the model parameter to maximize the P(O|model) 
 *      // BaumWelch method
 *      
 * 		References:
 *  		
 *      [1] Richard Durbin, Sean R. Eddy, Anders Krogh, Graeme Mitchison. Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids. Cambridge University Press, 1999
 *      [2] L. R. Rabiner, A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition, Proc. of the IEEE, Vol.77, No.2, pp.257--286, 1989.
 * 
 */
public class MHMM{
	
	public static void main(String[] args) throws Exception{		
//		if(args.length!=3)
//		{
//			System.out.println("Usage error");
//			System.out.println("Usage: testvit <model.hmm> <obs.seq> \n");
//			System.exit(0);
//		}
		
//		System.exit(1);
		
		HMM hmm = new HMM("E:\\umdhmm-v1.02\\000test.hmm");	
//		HMM hmm = new HMM(args[0]);
		
//		hmm.Output_HMMMode();		
		
		int[] obs = ReadSequence("E:\\umdhmm-v1.02\\test.seq");
//		int[] obs = ReadSequence(args[1]);
		
//		hmm.Output_HMMModel();
		
//		System.out.println(obs.length);
		
		int[] q = new int[obs.length];
		double[][] delta = new double[obs.length+1][hmm.getN()+1];
		int[][] psi = new int[obs.length+1][hmm.getN()+1];			
		double[][] alpha = new double[obs.length][hmm.getN()+1];
		double[] scale = new double[obs.length];		
		
		// Question 3		
//		gamma = dmatrix(1, T, 1, hmm.N);
		double[][] gamma = new double[obs.length][hmm.getN()+1];
		double[][] beta = new double[obs.length][hmm.getN()+1];		
		String[] pred = BaumWelch(hmm,obs.length-1,obs,alpha, beta,gamma).split("\\s+");				
//
//		
//		
		hmm.Output_HMMModel();			
		System.exit(0);		
		///// final log score divid the length of the sequence would be better
//		System.out.println("\nIteration=" + pred[0] + "initial Log" + pred[1] + " final Log" + pred[2]);
			
// Question 2
		double prob = Viterbi(hmm,obs.length-1,obs,delta,psi,q);		
		DecimalFormat df = new DecimalFormat("###.###");      // Correct the identity to 3 decimal places. 
		System.out.println("Viterbi  MLE log prob " + df.format(Math.log(prob)));		
		Output_status(q);
		System.out.println();
		prob=ViterbiLog(hmm,obs.length-1,obs,delta,psi,q);		
		Output_status(q);
		System.out.println();
		System.out.println("Viterbi  MLE log prob " + df.format(prob));
		
		// Question 1
		//////////  test for here		
//		double test = Forward(hmm, obs.length-1, obs,alpha);		
//		System.out.println(Math.log(test));
//		double test2 = ForwardWithScale(hmm, obs.length-1,obs,alpha,scale);		
//		System.out.println(test2);
		
		//////////  End of question 1	
//		for(int i=1;i<obs.length-1;i++)
//		{
//			for(int j=1;j<=hmm.getN();j++)
//			{
//				System.out.print(alpha[i][j] + " ");
//			}
//			System.out.println();
//		}
//		System.out.println("\nT=" + (obs.length-1));
//		System.out.println("\n%" + test);

		// testfor subrout
	}
	
	// start of random sequence
	
	
	
	// end of random sequence
	
	public static int num(char a){
		String order = "GAVLISTCMPDNEQKRHFYW";
		int i;
		for(i=0;i<order.length();i++){
			if(a==order.charAt(i)){
				return (i+1);
			}
		}
		return 1;
	}
    
	public static String read_fasta(String fasta) throws Exception {
		BufferedReader br = null;
		String line = "";
		String x = "";
		br = new BufferedReader(new FileReader(fasta));
		while ((line = br.readLine()) != null) {
			if (line.startsWith(">")) {

			} else {
				x = x + line;
			}
		}
		br.close();
		return ("*" + x);
	}
	
	public static int[] ReadSequence_fasta(String seq) throws Exception{
		int[] Obs = null;
		String ss = read_fasta(seq);
		Obs = new int[ss.length()];
		for(int i=1;i<ss.length();i++){
			Obs[i] = num(ss.charAt(i));
		}
		return Obs;
	}
	
	public static int[] ReadSequence(String seq) throws Exception{
		BufferedReader br;
		String line="";
		int[] Obs = null;
		br = new BufferedReader(new FileReader(seq));
			  while ((line = br.readLine()) != null) {
					if(!line.startsWith("T")){						
						String[] wds = line.split("\\s+");
						Obs = new int[wds.length+1];
						for(int k=1;k<=wds.length;k++){
							Obs[k] = Integer.valueOf(wds[k-1]).intValue();
							//System.out.println(Obs[k]);
						}
					}
			 }
		return Obs;
	}

	public  static String Output_status(int[] q){
		String ss = "";
		for(int i=1;i<q.length;i++){
			//System.out.print(q[i] + "\t");
			if(q[i]==1){
				//System.out.print('H');
				ss = ss + "H";
			}else if(q[i]==2){
				//System.out.print('E');
				ss = ss + "E";
			}else{
				//System.out.print('C');
				ss = ss + "C";
			}
		}		
		return ss;
	}
	
	public static double Viterbi(HMM hmm,int T,int[] O,double[][] delta,int[][] psi,int[] q){
		int i, j;	/* state indices */
		int t;		/* time index */	
		int	maxvalind;
		double	maxval, val;
		
		/* 1. Initialization  */
		
		for (i = 1; i <= hmm.getN(); i++) 
		{
			//delta[1][i] = phmm->pi[i] * (phmm->B[i][O[1]]);
			delta[1][i] = hmm.getPI()[i]*(hmm.getB()[i][O[1]]);			
			psi[1][i] = 0;
		}	
		
		/* 2. Recursion */		
		for (t = 2; t <= T; t++) 
		{
			//for (j = 1; j <= phmm->N; j++) {
			for (j = 1; j <= hmm.getN(); j++) 
			{
				maxval = 0.0;
				maxvalind = 1;	
				//for (i = 1; i <= phmm->N; i++) {
				for (i = 1; i <= hmm.getN(); i++) 
				{
					//val = delta[t-1][i]*(phmm->A[i][j]);
					val = delta[t-1][i]*(hmm.getA()[i][j]);
					if (val > maxval) 
					{
						maxval = val;	
						maxvalind = i;	
					}
				}				
				//delta[t][j] = maxval*(phmm->B[j][O[t]]);				
				int kkk = O[t];
				double temp = maxval*(hmm.getB()[j][kkk]);
				delta[t][j] = temp;
				psi[t][j] = maxvalind;
			}
		}
			
		/* 3. Termination */

		double pprob = 0.0;
		q[T] = 1;
		
		//System.out.println("\n\nDelta output:\n");
		//// for (i = 1; i <= phmm->N; i++) {
		for (i = 1; i <= hmm.getN(); i++) {
				//System.out.println(delta[T][i]);
	            if (delta[T][i] > pprob) {
				pprob = delta[T][i];	
				q[T] = i;
			}
		}
		//System.out.println("End\n");
		
		
		/* 4. Path (state sequence) backtracking */

		for (t = T - 1; t >= 1; t--)
			q[t] = psi[t+1][q[t+1]];
		
		return pprob;
	}
	
	
	public static double ViterbiLog(HMM hmm,int T,int[] O,double[][] delta,int[][] psi,int[] q)
	{
		 int     i, j;   /* state indices */
	     int     t;      /* time index */
	     
	     int     maxvalind;
	     double  maxval, val;
		 double[][]  biot = null;
		 
		/* 0. Preprocessing */

		// for (i = 1; i <= phmm->N; i++) 
		   for(i=1;i<=hmm.getN();i++)
				//phmm->pi[i] = log(phmm->pi[i]);
			       hmm.setPI_i(i, Math.log(hmm.getPI_i(i)));
			      
			
//		for (i = 1; i <= phmm->N; i++) 
//			for (j = 1; j <= phmm->N; j++) {
		   
//			phmm->A[i][j] = log(phmm->A[i][j]);
//			}
		   
		   for (i = 1; i <= hmm.getN(); i++)   
		   {
			   for (j = 1; j <= hmm.getN(); j++) 
			   {
				   hmm.setA_i_j(i, j, Math.log(hmm.getA_i_j(i, j)));
			   }
		   }
		   
		   // biot = dmatrix(1, phmm->N, 1, T);		   
		   biot = new double[hmm.getN()+1][T+1];
		   
			for (i = 1; i <= hmm.getN(); i++) 
			{
				for (t = 1; t <= T; t++) 
				{					
					biot[i][t] = Math.log(hmm.getB()[i][O[t]]);
				}
			}
			
			/* 1. Initialization  */
			 
	        for (i = 1; i <= hmm.getN(); i++) {	               
	        		delta[1][i] = hmm.getPI()[i] + biot[i][1];
	                psi[1][i] = 0;
	        }
	        
	        /* 2. Recursion */
	        
	        for (t = 2; t <= T; t++) 
	        {
	              //  for (j = 1; j <= phmm->N; j++) {
	        		for (j = 1; j <= hmm.getN(); j++) 
	        		{
	                        maxval = -10000000;
	                        maxvalind = 1;
	                        for (i = 1; i <= hmm.getN(); i++) 
	                        {
	                                val = delta[t-1][i] + (hmm.getA()[i][j]);
	                                if (val > maxval) 
	                                {
	                                        maxval = val;
	                                        maxvalind = i;
	                                }
	                        }	 
	                        delta[t][j] = maxval + biot[j][t]; 
	                        psi[t][j] = maxvalind;	 
	                }
	        }
			
			
	        /* 3. Termination */
	        
	        double pprob = -10000;
	        q[T] = 1;
	        for (i = 1; i <= hmm.getN(); i++) 
	        {
	                if (delta[T][i] > pprob) 
	                {
	                		pprob = delta[T][i];
	                        q[T] = i;
	                }
	        }
	 
	 
		/* 4. Path (state sequence) backtracking */

		for (t = T - 1; t >= 1; t--)
		{
			q[t] = psi[t+1][q[t+1]];
		}    
	        
		return pprob;
		
	}
	
//	void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,
//			double **gamma, int *pniter, 
//			double *plogprobinit, double *plogprobfinal)
	
	public static String BaumWelch(HMM hmm, int T, int[] O, double[][] alpha, double[][] beta,double[][] gamma)
	{
		int	i, j, k;
		int	t, l = 0;

		//  output 
		int pniter;
		double plogprobinit;
		double plogprobfinal;
		
		double	logprobf, logprobb,  threshold;
		double	numeratorA, denominatorA;
		double	numeratorB, denominatorB;	
		double[] scale;
		double delta, deltaprev, logprobprev;
		deltaprev = 10e-70;		
		double[][][] xi = new double[T+1][hmm.getN()+1][hmm.getN()+1] ;
		scale = new double[T+1]; 				
//		ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
//		*plogprobinit = logprobf; /* log P(O |intial model) */
//		BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
//		ComputeGamma(phmm, T, alpha, beta, gamma);
//		ComputeXi(phmm, T, O, alpha, beta, xi);
//		logprobprev = logprobf;		
		logprobf = ForwardWithScale(hmm,T,O,alpha,scale);  // log P(O |intial model)
		plogprobinit = logprobf;
		
		//System.out.println("plogprobinit:" + plogprobinit);		
		BackwardWithScale(hmm, T, O, beta, scale);		
		ComputeGamma(hmm, T, alpha, beta, gamma);
		ComputeXi(hmm, T, O, alpha, beta, xi);		
		logprobprev = logprobf;		
	    do{	
			/* reestimate frequency of state i in time t=1 */
			// for (i = 1; i <= phmm->N; i++) 
			for (i = 1; i <= hmm.getN(); i++) 
			{
				//phmm->pi[i] = .001 + .999*gamma[1][i];
				hmm.setPI_i(i,0.001 + 0.999*gamma[1][i]);
				//System.out.println(gamma[1][i]);
			}

			/* reestimate transition matrix  and symbol prob in
			   each state */
			// for (i = 1; i <= phmm->N; i++) { 
			for (i = 1; i <= hmm.getN(); i++) 
			{ 
				denominatorA = 0.0;
				for (t = 1; t <= (T - 1); t++) 
					denominatorA += gamma[t][i];

			//	for (j = 1; j <= phmm->N; j++) {
				for (j = 1; j <= hmm.getN(); j++) {
					numeratorA = 0.0;
					for (t = 1; t <= (T - 1); t++) 
						numeratorA += xi[t][i][j];
						// phmm->A[i][j] = .001 + .999*numeratorA/denominatorA;
						   hmm.setA_i_j(i, j, 0.001 + 0.999*numeratorA/denominatorA);
				}

				denominatorB = denominatorA + gamma[T][i]; 				
				
				//for (k = 1; k <= phmm->M; k++) {
				for (k = 1; k <= hmm.getM(); k++) 
				{
					numeratorB = 0.0;
					for (t = 1; t <= T; t++) 
					{
						if ((O[t]==k)) 
							numeratorB += gamma[t][i];
					}
					// phmm->B[i][k] = .001 + .999*numeratorB/denominatorB;
					// shold note that k is not j
					hmm.setB_i_j(i, k, 0.001 + 0.999*numeratorB/denominatorB);
				}
			}

			//ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
			logprobf = ForwardWithScale(hmm, T, O, alpha, scale);

//			for(int mmm1=1;mmm1<=T;mmm1++)
//			{
//				for(int mmm2=1;mmm2<=hmm.getN();mmm2++)
//				{
//						System.out.print("alpha [" + mmm1 + "][" + mmm2 + "]"+ alpha[mmm1][mmm2] + "\t");
//				}
//				System.out.println();
//			}
			
			//BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
			BackwardWithScale(hmm, T, O, beta, scale);
			
		
			//ComputeGamma(phmm, T, alpha, beta, gamma);
			ComputeGamma(hmm, T, alpha, beta, gamma);
			
			
			//ComputeXi(phmm, T, O, alpha, beta, xi);
			ComputeXi(hmm, T, O, alpha, beta, xi);

			
			/* compute difference between log probability of 
			   two iterations */
			delta = logprobf - logprobprev; 
			logprobprev = logprobf;
			l++;		
			
//			System.out.println("delta logprobf logprobprev" + delta + " " + logprobf + " " +logprobprev);
	
		}while (delta > 0.001); /* if log probability does not 
	                                  change much, exit */ 
		pniter = l;
		plogprobfinal = logprobf; /* log P(O|estimated model) */
		return (pniter + " " + plogprobinit + " " + plogprobfinal ); // return iteration,initial Log Prob(observation | init model)
																	 // "Log Prob(observation | estimated model)
	}
	
	public static double Forward(HMM hmm, int T, int[] O, double[][] alpha)
	{
	        int     i, j;   /* state indices */
	        int     t;      /* time index */
	 
	        double sum;     /* partial sum */
	 
	        /* 1. Initialization */
	 
	        System.out.println();
	        //for (i = 1; i <= phmm->N; i++)
	        for (i = 1; i <= hmm.getN(); i++)
	               // alpha[1][i] = phmm->pi[i]* phmm->B[i][O[1]];
	        {
	        		alpha[1][i] = hmm.getPI_i(i)*(hmm.getB()[i][O[1]]);
	        		System.out.println(alpha[1][i]);
	        }
	        
	        /* 2. Induction */
	 
	        for (t = 1; t < T; t++) {
	                //for (j = 1; j <= phmm->N; j++) {
	        		for (j = 1; j <= hmm.getN(); j++) {
	                        sum = 0.0;
	                        //for (i = 1; i <= phmm->N; i++)
	                        for (i = 1; i <= hmm.getN(); i++)
	                               // sum += alpha[t][i]* (phmm->A[i][j]);
	                        		sum += alpha[t][i]* (hmm.getA()[i][j]);
	                        //alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
	                        alpha[t+1][j] = sum*(hmm.getB()[j][O[t+1]]);
	                }
	        }
	 
	        /* 3. Termination */
	        double pprob = 0.0;
	        for (i = 1; i <= hmm.getN(); i++)
	                pprob += alpha[T][i];
	        System.out.println(pprob);
	        return pprob;
	 
	}
	
	
	////////////////////  forward with scale  #################################3
	
	
	public static double ForwardWithScale(HMM hmm, int T, int[] O, double[][] alpha,double[] scale)
		/*  pprob is the LOG probability */
		{
			int	i, j; 	/* state indices */
			int	t;	/* time index */

			double sum;	/* partial sum */

			/* 1. Initialization */

			scale[1] = 0.0;	
			//for (i = 1; i <= phmm->N; i++) {
			for (i = 1; i <= hmm.getN(); i++) {
				//alpha[1][i] = phmm->pi[i]* (phmm->B[i][O[1]]);
				alpha[1][i] = hmm.getPI_i(i)* (hmm.getB()[i][O[1]]);
				scale[1] += alpha[1][i];
			}
			//for (i = 1; i <= phmm->N; i++) 
			for (i = 1; i <= hmm.getN(); i++) 
				alpha[1][i] /= scale[1]; 
			
			/* 2. Induction */

			for (t = 1; t <= T - 1; t++) {
				scale[t+1] = 0.0;
				//for (j = 1; j <= phmm->N; j++) {
				for (j = 1; j <= hmm.getN(); j++) {
					sum = 0.0;
					//for (i = 1; i <= phmm->N; i++) 
					for (i = 1; i <= hmm.getN(); i++) 
						//sum += alpha[t][i]* (phmm->A[i][j]); 
						sum += alpha[t][i]* (hmm.getA()[i][j]);

					//alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);
					alpha[t+1][j] = sum*(hmm.getB()[j][O[t+1]]);
					scale[t+1] += alpha[t+1][j];
				}
				//for (j = 1; j <= phmm->N; j++) 
				for (j = 1; j <= hmm.getN(); j++) 
					alpha[t+1][j] /= scale[t+1]; 
			}

			/* 3. Termination */
			double pprob = 0.0;

			for (t = 1; t <= T; t++)
				pprob += Math.log(scale[t]);
			
			return pprob;
			
		}
	
	///////////////////##  end of forward with scale ###########################
	
	
	public double Backward(HMM hmm, int T, int[] O, double[][] beta)
	{
	        int     i, j;   /* state indices */
	        int     t;      /* time index */
	        double sum; 
	        /* 1. Initialization */	 
	       // for (i = 1; i <= phmm->N; i++)
	        for (i = 1; i <= hmm.getN(); i++)
	                beta[T][i] = 1.0;
	 
	        /* 2. Induction */	 
	        for (t = T - 1; t >= 1; t--) {
	                //for (i = 1; i <= phmm->N; i++) {
	        			for (i = 1; i <= hmm.getN(); i++) {
	                        sum = 0.0;
	                        //for (j = 1; j <= phmm->N; j++)
	                        for (j = 1; j <= hmm.getN(); j++)
	                        //sum += phmm->A[i][j] *(phmm->B[j][O[t+1]])*beta[t+1][j];
	                        sum += hmm.getA()[i][j] *(hmm.getB()[j][O[t+1]])*beta[t+1][j];
	                        beta[t][i] = sum;
	 
	                }
	        }
	 
	        /* 3. Termination */
	        double pprob = 0.0;
	        //for (i = 1; i <= phmm->N; i++)
	        for (i = 1; i <= hmm.getN(); i++)
	                pprob += beta[1][i];
	        return pprob;
	 
	}
	
	public static void BackwardWithScale(HMM hmm, int T, int[] O, double[][] beta,double[] scale)
	{
		        int     i, j;   /* state indices */
		        int     t;      /* time index */
			
		        double sum;
		 
		 
		        /* 1. Initialization */
		 
		       // for (i = 1; i <= phmm->N; i++)
		        for (i = 1; i <= hmm.getN(); i++)
		                beta[T][i] = 1.0/scale[T]; 
		 
		        /* 2. Induction */
		 
		        for (t = T - 1; t >= 1; t--) 
		        {
		              //  for (i = 1; i <= phmm->N; i++) {
		        	 for (i = 1; i <= hmm.getN(); i++) 
		        	 {
		        		 		sum = 0.0;
		                        //for (j = 1; j <= phmm->N; j++)
		        		 		for (j = 1; j <= hmm.getN(); j++)
		                        	//sum += phmm->A[i][j] *(phmm->B[j][O[t+1]])*beta[t+1][j];
		        		 			sum += hmm.getA()[i][j] *(hmm.getB()[j][O[t+1]])*beta[t+1][j];
		                        beta[t][i] = sum/scale[t];
		 
		             }
		        }		        
	}
	
	public static void ComputeGamma(HMM hmm, int T, double[][] alpha, double[][] beta,double[][] gamma)
	{
			int 	i, j;
			int	t;
			double	denominator;

			for (t = 1; t <= T; t++) 
			{
				denominator = 0.0;
				//for (j = 1; j <= phmm->N; j++) {
				for (j = 1; j <= hmm.getN(); j++) {
					gamma[t][j] = alpha[t][j]*beta[t][j];
					denominator += gamma[t][j];
				}

				//for (i = 1; i <= phmm->N; i++) 
				for (i = 1; i <= hmm.getN(); i++) 
					gamma[t][i] = gamma[t][i]/denominator;
			}
	}
	
	
	public static void ComputeXi(HMM hmm, int T, int[] O, double[][] alpha, double[][] beta,double[][][] xi)
	{
			int i, j;
			int t;
			double sum;

			for (t = 1; t <= T - 1; t++) {
				sum = 0.0;	
				//for (i = 1; i <= phmm->N; i++) 
				for (i = 1; i <= hmm.getN(); i++)
					//for (j = 1; j <= phmm->N; j++) {
					for (j = 1; j <= hmm.getN(); j++) {
						//xi[t][i][j] = alpha[t][i]*beta[t+1][j]*(phmm->A[i][j])*(phmm->B[j][O[t+1]]);
						xi[t][i][j] = alpha[t][i]*beta[t+1][j]*(hmm.getA()[i][j])*(hmm.getB()[j][O[t+1]]);
						sum += xi[t][i][j];
					}

				//for (i = 1; i <= phmm->N; i++) 
				for (i = 1; i <= hmm.getN(); i++) 
					//for (j = 1; j <= phmm->N; j++) 
					for (j = 1; j <= hmm.getN(); j++) 
						xi[t][i][j]  /= sum;
			}
	}

}

class HMM{
	/* A[1..N][1..N]. a[i][j] is the transition prob
	   of going from state i at time t to state j
	   at time t+1.
	
	   B[1..N][1..M]. b[j][k] is the probability of
	   of observing symbol k in state j. 
	   
	   pi[1..N] pi[i] is the initial state distribution.
	   */	
	double[][] A;
	double[][] B; 
	double[]   PI;
	int M = -1;
	int N = -1;
	int A_count=0;
	int B_count=0;
	int k=-1;
	String[] wds = null;

	public HMM(String model){
		BufferedReader br;
		String line;		
		M=N=-1;
		int s_A = 0;
		int s_B = 0;
		int s_Pi= 0;
		try {
			br = new BufferedReader(new FileReader(model));
			while ((line = br.readLine()) != null) {												
				if(line.startsWith("M=")){					
					wds = line.split("\\s+");
					M = Integer.valueOf(wds[1]).intValue();
				}else if(line.startsWith("N=")){
					wds = line.split("\\s+");
					N = Integer.valueOf(wds[1]).intValue();
					A = new double[N+1][N+1];
					B = new double[N+1][M+1];	
					PI= new double[N+1];
				}			
				
				if((1==s_A)&&(!line.startsWith("B:"))){
					A_count++;
					wds = line.split("\\s+");
					for(k=1;k<=wds.length;k++)
					{
						A[A_count][k] = Double.valueOf(wds[k-1]).doubleValue();
					}
				}else if((1==s_B)&&(!line.toUpperCase().startsWith("PI"))){
					B_count++;
					wds = line.split("\\s+");
					for(k=1;k<=wds.length;k++)
					{
						B[B_count][k] = Double.valueOf(wds[k-1]).doubleValue();
					}
				}else if(1==s_Pi){
					wds = line.split("\\s+");
					for(k=1;k<=wds.length;k++)
					{
						PI[k] = Double.valueOf(wds[k-1]).doubleValue();
					}
				}			
				
				if(line.startsWith("A:")){
					s_A = 1;
					s_B = 0;
					s_Pi= 0;
				}else if(line.startsWith("B:")){
					s_A = 0;
					s_B = 1;
					s_Pi= 0;
				}else if(line.toUpperCase().startsWith("PI")){
					s_A = 0;
					s_B = 0;
					s_Pi= 1;
				}	
				
//				else if(line.startsWith("A"))
//				{
//					A_count++;
//					wds = line.split("\\s+");
//					for(k=1;k<wds.length;k++)
//					{
//						A[A_count][k] = Double.valueOf(wds[k]).doubleValue();
//					}
//				}
//				else if(line.startsWith("B"))
//				{
//					B_count++;
//					wds = line.split("\\s+");
//					for(k=1;k<wds.length;k++)
//					{
//						B[B_count][k] = Double.valueOf(wds[k]).doubleValue();
//					}
//				}
//				else if(line.toUpperCase().startsWith("PI"))
//				{
//					wds = line.split("\\s+");
//					for(k=1;k<wds.length;k++)
//					{
//						PI[k] = Double.valueOf(wds[k]).doubleValue();
//					}
//				}
				
				
			}
		} catch (FileNotFoundException e) {			
			e.printStackTrace();
		} catch (IOException e) {			
			e.printStackTrace();
		}			
		
	}
	
	public void Output_HMMMode() {
		// TODO Auto-generated method stub
		
	}

	public HMM() 
	{
		// TODO Auto-generated constructor stub
	}

	public void CreatABPIMatrix(int M,int N)
	{
		A = new double[N+1][N+1];
		B = new double[N+1][M+1];
		PI = new double[N+1];
	}
	
	public void Output_HMMModel() throws Exception
	{
		FileWriter fw = new FileWriter(new File("C:\\MHMM\\MHMM\\MHMM\\16vpA\\SS2.hmm"));
		int i,j;
		System.out.println("M= " + M);
		System.out.println("N= " + N);
		fw.append("M= " + M + "\n");
		fw.append("N= " + N + "\n");
		System.out.println("A: ");
		fw.append("A: \n");
		for(i=1;i<=N;i++)
		{
			for(j=1;j<=N;j++)
			{
				System.out.print(A[i][j] + "\t");
				fw.append(A[i][j] + "\t");
			}
			System.out.println();
			fw.append("\n");
		}
		System.out.println("B: ");
		fw.append("B: \n");
		for(i=1;i<=N;i++)
		{
			for(j=1;j<=M;j++)
			{
				System.out.print(B[i][j] + "\t");
				fw.append(B[i][j] + "\t");
			}
			System.out.println();
			fw.append("\n");
		}
		System.out.println("pi:");
		fw.append("pi:\n");
		for(i=1;i<=N;i++)
		{
			System.out.print(PI[i] + "\t");
			fw.append(PI[i] + "\t");
		}
		fw.flush();
	}

	public int getM(){
		return M;
	}


	public void setM(int m) 
	{
		M = m;
	}


	public int getN() 
	{
		return N;
	}


	public void setN(int n) 
	{
		N = n;
	}


	public double[][] getA() 
	{
		return A;
	}

	public double getA_i_j(int i,int j)
	{
		return A[i][j];
	}

	public void setA_i_j(int i,int j,double new_value)
	{
		A[i][j] = new_value;
	}
	
	public void setB_i_j(int i,int j,double new_value)
	{
		B[i][j] = new_value;
	}
	
	public void setA(double[][] a) 
	{
		A = a;
	}


	public double[][] getB() 
	{
		return B;
	}


	public void setB(double[][] b) 
	{
		B = b;
	}


	public double[] getPI() 
	{
		return PI;
	}


	public void setPI(double[] pi) 
	{
		PI = pi;
	}
	
	public void setPI_i(int i,double pi_value)
	{
		PI[i] = pi_value;
	}
	
	public double getPI_i(int i)
	{
		return PI[i];
	}

	public void cloneA(double[][] a2)
	{
		int i,j;
		for(i=1;i<=N;i++)
		{
			for(j=1;j<=N;j++)
			{
				A[i][j] = a2[i][j];
			}
			System.out.println();
		}		
	}
	
	public void cloneB(double[][] a2)
	{
		int i,j;
		for(i=1;i<=N;i++)
		{
			for(j=1;j<=M;j++)
			{
				B[i][j] = a2[i][j];
			}
			System.out.println();
		}
	}

	public void clonePI(double[] pi2) 
	{
		int i;
		for(i=1;i<=N;i++)
		{
			PI[i]=pi2[i];
		}		
	}
	
}

class TOOLS
{
	void CLONE(HMM hmm1,HMM hmm2)
	{
		 int i, j, k;
		 
	     // phmm2->M = phmm1->M;
		    hmm2.setM(hmm1.getM());
	 
	     // phmm2->N = phmm1->N;
		    hmm2.setN(hmm1.getN());	 
		    
	     // phmm2->A = (double **) dmatrix(1, phmm2->N, 1, phmm2->N);		    
	 
//	        for (i = 1; i <= phmm2->N; i++)
//	                for (j = 1; j <= phmm2->N; j++)
//	                        phmm2->A[i][j] = phmm1->A[i][j];
//	 
//	        phmm2->B = (double **) dmatrix(1, phmm2->N, 1, phmm2->M);
		    
		    hmm2.CreatABPIMatrix(hmm2.getM(), hmm2.getN());
		    
//	        for (j = 1; j <= phmm2->N; j++)
//	                for (k = 1; k <= phmm2->M; k++)
//	                        phmm2->B[j][k] = phmm1->B[j][k];
	 
		    hmm2.cloneA(hmm1.getA());
		    hmm2.cloneB(hmm1.getB());
		    hmm2.clonePI(hmm1.getPI());	       
	}
}