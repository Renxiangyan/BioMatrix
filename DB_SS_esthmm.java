package machinelearning;

import java.text.DecimalFormat;

public class DB_SS_esthmm {
	// Estimates the HMM from a given symbol sequence using BaumWelch
	// ./esthmm test.hmm test.seq
	// output: HMM model
	// Problem 3 How do we adjust the model parameter to maximize the P(O|model) 
	// BaumWelch method	
	public static void main(String[] args) throws Exception {
//		if(args.length!=2)
//		{
//			System.out.println("Usage error");
//			System.out.println("Usage: -jar I-HMM_esthmm.jar <model.hmm> <obs.seq> \n");
//			System.exit(0);
//		}		
//		System.exit(1);				
		String mpth = "C:\\MHMM\\MHMM\\MHMM\\16vpA\\";
		HMM hmm = new HMM("C:\\MHMM\\MHMM\\MHMM\\16vpA\\SS2.hmm");	
//		HMM hmm = new HMM(args[0]);
		
//		hmm.Output_HMMMode();				
		int[] obs = MHMM.ReadSequence_fasta(mpth + "seq.txt");
//		int[] obs = ReadSequence(args[1]);
//		hmm.Output_HMMModel();
//		System.out.println(obs.length);
		
		int[] q = new int[obs.length];
		double[][] delta = new double[obs.length+1][hmm.getN()+1];
		int[][] psi = new int[obs.length+1][hmm.getN()+1];			
		double[][] alpha = new double[obs.length][hmm.getN()+1];
		double[] scale = new double[obs.length];		
		
		// Question 3		
		// gamma = dmatrix(1, T, 1, hmm.N);
		
		//String std = "CCCCCCCCCCCHHHHHHHHHHHHCCCHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHCHHHHHCEEEEECCHHHHHCCCCCCCCCCCEEEECCEECCHHHHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCHHHHHHHHHHHHHHCCCCCCCEEEEC";

		String std = MHMM.read_fasta(mpth + "16vpA.sec");
		std = std.replaceAll("\\*", "");
		double[][] gamma = new double[obs.length][hmm.getN()+1];
		double[][] beta = new double[obs.length][hmm.getN()+1];		
		String[] pred = MHMM.BaumWelch(hmm,obs.length-1,obs,alpha, beta,gamma).split("\\s+");				
		hmm.Output_HMMModel();							
		///////////////////////////////////////////////////////		
		double prob = MHMM.Viterbi(hmm,obs.length-1,obs,delta,psi,q);		
//		DecimalFormat df = new DecimalFormat("###.###");      // Correct the identity to 3 decimal places. 
//		System.out.println("Viterbi  MLE log prob " + df.format(Math.log(prob)));		
		String ss = MHMM.Output_status(q);
		double s = judge(ss,std);
		System.out.println("Q3=" + s);
//		System.out.println();
//		prob=MHMM.ViterbiLog(hmm,obs.length-1,obs,delta,psi,q);		
//		MHMM.Output_status(q);
		//System.exit(1);
//		System.out.println();
//		System.out.println("Viterbi  MLE log prob " + df.format(prob));
		
		
		///////////////////////////////////////////////////////	
		//System.exit(0);		
		///// final log score divid the length of the sequence would be better
//		System.out.println("\nIteration=" + pred[0] + "initial Log" + pred[1] + " final Log" + pred[2]);
			
// Question 2
		prob = MHMM.Viterbi(hmm,obs.length-1,obs,delta,psi,q);		
		DecimalFormat df = new DecimalFormat("###.###");      // Correct the identity to 3 decimal places. 
		System.out.println("Viterbi  MLE log prob " + df.format(Math.log(prob)));		
		MHMM.Output_status(q);
		System.out.println();
		prob=MHMM.ViterbiLog(hmm,obs.length-1,obs,delta,psi,q);		
		MHMM.Output_status(q);
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
	
	public static double judge(String s1, String s2){
		if(s1.length()!=s2.length()){
			System.out.println("Length Error ");
		}
		double s = 0;
		for(int i=0;i<s1.length();i++){
			char a = s1.charAt(i);
			char b = s2.charAt(i);
			if(a==b){
				s++;
			}
		}
		return s/s1.length();
	}

}



















